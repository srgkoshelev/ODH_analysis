import math
import heat_transfer as ht
from copy import copy

# Setting up the units
ureg = ht.ureg
Q_ = ureg.Quantity

# Loading FESHM 4240 Failure rates
from .FESHM4240_TABLES import TABLE_1, TABLE_2

# Probability of failure on demand for main cases
lambda_odh = 2.3e-6/ureg.hr  # from CMTF Hi Bay ODH EN01878 pp. 27-28. That is the most correct value I have seen in use
# TODO Update to value from J. Anderson's document
PFD_ODH = Q_('2 * 10^-3')
PFD_SOLENOID = TABLE_2['Valve, solenoid']['Failure to operate']
PFD_POWER = TABLE_1['Electrical Power Failure']['Demand rate']
PFD_DEWAR_INSULATION = TABLE_1['Dewar']['Loss of vacuum']
LAMBDA_FAN = TABLE_2['Fan']['Failure to run']
TRANSFER_LINE_LEAK_AREA = Q_('10 mm^2')
SHOW_SENS = 1e-7/ureg.hr


class Source:
    """Define the possible source of inert gas"""
    def __init__(self, name, Fluid, volume, N=1):
        self.name = name
        self.Fluid = Fluid
        self.leaks = {}
        # Number of sources if multiple exist, e.g. gas cylinders
        # Increases probability of failure by N.
        self.N = N
        # Calculating volume at standard conditions
        TempState = ht.ThermState(Fluid.name, backend=Fluid.backend)
        TempState.update('T', ht.T_NTP, 'P', ht.P_NTP)
        self.volume = volume*Fluid.Dmass/TempState.Dmass
        self.volume.ito(ureg.feet**3)
        # By default assume there is no isolation valve
        # that is used by ODH system
        self.isol_valve = False

    def gas_pipe_failure(self, Pipe, Fluid=None, N_welds=1, max_flow=None):
        """
        Calculate failure rate, flow rate and expected time duration of the
        event for gas pipe failure. Based on FESHM 4240.
        """
        # If Fluid not defined use Fluid of the Source
        Fluid = Fluid or self.Fluid
        # Failure rate coefficients; Piping failure rate is per unit of length,
        # weld is dependent on number of welds, pipe OD and wall thickness
        failure_rate_coeff = {'Piping': Pipe.L, 'Pipe weld': N_welds *
                              Pipe.OD / Pipe.wall}
        # Piping and weld leaks as per Table 2
        for cause in ['Piping', 'Pipe weld']:
            for mode in TABLE_2[cause].keys():
                if Pipe.D > 2 or mode != 'Large leak':  # Large leak only for D > 2"
                    name = f'{cause} {mode.lower()}: {Pipe}'
                    TempPipe = copy(Pipe)
                    # Average path for the flow will be half of piping length
                    # for gas piping
                    TempPipe.L = Pipe.L / 2
                    if mode == 'Rupture':
                        failure_rate = failure_rate_coeff[cause] * \
                            TABLE_2[cause][mode]
                        # For rupture calculate flow through available
                        # pipe area
                        area = Pipe.area
                    else:
                        failure_rate = failure_rate_coeff[cause] * \
                            TABLE_2[cause][mode]['Failure rate']
                        area = TABLE_2[cause][mode]['Area']
                    if max_flow is not None:
                        q_std_max = ht.piping.to_standard_flow(max_flow, Fluid)
                        q_std = min(self._leak_flow(TempPipe, area, Fluid),
                                    q_std_max)
                    else:
                        q_std = self._leak_flow(TempPipe, area, Fluid)
                    tau = self.volume/q_std
                    self.leaks[name] = (failure_rate.to(1/ureg.hr), q_std,
                                        tau.to(ureg.min))

    def transfer_line_failure(self, Pipe, Fluid=None, N_lines=1):
        """
        Calculate failure rate, flow rate and expected time duration of the event for transfer line failure. Based on FESHM 4240.
        """
        area_cases = {'Leak': TRANSFER_LINE_LEAK_AREA,  # Leak area, assumed based on piping leak areas in Table 2
                      'Rupture': Pipe.area,  # Full flow through the transfer line
                }
        for mode in TABLE_1['Fluid line']:
            name = f'Fluid line {mode.lower()}: {Pipe}'
            failure_rate = N_lines * TABLE_1['Fluid line'][mode]
            area = area_cases[mode]
            Fluid = Fluid or self.Fluid  # If Fluid not defined use Fluid of the Source
            q_std = self._leak_flow(Pipe, area, Fluid)
            tau = self.volume/q_std
            self.leaks[name] = (failure_rate.to(1/ureg.hr), q_std, tau.to(ureg.min))

    def dewar_insulation_failure(self, flow_rate, Fluid=None):
        """Calculate failure rate, flow rate and expected time duration of the
        failure event for the dewar insulation failure.

        Based on FESHM4240."""
        failure_rate = PFD_DEWAR_INSULATION
        if Fluid is None:
            Fluid = self.Fluid
        tau = self.volume/flow_rate
        self.leaks['Dewar insulation failure'] = (failure_rate.to(1/ureg.hr),
                                                  flow_rate, tau.to(ureg.min))

    def failure_mode(self, name, failure_rate, flow_rate, Fluid=None):
        """General failure mode."""
        Fluid = Fluid or self.Fluid  # If Fluid not defined use Fluid of the Source
        flow_rate = ht.piping.to_standard_flow(flow_rate, Fluid)
        tau = self.volume/flow_rate
        self.leaks[name] = (failure_rate.to(1/ureg.hr), flow_rate,
                            tau.to(ureg.min))

    def constant_leak(self, name, flow_rate):
        tau = self.volume/flow_rate
        self.leaks[name] = (None, flow_rate, tau.to(ureg.min))

    def _leak_flow(self, Pipe, Area, Fluid):
        d = (4*Area/math.pi)**0.5  # diameter for the leak opening
        Entrance = ht.piping.Entrance(d)
        Exit = ht.piping.Exit(d)
        TempPiping = ht.piping.Piping(Fluid)
        TempPiping.add(Entrance,
                       Pipe,
                       Exit,)
        m_dot = TempPiping.m_dot(ht.P_NTP)
        return ht.piping.to_standard_flow(m_dot, Fluid)

    @property
    def sol_PFD(self):
        """Calculate probability of failure on demand (PFD) for solenoid valve.

        If the source doesn't have isolating solenoid valve this probability is 1.
        """
        return (not self.isol_valve) or PFD_SOLENOID

    def __add__(self, other):
        if self.Fluid.name == other.Fluid.name:
            total_volume = self.volume + other.volume
            return Source(f'self.name+other.name', self.Fluid, total_volume)
        else:
            print('\nBoth volumes should contain the same fluid')
            return None

#    def __mul__(self, num):
#        Copy = copy(self)
#        Copy.N *= num
#        return Copy
#
#    def __rmul__(self, num):
#        Copy = copy(self)
#        Copy.N *= num
#        return Copy

    def __str__(self):
        return f'{self.name}, ODH source with ' + \
            f'{self.volume.to(ureg.ft**3):.3g~} ' + \
            f'of {self.Fluid.name} gas.'

    def print(self):
        print('{} is an ODH source of {} with volume of {:.3~}'.format(
            self.name, self.Fluid.name, self.volume))

    def print_leaks(self):
        for key in sorted(self.leaks.keys()):
            print('Failure mode: '+key)
            print('Failure rate: {:.2~}'.format(self.leaks[key][0]))
            print('Flow rate: {:.2~}'.format(
                self.leaks[key][1].to(ureg.ft**3/ureg.min)))
            print('Event duration: {:.2~}'.format(self.leaks[key][2]))
            print()


class Volume:
    """Volume/building affected by inert gases.
    """
    def __init__(self, name, volume, Q_fan, N_fans, Test_period):
        """
        Main method of protection against ODH is a complex system involving ODH heads and chassis or PLCs that power louvers and fans (which may fail separately). PFD_ODH describes the probability of failure of this system to register, transmit and respond to the reduction of oxygen concetration of the volume. Default value is based on analysis performed by J. Anderson and presented ... When necessary, the value can be redefined. One must take care calculating PFD_ODH as it may be complicated to properly add probabilities.
"""
        self.name = name
        self.volume = volume
        self.PFD_ODH = PFD_ODH  #  Default value for ODH system failure
        self.lambda_fan = LAMBDA_FAN  #  Table 2 value is default
        self.Q_fan = Q_fan
        self.N_fans = N_fans
        self.Test_period = Test_period

    def odh(self, sources, power_outage=False):
        """
        Calculate ODH fatality rate and recommend ODH class designation.
        Uses list of odh.source instances.
        Power specifies whether there is a power outage. Default is no outage.
        """
        self.phi = 0  # Recalculate fatality rate
        self.failure_modes = []
        # Probability of power failure in the building:
        # PFD_power if no outage, 1 if there is outage
        PFD_power_build = power_outage or PFD_POWER
        # Calculate fan probability of failure
        self._fan_fail(self.Test_period, self.Q_fan, self.N_fans)
        # Calculate fatality rates for each source
        for source in sources:
            for failure_mode_name, leak in source.leaks.items():
                leak_failure_rate = leak[0]
                if leak_failure_rate is not None:  # None for constant leak
                    self._fatality_no_response(source, failure_mode_name, leak,
                                               source.sol_PFD, PFD_power_build)
                    self._fatality_fan_powered(source, failure_mode_name, leak,
                                               PFD_power_build)
                else:
                    # TODO rework constant leak (will throw errors)
                    raise NotImplementedError('Constant leak analysis is not '
                                              'implemented yet.')
                    # O2_conc = conc_vent(self.volume, q_leak, 0*ureg('ft^3/min'), tau)
                    # F_i = self._fatality_prob(O2_conc)
                    # if F_i > 0:
                        # print('Constant leak can be fatal: {source}, \
                        # Leak: {leak_rate:.3g~}, tau: {tau:.3g~}')
                    # phi_i = F_i / tau  #  Assessing fatality rate of the cont. leak
                    # self.failure_modes.append((phi_i, source, failure_mode_name, O2_conc,
                                               # 1.0, 1.0,
                                               # F_i, power_outage, q_leak, tau, Q_fan))
        # Sort failure modes by fatality rate
        self.failure_modes.sort(key=lambda x: x[0], reverse=True)

    def _fatality_no_response(self, source, failure_mode_name, leak, sol_PFD,
                              PFD_power_build):
        """Calculate fatality rate in the volume for ODH protection failure.

        Calculate failure rate of leak occuring and no safety response
        occuring due to power failure and isolation solenoid failure,
        or power on and ODH system failure.
        O2 concentration is limited only by amount of inert gas the source has.
        Fans are not operational.

        Parameters
        ----------
        leak : element Source.leaks.values
        """
        (leak_failure_rate, q_leak, tau) = leak
        P_no_response = float(PFD_power_build)*float(sol_PFD) +\
            (1-PFD_power_build)*self.PFD_ODH
        P_i = leak_failure_rate * P_no_response
        Q_fan = 0 * ureg('ft^3/min')
        O2_conc = conc_vent(self.volume, q_leak, Q_fan, tau)
        F_i = self._fatality_prob(O2_conc)
        phi_i = P_i*F_i
        self.phi += phi_i
        self.failure_modes.append((phi_i, source, failure_mode_name, O2_conc,
                                   leak_failure_rate, P_i, F_i,
                                   PFD_power_build == 1, q_leak, tau, Q_fan))

    def _fatality_fan_powered(self, source, failure_mode_name, leak,
                              PFD_power_build):
        """Calculate fatality rates for fan failure on demand.

        Calculate fatality rates for the case of ODH system responding and
        fans powered.
        """
        (leak_failure_rate, q_leak, tau) = leak
        for (P_fan, Q_fan) in self.Fan_flowrates:
            # Probability of power on, ODH system working, and m number of fans
            # with flow rate Q_fan on.
            P_response = (1-PFD_power_build) * (1-self.PFD_ODH) * P_fan
            P_i = leak_failure_rate * P_response
            O2_conc = conc_vent(self.volume, q_leak, Q_fan, tau)
            F_i = self._fatality_prob(O2_conc)
            phi_i = P_i*F_i
            self.phi += phi_i
            self.failure_modes.append((phi_i, source, failure_mode_name,
                                       O2_conc, leak_failure_rate, P_i, F_i,
                                       PFD_power_build == 1, q_leak, tau,
                                       Q_fan))

    def _fan_fail(self, Test_period, Q_fan, N_fans):
        """
        Calculate (Probability, flow) pairs for all combinations of fans
        working. All fans are expected to have same volume flow.
        """
        # TODO add fans with different volumetric rates (see report as well)
        Fail_rate = self.lambda_fan
        Fan_flowrates = []
        for m in range(N_fans+1):
            # Probability of exactly m units starting
            P_m_fan_work = prob_m_of_n(m, N_fans, Test_period, Fail_rate)
            Fan_flowrates.append((P_m_fan_work, Q_fan*m))
        self.Fan_flowrates = Fan_flowrates

    def _fatality_prob(self, O2_conc):
        if O2_conc >= 0.18:  # Lowest oxygen concentration above 18%
            Fi = 0
        elif O2_conc <= 0.088:  # 8.8% of oxygen is assumed to be 100% fatal
            Fi = 1
        else:
            # Fi formula, reverse engineered using 8.8% and 18% thresholds
            Fi = 10**(6.5-76*O2_conc)
        return Fi

    def odh_class(self):
        if self.phi < 1e-7/ureg.hr:
            return 0
        elif self.phi < 1e-5/ureg.hr:
            return 1
        elif self.phi < 1e-3/ureg.hr:
            return 2
        else:
            print('ODH fatality rate is too high. Please, check calculations')
            return None

    def report(self, brief=True):
        title = f'ODH report for {self}'
        padding = len(title) + 10
        print('#'*padding)
        print(title)
        print('-'*padding)
        if brief:
            print('Printing brief ODH report')
            print(f'Only leaks with Fatality rate > {SHOW_SENS} are shown')
        for failure_mode in self.failure_modes:
            phi_i = failure_mode[0]
            source = failure_mode[1]
            failure_mode_name = failure_mode[2]
            O2_conc = failure_mode[3]
            leak_failure_rate = failure_mode[4]
            P_i = failure_mode[5]
            F_i = failure_mode[6]
            power_outage = failure_mode[7]
            q_leak = failure_mode[8]
            tau = failure_mode[9]
            N_fan = int(failure_mode[10] / self.Q_fan)
            Q_fan = failure_mode[10]
            if phi_i >= SHOW_SENS or not brief:
                print(f'\n Source: {source.name}')
                print(f' Failure: {failure_mode_name}')
                print(f' Fatality rate: {phi_i.to(1/ureg.hr):.2~}')
                print(f' Building is powered: {not power_outage}')
                print(f' Oxygen concentration: {O2_conc:.0%}, '
                      f'{O2_conc/0.21:.0%} percent of norm')
                print(f' Leak failure rate: {leak_failure_rate:.3g~}')
                print(' ODH protection PFD: '
                      f'{(P_i/leak_failure_rate).to(ureg.dimensionless):.2~}')
                print(f' Total failure rate: {P_i.to(1/ureg.hr):.2~}')
                print(f' Leak rate: {q_leak:.2~}')
                print(f' Event duration: {tau:.2~}')
                print(f' Fans working: {N_fan:}')
                print(f' Fan rate: {Q_fan:.2~}')
                print(f' Fatality prob: {F_i:.2g}')

    def __str__(self):
        return (f'Volume: {self.name}, {self.volume.to(ureg.ft**3):.2~}')

    # def source_safe(self, source, escape = True):
    #    """
    #    Estimate the impact of the Source volume on oxygen concetration. Smaller sources might not be able to drop oxygen concentration to dangerous levels.
    #    """
    #    if escape == True:  # if mixed air is allowed to escape within considered volume
    #        O2_conc = 0.21*self.volume/(self.volume+source.volume)
    #    else:  # worst case; inert gas is trapped and expells the air outside the considered volume
    #        O2_conc = 0.21*(1-source.volume/self.volume)
    #    return self._fatality_prob(O2_conc) == 0


def prob_m_of_n(m, n, T, l):
    """
    Calculate the probability of m out of n units working.
    inputs:
        T - test period, hr
        l = lambda = failure rate for fan, 1/hr
    """
    C_n_m = math.factorial(n)/(math.factorial(n-m)*math.factorial(m))
    PFD_one_unit = l*T
    # Adjustment coefficient: T/(n-m+1) will be average failure reveal time
    # (D. Smith, Reliability..., p. 108)
    F_adj = 1/(n-m+1)
    m_of_n = C_n_m*(PFD_one_unit)**(n-m)*(1-PFD_one_unit)**m*F_adj
    return m_of_n


def conc_vent(V, R, Q, t):
    """
    V - volume of the confined space (ft3 or m3)
    R - spill rate into confined space (scfm or m3/s)
    Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside
    t = time, (minutes or seconds) beginning of release is at t=0
    C - oxygen concentration in confined space
    Case B
    """
    if Q > 0:
        C = 0.21/(Q+R) * (Q+R*math.e**-(Q+R)/V*t)
    elif abs(Q) <= R:
        C = 0.21*math.e**-(R/V*t)
    elif abs(Q) > R:
        C = 0.21*(1-R/abs(Q)*(1-math.e**-(abs(Q)*t/V)))
    return C


def conc_final(V, R, Q):
    """
    V - volume of the confined space (ft3 or m3)
    R - spill rate into confined space (scfm or m3/s)
    Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds to blowing air into the confined space, negative - drawing contaminated air outside
    C - oxygen concentration in confined space
    Case B
    """
    if Q > 0*ureg('ft^3/min'):
        C = 0.21/(Q+R)*Q
    elif abs(Q) <= abs(R):
        C = 0
    elif abs(Q) > abs(R):
        C = 0.21*(1-R/abs(Q))
    return C

def conc_after(V, C_e, Q, t, t_e):
    """
    V - volume of the confined space (ft3 or m3)
    R - spill rate into confined space (scfm or m3/s)
    Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds to blowing air into the confined space, negative - drawing contaminated air outside
    t = time, (minutes or seconds) beginning of release is at t=0
    C - oxygen concentration in confined space
    C_e = oxygen concentration when the release has ended
    """
    C = 0.21-(0.21-C_e)*math.e**-(abs(Q)/V*(t-t_e))
    return C


def print_result(*Volumes):
    """
    Print the results of the ODH analysis for a volume. If several volumes given (in case of interlapping volumes) the worst case will be printed.
    """
    max_phi = -1/ureg.hr
    for volume in Volumes:
        if volume.phi > max_phi:
            max_volume = volume
    line_1 = '# Fatality rate for {} is {:.1e}  # '.format(max_volume, volume.phi)
    pad = len(line_1)
    line_2 = '# Recommended ODH class {}'.format(max_volume.odh_class()).ljust(pad-1)+'#'
    print('#'*pad)
    print(line_1)
    print(line_2)
    print('#'*pad)
