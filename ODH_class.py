"""odh - Oxygen Deficiency Hazard analysis tool.

Based on Fermilab ES&H Manual chapter 4240 (see
https://eshq.fnal.gov/manuals/feshm/)

This package is maintained at https://github.com/srgkoshelev/ODH_analysis"""

import math
import heat_transfer as ht
from copy import copy

# Setting up the units
ureg = ht.ureg
Q_ = ureg.Quantity

# Loading FESHM 4240 Failure rates
from .FESHM4240_TABLES import TABLE_1, TABLE_2

# Probability of failure on demand for main cases
PFD_ODH = Q_('2 * 10^-3')
# TODO Update to value from J. Anderson's document
TRANSFER_LINE_LEAK_AREA = Q_('10 mm^2')
SHOW_SENS = 5e-8/ureg.hr


class Source:
    """Source of inert gas

    Attributes
    ----------
    sol_PFD : float
        Probability of failure on demand (PFD) for solenoid valve.
    """
    def __init__(self, name, fluid, volume, N=1):
        """Define the possible source of inert gas.

        Parameters
        ----------
        name : str
            Name of the source.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        volume : ureg.Quantity {length: 3}
            Volume of the fluid stored.
        N : int
            Quantity of the sources if several similar sources exist,
            e.g. gas bottles.
        """
        self.name = name
        self.fluid = fluid
        self.leaks = {}
        # Number of sources if multiple exist, e.g. gas cylinders
        # Increases probability of failure by N.
        self.N = N
        # Calculating volume at standard conditions
        temp_state = fluid.copy()
        temp_state.update('T', ht.T_NTP, 'P', ht.P_NTP)
        self.volume = volume*fluid.Dmass/temp_state.Dmass
        self.volume.ito(ureg.feet**3)
        # By default assume there is no isolation valve
        # that is used by ODH system
        self.isol_valve = False

    def gas_pipe_failure(self, Pipe, fluid=None, N_welds=1, max_flow=None):
        """Add gas pipe failure to the leaks dict.

        Store failure rate, flow rate and expected time duration of the
        event for gas pipe failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        Pipe : heat_transfer.piping.Pipe
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        N_welds : int
            Number of welds on the gas pipe.
        max_flow : ureg.Quantity {mass: 1, time: -1} or {length: 3, time: -1}
            Max mass or volumetric flow through if limited,
            e.g. by compressor output.
        """
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
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
                        q_std_max = ht.piping.to_standard_flow(max_flow, fluid)
                        q_std = min(self._leak_flow(TempPipe, area, fluid),
                                    q_std_max)
                    else:
                        q_std = self._leak_flow(TempPipe, area, fluid)
                    tau = self.volume/q_std
                    self.leaks[name] = (failure_rate.to(1/ureg.hr), q_std,
                                        tau.to(ureg.min))

    def transfer_line_failure(self, Pipe, fluid=None, N=1):
        """Add transfer line failure to leaks dict.

        Store failure rate, flow rate and expected time duration of
        the event for transfer line failure. Based on FESHM 4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        Pipe : heat_transfer.Pipe
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        N : int
            Number of bayonets/soft seals on the transfer line
        """
        # TODO Make leak and rupture areas adjustable, add info to docstring
        area_cases = {'Leak': TRANSFER_LINE_LEAK_AREA,
                      'Rupture': Pipe.area}
        for mode in TABLE_1['Fluid line']:
            name = f'Fluid line {mode.lower()}: {Pipe}'
            failure_rate = N * TABLE_1['Fluid line'][mode]
            area = area_cases[mode]
            # If fluid not defined use fluid of the Source
            fluid = fluid or self.fluid
            q_std = self._leak_flow(Pipe, area, fluid)
            tau = self.volume/q_std
            self.leaks[name] = (failure_rate.to(1/ureg.hr), q_std,
                                tau.to(ureg.min))

    def dewar_insulation_failure(self, flow_rate, fluid=None):
        """Add dewar insulation failure to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for the dewar insulation failure. Based on FESHM4240.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        flow_rate : ureg.Quantity {mass: 1, time: -1} or {length: 3, time: -1}
            Relief flow rate for the case of dewar insulation failure.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        """
        failure_rate = TABLE_1['Dewar']['Loss of vacuum']
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
        q_std = ht.piping.to_standard_flow(flow_rate, fluid)
        tau = self.volume/q_std
        self.leaks['Dewar insulation failure'] = (failure_rate.to(1/ureg.hr),
                                                  q_std, tau.to(ureg.min))

    def failure_mode(self, name, failure_rate, flow_rate, fluid=None, N=1):
        """Add general failure mode to leaks dict.

        Store failure rate, flow rate and expected time duration of the
        failure event for general failure mode.
        Failure modes are analyzed by `Volume.odh` method.

        Parameters
        ----------
        name : str
            Name of the failure mode
        failure rate : ureg.Quantity {time: -1}
            Failure rate of the failure mode,
            i.e. how often the failure occurs
        flow_rate : ureg.Quantity {mass: 1, time: -1} or {length: 3, time: -1}
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.
        N : int
            Quantity of similar failure modes.
        """
        # If fluid not defined use fluid of the Source
        fluid = fluid or self.fluid
        q_std = ht.piping.to_standard_flow(flow_rate, fluid)
        tau = self.volume/q_std
        self.leaks[name] = (N*failure_rate.to(1/ureg.hr), q_std,
                            tau.to(ureg.min))

    # def constant_leak(self, name, flow_rate):
    #     tau = self.volume/flow_rate
    #     self.leaks[name] = (None, flow_rate, tau.to(ureg.min))

    def _leak_flow(self, Pipe, Area, fluid):
        """Calculate leak flow through a pipe

        Parameters
        ----------
        Pipe : heat_transfer.Pipe
        Area : ureg.Quantity {length: 2}
            Area of the leak.
        fluid : heat_transfer.ThermState
            Thermodynamic state of the fluid stored in the source.

        Returns
        -------
        ureg.Quantity {length: 3, time: -1}
            Standard volumetric flow. Conditions are defined in
            heat_transfer package (generally NTP).
        """
        d = (4*Area/math.pi)**0.5  # diameter for the leak opening
        Entrance = ht.piping.Entrance(d)
        Exit = ht.piping.Exit(d)
        TempPiping = ht.piping.Piping(fluid)
        TempPiping.add(Entrance,
                       Pipe,
                       Exit,
        )
        m_dot = TempPiping.m_dot(ht.P_NTP)
        return ht.piping.to_standard_flow(m_dot, fluid)

    @property
    def sol_PFD(self):
        """Probability of failure of a solenoid device
            If the source doesn't have isolating solenoid valve
            the probability is 1.
        """
        return ((not self.isol_valve) or
                TABLE_2['Valve, solenoid']['Failure to operate'])

    @staticmethod
    def combine(name, sources):
        """Combine several ODH sources sharing volume.

        Can be used for failure modes affecting several sources in parallel.

        Parameters
        ----------
        name : str
            Name of the new combined source.
        sources : list of Source
            Sources connected together.

        Returns
        -------
        Source
            Combined source of inert gas.
        """
        fluid = ht.ThermState(sources[0].fluid.name, T=ht.T_NTP, P=ht.P_NTP)
        if all([source.fluid.name == fluid.name for source in sources]):
            total_volume = sum([source.volume for source in sources])
            return Source(name, fluid, total_volume)
        else:
            print('\nAll volumes should contain the same fluid')
            return None

    def __str__(self):
        return f'{self.name}, ODH source with ' + \
            f'{self.volume.to(ureg.ft**3):.3g~} ' + \
            f'of {self.fluid.name} gas.'

    def print_leaks(self):
        """Print information on the leaks defined for the source."""
        for key in sorted(self.leaks.keys()):
            print('Failure mode: '+key)
            print('Failure rate: {:.2~}'.format(self.leaks[key][0]))
            print('Flow rate: {:.2~}'.format(
                self.leaks[key][1].to(ureg.ft**3/ureg.min)))
            print('Event duration: {:.2~}'.format(self.leaks[key][2]))
            print()


class Volume:
    """Volume/building affected by inert gases."""
    def __init__(self, name, volume, Q_fan, N_fans, Test_period):
        """Define a volume affected by inert gas release from  a `Source`.

        Parameters
        ----------
        name : str
            Name of the volume.
        volume : ureg.Quantity {length: 3}
            Volume of the building or part of the building.
        Q_fan : ureg.Quantity {length: 3, time: -1}
            Volumetric flow of a single ODH fan installed in the volume.
        N_fans : int
            Number of fans installed.
        Test_period : ureg.Quantity {time: 1}
            Test period of the fans.
        """
        self.name = name
        self.volume = volume
        self.PFD_ODH = PFD_ODH  # Default value for ODH system failure
        self.lambda_fan = TABLE_2['Fan']['Failure to run']
        self.Q_fan = Q_fan
        self.N_fans = N_fans
        self.Test_period = Test_period

    def odh(self, sources, power_outage=False):
        """Calculate ODH fatality rate for given `Source`s.

        For each leak of each source ODH conditions are analyzed and
        fatality rates are calculated. The results are collected in
        failure_modes list.

        Parameters
        ----------
        sources : list
            Sources affecting the volume.
        power_outage : bool
            Shows whether there is a power outage is in effect.
            Default is no outage.
        """
        self.phi = 0  # Recalculate fatality rate
        self.failure_modes = []
        # Probability of power failure in the building:
        # PFD_power if no outage, 1 if there is outage
        PFD_power_build = (power_outage or
                           TABLE_1['Electrical Power Failure']['Demand rate'])
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
                    """O2_conc = conc_vent(self.volume, q_leak, 0*ureg('ft^3/min'), tau)
                    F_i = self._fatality_prob(O2_conc)
                    if F_i > 0:
                        print('Constant leak can be fatal: {source}, \
                        Leak: {leak_rate:.3g~}, tau: {tau:.3g~}')
                    phi_i = F_i / tau  #  Assessing fatality rate of the cont. leak
                    self.failure_modes.append((phi_i, source, failure_mode_name, O2_conc,
                                               1.0, 1.0,
                                               F_i, power_outage, q_leak, tau, Q_fan))
                    """
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
        Adds calculation results to the failure_modes list.

        Parameters
        ----------
        source : Source
        failure_mode_name : str
            Name of the failure mode
        leak : tuple (ureg.Quantity {time: -1},
                      ureg.Quantity {length: 3, time: -1},
                      ureg.Quantity {time: 1})
            Leak failure rate, volumetric flow rate, and time of the event.
        sol_PFD : float
            Probability of source solenoid failure.
        PFD_power_building : float
            Probability of power failure.
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
        fans powered but some of the fans failing on demand.
        See wiki for further explanation.
        Adds calculation results to the failure_modes list.

        Parameters
        ----------
        source : Source
        failure_mode_name : str
            Name of the failure mode
        leak : tuple (ureg.Quantity {time: -1},
                      ureg.Quantity {length: 3, time: -1},
                      ureg.Quantity {time: 1})
            Leak failure rate, volumetric flow rate, and time of the event.
        sol_PFD : float
            Probability of source solenoid failure.
        PFD_power_building : float
            Probability of power failure.
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
        """Calculate (Probability, flow) pairs for all combinations of fans
        working.

        All fans are expected to have same volume flow.

        Parameters
        ----------
        Test_period : ureg.Quantity {time: 1}
            Test period of the fans.
        Q_fan : ureg.Quantity {length: 3, time: -1}
            Volumetric flow of a single ODH fan installed in the volume.
        N_fans : int
            Number of fans installed.
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
        """Calculate fatality probability for given oxygen concentration.

        The equation is fitted from the FESHM 4240 plot.

        Parameters
        ----------
        O2_conc : float
            Oxygen concentration.

        Returns
        -------
        float
            Fatality rate.
        """
        if O2_conc >= 0.18:  # Lowest oxygen concentration above 18%
            Fi = 0
        elif O2_conc <= 0.088:  # 8.8% of oxygen is assumed to be 100% fatal
            Fi = 1
        else:
            # Fi formula, reverse engineered using 8.8% and 18% thresholds
            Fi = 10**(6.5-76*O2_conc)
        return Fi

    def odh_class(self):
        """Calculate ODH class as defined in FESHM 4240.

        Returns
        -------
        int
            ODH class.
        """
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
        """Print a report for failure modes and effects.

        The report is sorted by fatality rate descending."""
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
                print()
                print(f' Source:               {source.name}')
                print(f' Failure:              {failure_mode_name}')
                print(f' Fatality rate:        {phi_i.to(1/ureg.hr):.2~}')
                print(f' Building is powered:  {not power_outage}')
                print(f' Oxygen concentration: {O2_conc:.0%}, '
                      f'{O2_conc/0.21:.0%} percent of norm')
                print(f' Leak failure rate:    {leak_failure_rate:.3g~}')
                print(' ODH protection PFD:    '
                      f'{(P_i/leak_failure_rate).to(ureg.dimensionless):.2~}')
                print(f' Total failure rate:   {P_i.to(1/ureg.hr):.2~}')
                print(f' Leak rate:            {q_leak:.2~}')
                print(f' Event duration:       {tau:.2~}')
                print(f' Fans working:         {N_fan:}')
                print(f' Fan rate:             {Q_fan:.2~}')
                print(f' Fatality prob:        {F_i:.2g}')

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
    """Calculate the probability of m out of n units working.

    Calculation is done using binomial distribution.

    Parameters
    ----------
    m : int
    Number of units working.
    n : int
    Total number of units.
    T : ureg.Quantity {time: 1}
        Test period
    l : ureg.Quantity {time: -1}
        Failure rate (\\lambda) of a fan

    Returns
    -------
    float
        Probability of m out of n units working.
    """
    C_n_m = math.factorial(n)/(math.factorial(n-m)*math.factorial(m))
    PFD_one_unit = l*T
    # Adjustment coefficient: T/(n-m+1) will be average failure reveal time
    # (D. Smith, Reliability..., p. 108)
    F_adj = 1/(n-m+1)
    m_of_n = C_n_m*(PFD_one_unit)**(n-m)*(1-PFD_one_unit)**m*F_adj
    return m_of_n


def conc_vent(V, R, Q, t):
    """Calculate the oxygen concentration at the end of the event.

    As defined by FESHM 4240 6.1.A, Cases A, B, and C.

    Parameters
    ----------
    V : ureg.Quantity {length: 3}
        Volume of the confined space.
    R : ureg.Quantity {length: 3, time: -1}
        Volumetric spill rate into confined space.
    Q : ureg.Quantity {length: 3, time: -1}
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.
    t : ureg.Quantity {time: 1}
        time, beginning of release is at t=0.

    Returns
    -------
    float
        Oxygen concentration.
    """
    if Q > 0:
        C = 0.21/(Q+R) * (Q+R*math.e**-(Q+R)/V*t)
    elif abs(Q) <= R:
        C = 0.21*math.e**-(R/V*t)
    elif abs(Q) > R:
        C = 0.21*(1-R/abs(Q)*(1-math.e**-(abs(Q)*t/V)))
    return C


def conc_final(V, R, Q):
    """Calculate the final oxygen concentration for continuous flow.

    Equivalent to conc_vent(V, R, Q, float('inf')).

    Parameters
    ----------
    V : ureg.Quantity {length: 3}
        Volume of the confined space.
    R : ureg.Quantity {length: 3, time: -1}
        Volumetric spill rate into confined space.
    Q : ureg.Quantity {length: 3, time: -1}
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.

    Returns
    -------
    float
        Oxygen concentration.
    """
    if Q > 0:
        C = 0.21/(Q+R)*Q
    elif abs(Q) <= abs(R):
        C = 0
    elif abs(Q) > abs(R):
        C = 0.21*(1-R/abs(Q))
    return C


def conc_after(V, C_e, Q, t, t_e):
    """Calculate the oxygen concentration in the confined volume after
    the release has ended.

    As defined by FESHM 4240 6.1.A, Case D.

    Parameters
    ----------
    V : ureg.Quantity {length: 3}
        Volume of the confined space.
    C_e : float
        Oxygen concentration at the end of the release.
    Q : ureg.Quantity {length: 3, time: -1}
        Volumetric ventilation rate of fan(s); positive value corresponds
        to blowing air into the confined space, negative - drawing contaminated
        air outside.
    t : ureg.Quantity {time: 1}
        time, beginning of release is at t=0.
    t_e : ureg.Quantity {time: 1}
        time when release ended.

    Returns
    -------
    float
        Oxygen concentration.
    """
    C = 0.21-(0.21-C_e)*math.e**-(abs(Q)/V*(t-t_e))
    return C


def print_result(*Volumes):
    """Print the results of the ODH analysis for a volume.

    If several volumes given (in case of interlapping volumes) the worst case
    will be printed.
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
