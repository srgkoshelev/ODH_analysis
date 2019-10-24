#!python3
#Defining base classes for ODH analysis

import math, sys, logging
import heat_transfer as ht
from copy import copy

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)

#Setting up the units
ureg = ht.ureg
Q_ = ureg.Quantity
#ureg.auto_reduce_dimensions = True

#Loading FESHM 4240 Failure rates
from .FESHM4240_TABLES import *

#Probability of failure on demand for main cases
lambda_odh = 2.3e-6/ureg.hr #from CMTF Hi Bay ODH EN01878 pp. 27-28. That is the most correct value I have seen in use
#TODO State the source for ODH system PFD

class Source:
    """Define the possible source of inert gas"""
    def __init__(self, name, Fluid, volume, N=1):
        self.name = name
        self.Fluid = Fluid
        self.Leaks = {}
        self.N = N #Number of sources if multiple exist, e.g. gas cylinders. Increases probability of failure by N.
        #Calculating volume at standard conditions
        TempState = copy(Fluid)
        TempState.update('T', ht.T_NTP, 'P', ht.P_NTP)
        self.volume = volume*Fluid.Dmass/TempState.Dmass
        self.volume.ito(ureg.feet**3)
        self.isol_valve = False #By default assume there is no isolation valve that is used by ODH system

    def gas_pipe_failure(self, Pipe, Fluid=None, N_welds=1, max_flow=None):
        """
        Calculate failure rate, flow rate and expected time duration of the event for gas pipe failure. Based on FESHM 4240.
        """
        Fluid = Fluid or self.Fluid #If Fluid not defined use Fluid of the Source
        failure_rate_coeff = {'Piping': Pipe.L, 'Pipe weld': N_welds * Pipe.OD / Pipe.wall} #Failure rate coefficients; Piping failure rate is per unit of length, weld is dependent on number of welds, pipe OD and wall thickness
        for cause in ['Piping', 'Pipe weld']: #Piping and weld leaks as per Table 2
            for mode in ['Small leak', 'Large leak', 'Rupture']:
                if Pipe.D > 2 or mode is not 'Large leak': #Large leak only for D > 2"
                    name = f'{cause} {mode.lower()}: {Pipe}'
                    TempPipe = copy(Pipe)
                    TempPipe.L = Pipe.L / 2 #Average path for the flow will be half of piping length for gas piping
                    if mode == 'Rupture':
                        failure_rate = failure_rate_coeff[cause] * TABLE_2[cause][mode] #Leak failure rate from FESHM Table 2
                        area = Pipe.Area #for rupture calculate flow through available pipe area
                    else:
                        failure_rate = failure_rate_coeff[cause] * TABLE_2[cause][mode]['Failure rate'] #Leak failure rate from FESHM Table 2
                        area = TABLE_2[cause][mode]['Area'] #Leak area from FESHM Table 2
                    if max_flow is not None:
                        q_std_max = ht.piping.to_standard_flow(max_flow, Fluid)
                        q_std = min(self._leak_flow(TempPipe, area, Fluid), q_std_max)
                    else:
                        q_std = self._leak_flow(TempPipe, area, Fluid)
                    tau = self.volume/q_std
                    self.Leaks[name] = (failure_rate.to(1/ureg.hr), q_std, tau.to(ureg.min))

    def transfer_line_failure(Pipe, Fluid=None, N_lines=1):
        """
        Calculate failure rate, flow rate and expected time duration of the event for transfer line failure. Based on FESHM 4240.
        """
        area_cases = {'Leak': Q_('10 mm^2'), #Leak area, assumed based on piping leak areas in Table 2
                'Rupture': Pipe.Area, #Full flow through the transfer line
                }
        for mode in ['Leak', 'Rupture']:
            name = f'Fluid line {mode.lower()}: {Pipe}'
            failure_rate = N_lines * TABLE_1['Fluid line'][mode]
            area_cases = Area[mode]
            Fluid = Fluid or self.Fluid #If Fluid not defined use Fluid of the Source
            q_std = self._leak_flow(Pipe, area, Fluid)
            tau = self.volume/q_std
            self.Leaks[name] = (failure_rate.to(1/ureg.hr), q_std, tau.to(ureg.min))

    def failure_mode(name, failure_rate, flow_rate):
            tau = self.volume/flow_rate
            self.Leaks[name] = (failure_rate.to(1/ureg.hr), flow_rate, tau.to(ureg.min))

    def _leak_flow(self, Pipe, Area, Fluid):
        d = (4*Area/math.pi)**0.5 #diameter for the leak opening
        Entrance = ht.piping.Entrance(d)
        Exit = ht.piping.Exit(d)
        TempPiping = ht.piping.Piping(Fluid)
        TempPiping.add(Entrance,
                       Pipe,
                       Exit,
        )
        m_dot = TempPiping.m_dot(ht.P_NTP)
        return ht.piping.to_standard_flow(m_dot, Fluid)

    def __add__ (self, other):
        if self.Fluid.name == other.Fluid.name:
            total_volume = self.volume + other.volume
            TempState = ht.ThermState(name=self.Fluid.name, backend=self.Fluid.backend)
            TempState.update('T', ht.T_NTP, 'P', ht.P_NTP)
            return Source(None, '', TempState, total_volume) 
        else:
            logger.error ('\nBoth volumes should contain the same fluid')

    def __mul__(self, num):
        self.N *= num
        return self

    def __rmul__(self, num):
        self.N *= num
        return self

    def __str__ (self):
        if self.name:
            return self.name
        else:
            return ''

    def print (self):
        if self.name:
            print('{} is an ODH source of {} with volume of {:.3~}'.format(self.name, self.Fluid.name, self.volume))

    def print_leaks (self):
        for key in sorted(self.Leaks.keys()):
            print ('Failure mode: ' + key + ' pipe')
            print ('Failure rate: {:.2~}'.format(self.Leaks[key][0]))
            print ('Flow rate: {:.2~}'.format(self.Leaks[key][1].to(ureg.ft**3/ureg.min)))
            print ('Event duration: {:.2~}'.format(self.Leaks[key][2]))
            print ()


class Volume:
    '''
    Volume/building affected by inert gases.
    '''
    show_sens = 1e-7/ureg.hr
    def __init__  (self, name, Fluids, volume):
        self.name = name
        self.Fluids = Fluids #For gas stratification gases that affect the layer
        self.volume = volume

    def __str__ (self):
            return self.name

    def fan_fail (self, Test_period, Fail_rate, Q_fan, N_fans, Mean_repair_time=0*ureg.hr):
        '''
        Calculate (Probability, flow) pairs for all combinations of fans working. All fans are expected to have same volume flow
        '''
        Fan_flowrates = []
        for m in range(N_fans+1):
            P_m_fan_work = prob_m_of_n(m, N_fans, Test_period, Fail_rate, Mean_repair_time) #Probability of exactly m units starting
            Fan_flowrates.append((P_m_fan_work, Q_fan*m))
        self.Fan_flowrates = Fan_flowrates
        self.mdt = Test_period/(N_fans+1)+Mean_repair_time #Mean Down Time for fan system, only 1 fan required to work - used for constant leak calculation

    def PFD_system (self, *PFDs):
        '''
        Calculate total probability of an ODH system failure for given probabilities of different element failures.
        '''
        PFD_system_on = 1
        for PFD in PFDs:
            PFD_system_on *= 1-PFD 
        return 1-PFD_system_on 

    def fatality_prob(self, O2_conc):
        if O2_conc >= 0.18: #Lowest oxygen concentration above 18%
            Fi = 0
        elif O2_conc <= 0.088: #8.8% of oxygen is assumed to be 100% fatal
            Fi = 1
        else:
            Fi = 10**(6.5-76*O2_conc) #Fi formula, reverse engineered using 8.8% and 18% thresholds; These values are used in FESHM chapter to approximate O2 partial pressure
        return Fi

    def odh_class(self):
        if self.phi < 1e-7/ureg.hr:
            return 0
        elif self.phi < 1e-5/ureg.hr:
            return 1
        elif self.phi < 1e-3/ureg.hr:
            return 2
        else:
            logger.error ('ODH fatality rate is too high. Please, check calculations')

    def odh (self, Power = True):
        '''
        Calculate ODH fatality rate and recommend ODH class designation.
        Uses list of odh.source instances.
        Power specifies whether there is a power outage. Default is no outage.
        '''
        self.phi = 0 #fatality rate
        PFD_power_build = Power*PFD['power']+(not Power) #Probability of power failure in the building: PFD_power if no outage, 1 if there is outage
        for source in Source.instances: #Sources is a list of Source objects
            if source.fluid not in self.Fluids:
                continue
            if len (source.Leaks) == 0:
                if self.source_safe(source):
                    if Power: #to avoid excess output for power outage case
                        print ("{} doesn't contain enough gas to cause fatality in {}\n".format(source.name, self.name))
                    continue
                else:
                    logger.warning('Potentially dangerous for {} source "{}" has no defined leaks!\n'.format(self.name, source.name))
            PFD_sol_source = source.isol_valve*PFD['sol']+(not source.isol_valve) #If there is an isolation valve, probability = PFD_sol; if there are no valve the probability to fail = 1 i.e. the flow cannot be stopped
            for cause, leak in source.Leaks.items():
                (P_leak, q_leak, tau) = leak
                if hasattr(P_leak, 'dimensionality'):
                    P_i = P_leak*(PFD_power_build*PFD_sol_source+(1-PFD_power_build)*PFD['odh']) #power and solenoid failure or power is on and ODH system fails: situation, when leak continues and fans do not start
                    #The previous equation is equal to P_i = P_leak*self.PFD_system(PFD_power, PFD_odh) probability of power or ODH system failure: fans do not start
                    O2_conc = conc_vent (self.volume, q_leak, 0*ureg('ft^3/min'), tau) #is limited by ammount of inert gas the source has; fans are not operational
                    F_i = self.fatality_prob(O2_conc)
                    self.phi += P_i*F_i
                    self.info(source, self.show_sens, cause, O2_conc, P_leak, P_i, F_i, Power, q_leak)
                    if hasattr(self, 'Fan_flowrates'):
                        for (P_fan, Q_fan) in self.Fan_flowrates:
                            P_i = P_leak*(1-self.PFD_system(PFD_power_build, PFD['odh']))*PFD_sol_source*P_fan #Probability of leak occuring, ODH system/power working, solenoid not closing and m number of fans working; closed solenoid cuts off the supply thus phi = 0
                            #Above is equal to P_i = P_leak*(1-self.PFD_system(PFD_power, PFD_odh))*P_fan - Probability of leak occuring, ODH system/power working and m number of fans working
                            O2_conc = conc_vent (self.volume, q_leak, Q_fan, tau)
                            F_i = self.fatality_prob(O2_conc)
                            self.phi += P_i*F_i
                            self.info(source, self.show_sens, cause, O2_conc, P_leak, P_i, F_i, Power, q_leak, Q_fan)
                    else:
                        raise Exception ('Need to calculate Fan flowrates first')
                else:
                    #Assuming constant leak is small enough to be completely dominated by a single fan
                    lambda_const = lambda_power + lambda_odh + self.Fan_flowrates[0][0]/self.mdt #Fan failure rate = PFD_0fans/MDT_fans
                    P_i = lambda_const #Failure rate for protection against constant leak
                    O2_conc = conc_vent (self.volume, q_leak, 0*ureg('ft^3/min'), tau)
                    F_i = self.fatality_prob(O2_conc)
                    self.phi += P_i*F_i
                    self.info(source, self.show_sens, 'Constant leak', 0, 'Constant', P_i, F_i, Power, q_leak, 0*ureg('ft^3/min'))
                    #TODO Constant leak analysis using HVAC ACH

    def info(self, source, sens, cause, O2_conc, P_leak, P_i, F_i, Power, q_leak, Q_fan=None):
        phi = P_i*F_i
        if phi >= sens and Power:
            print ('Major ODH cause with fatality rate > {:.2~}'.format(sens))
            print ('Volume: {}, {:.2~}'.format(self, self.volume.to(ureg.ft**3)))
            print ('Source: '+source.name)
            if q_leak and Q_fan:
                print ('Failure: '+cause)
            else:
                print ('Failure: '+cause+' and No power')
            print ('Oxygen concentration: {:.0%}, {:.0%} percent of norm'.format(O2_conc, O2_conc/0.21))
            if type(P_leak) == str:
                print ('Leak prob rate: {}'.format(P_leak))
            else:
                print ('Leak prob rate: {:.2~}'.format(P_leak))
            print ('System prob: {:.2%}'.format(P_i/P_leak))
            print ('Failure rate: {:.2~}'.format(P_i))
            print ('Leak rate: {:.2~}'.format(q_leak))
            if Q_fan:
                print ('Fan rate: {:.2~}'.format(Q_fan))
            print ('Fatality prob: {:.2g}'.format(F_i))
            print ('Fatality rate: {:.2~}\n'.format(P_i*F_i))


    def source_safe(self, source, escape = True):
        '''
        Estimate the impact of the Source vokume on oxygen concetration. Smaller sources might not be able to drop oxygen concentration to dangerous levels.
        '''
        if escape == True: #if mixed air is allowed to escape within considered volume
            O2_conc = 0.21*self.volume/(self.volume+source.volume)
        else: #worst case; inert gas is trapped and expells the air outside the considered volume
            O2_conc = 0.21*(1-source.volume/self.volume)
        return self.fatality_prob(O2_conc) == 0


def prob_m_of_n (m, n, T, l, MRT=0*ureg.hr):
    '''
    Calculate the probability of m out of n units working.
    inputs:
        T - test period, hr
        MRT - mean repair time (when not negligible compared to test period), hr
        l = lambda = failure rate for fan, 1/hr
    '''
    C_n_m = math.factorial(n)/(math.factorial(n-m)*math.factorial(m))
    PFD_one_unit = l*T
    F_adj = 1/(n-m+1)+MRT/T #Adjustment coefficient: T/(n-m+1) will be average failure reveal time (D. Smith, Reliability..., p. 108); MRT is added to reveal time to get total unavailability time
    m_of_n = C_n_m*(PFD_one_unit)**(n-m)*(1-PFD_one_unit)**m*F_adj #see ED00007314 for details
    return m_of_n

def conc_vent (V, R, Q, t):
    #V - volume of the confined space (ft3 or m3)
    #R - spill rate into confined space (scfm or m3/s)
    #Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds to blowing air into the confined space, negative - drawing contaminated air outside
    #t = time, (minutes or seconds) beginning of release is at t=0
    #C - oxygen concentration in confined space
    #Case B
    V = V.to(ureg.m**3).magnitude #There seems to be a bug with unit package regarding epxonentiation
    R = R.to(ureg.m**3/ureg.s).magnitude
    Q = Q.to(ureg.m**3/ureg.s).magnitude
    t = t.to(ureg.s).magnitude
    if Q > 0:
        C = 0.21/(Q+R)*(Q+R*math.e**-((Q+R)/V*t))
    elif abs(Q) <= R:
        C = 0.21*math.e**-(R/V*t)
    elif abs(Q) > R:
        C = 0.21*(1-R/abs(Q)*(1-math.e**-(abs(Q)*t/V)))
    return C


def conc_final (V, R, Q):
    #V - volume of the confined space (ft3 or m3)
    #R - spill rate into confined space (scfm or m3/s)
    #Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds to blowing air into the confined space, negative - drawing contaminated air outside
    #C - oxygen concentration in confined space
    #Case B
    if Q > 0*ureg('ft^3/min'):
            C = 0.21/(Q+R)*Q
    elif abs(Q) <= abs(R):
            C = 0
    elif abs(Q) > abs(R):
            C = 0.21*(1-R/abs(Q))
    return C

def conc_after (V, C_e, Q, t, t_e):
    #V - volume of the confined space (ft3 or m3)
    #R - spill rate into confined space (scfm or m3/s)
    #Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds to blowing air into the confined space, negative - drawing contaminated air outside
    #t = time, (minutes or seconds) beginning of release is at t=0
    #C - oxygen concentration in confined space
    #C_e = oxygen concentration when the release has ended
    C = 0.21-(0.21-C_e)*np.e**-(abs(Q)/V*(t-t_e))
    return C



#Printing results - may be moved elsewhere later
def print_result(*Volumes):
    '''
    Print the results of the ODH analysis for a volume. If several volumes given (in case of interlapping volumes) the worst case will be printed.
    '''
    max_phi = -1/ureg.hr
    for volume in Volumes:
        if volume.phi > max_phi:
            max_volume = volume
    line_1 = '# Fatality rate for {} is {:.1e} #'.format(max_volume, volume.phi)
    pad = len(line_1)
    line_2 = '# Recommended ODH class {}'.format(max_volume.odh_class()).ljust(pad-1)+'#'
    print ('#'*pad)
    print (line_1)
    print (line_2)
    print ('#'*pad)
