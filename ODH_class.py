#!python3
#Defining base classes for ODH analysis

import math, sys, logging
sys.path.append('D:/Personal/Python repo/')
from heat_transfer import functions as ht
from heat_transfer.piping import *

logging.basicConfig(format='%(levelname)s:%(message)s', level=logging.WARNING)
logger = logging.getLogger(__name__)

#Setting up the units
ureg = ht.ureg
Q_ = ureg.Quantity
#ureg.auto_reduce_dimensions = True

#Setting units for "standard" flow
T_NTP = Q_(68, ureg.degF) #Normal Temperature (NIST)
P_NTP = Q_(14.7, ureg.psi) #Normal Pressure (NIST)

T_MSC = Q_(15, ureg.degC) #Metric Standard Conditions (used by Crane TP-410)
P_MSC = Q_(101325, ureg.Pa) #Metric Standard Conditions (used by Crane TP-410)

#Probability of failure on demand for main cases
PFD_power = 1e-4 #For some reason FESHM chapter lists 3e-4 as demand rate and 1e-4/hr failure rate + 1 hr time off. D. Smith's Reliability... states that PFD = lambda*MTD with lambda = 1e-4 1/hr and MDT = 1 hr gives 1e-4 value.
lambda_power = 1e-4/ureg.hr #Electrical power failure rate; Used for constant leaks
PFD_odh = 2e-3 #Conservative etimate; from CMTF Hi Bay ODH EN01878 pp. 27-28. That is the most correct value I have seen in use
lambda_odh = 2.3e-6/ureg.hr #from CMTF Hi Bay ODH EN01878 pp. 27-28. That is the most correct value I have seen in use
PFD_sol = 1e-3 #Solenoid valve, failure to operate; Used as failure to close when not powered (usually solenoids used for isolating the source)

#Gas pipe safety factor
Fs_gas = 3 #Gas pipe leak probability is calculated using length of piping and number of welds. Unfortunately both values are hard to estimate accurately thus a safety factor is used

class odh_source:
    """Define the possible source of inert gas"""
    instances = [] #Keeping information on all the sources within the class
    def __init__ (self, name, fluid, Volume, phase = 'gas', pressure = P_NTP): 
        odh_source.instances.append(self) #adding initialized instance to instances list
        self.name = name
        self.fluid = fluid
        self.Leaks = {}
        assert phase in ['vapor','liquid', 'gas'], 'Phase can only be liquid or vapor: %r' % phase
        if phase == 'vapor' or phase == 'gas':
            self.volume = Volume*pressure/Q_(14.7,ureg.psi)
        elif phase == 'liquid':
            Fluid_data = ht.pack_fluid(fluid, T_NTP, P_NTP) #standard conditions
            (x, M, D_fluid_std) = ht.rp_init(Fluid_data)
            satur = ht.satp(P_NTP, x) #assuming incompressible liquid
            D_fluid_sat = satur ['Dliq']*ureg('mol/L')
            self.volume = Volume*D_fluid_sat/D_fluid_std
        self.volume.ito(ureg.feet**3)
        self.isol_valve = False #Shows if there is an isolation valve that is used by ODH system

    def __add__ (self, other):
        if self.fluid == other.fluid:
            total_volume = self.volume + other.volume
            return odh_source(None, self.fluid, total_volume, 'gas', P_NTP) 
        else:
            logger.error ('\nBoth volumes should contain the same fluid')

    def __str__ (self):
        if self.name:
            return self.name
        else:
            return ''

    def print (self):
        if self.name:
            print('{} is an ODH source of {} with volume of {:.3~}'.format(self.name, self.fluid, self.volume))
    
    def delete(*Sources):
        '''
        Helper function. Remove unused and autogenerated instances of sources.
        '''
        odh_source.instances = [x for x in odh_source.instances if not (x in Sources or x.name == None)] 

    def calculate_gas_leak (self, cause, failure_mode):
        '''
        Helper function for different cases of gas pipe/weld leaks and ruptures
        '''
        pipe = failure_mode['Pipe']
        N_welds = failure_mode.get('N_welds')
        if 'large' in cause and pipe.D<=2:
            return
        cause += ' {}'.format(str(pipe.D))
        Prob = pipe.L/(ureg.m*ureg.hr)
        if 'small' in cause:
            Area = 10*ureg.mm**2
            Prob *= 1e-9
        elif 'large' in cause:
            Area = 1000*ureg.mm**2
            Prob *= 1e-10
        elif 'rupture' in cause:
            Area = pipe.Area
            Prob *= 3e-11
        else: 
            logger.error ('Gas pipe failure cause is not recognized!:{}'.format(cause))

        if 'weld' in cause:
            Prob = N_welds*pipe.OD/(pipe.wall*ureg.hr)
            if 'small' in cause:
                Prob *= 2e-11
            elif 'large' in cause:
                Prob *= 2e-12
            elif 'rupture' in cause:
                Prob *= 6e-13

        q_std = self.leak_flow(failure_mode, Area)
        tau = self.volume/q_std
        Prob *= Fs_gas #Including safety factor due to uncertainty in piping length/weld number estimation
        Prob.ito(1/ureg.hr)
        self.Leaks[cause] = (Prob, q_std, tau.to(ureg.min))


    def leak (self, failure_mode = {'mode': 'gas line', 'Pipe':Pipe(3, 5, 10*ureg.m), 'Fluid_data': {'P':2.33*ureg('bar'), 'T':Q_(40,ureg.degC)}, 'max_flow':0.01*ureg('m^3/s')}):
        '''Calculating leak and probability for different cases of equipment failure
        '''
        Fluid_data = failure_mode.get('Fluid_data')
        if Fluid_data:
            (fluid, T_fluid, P_fluid) = ht.unpack_fluid(failure_mode.get('Fluid_data'))
            fluid = self.fluid #is defined by the source
            Fluid_data = ht.pack_fluid(fluid, T_fluid, P_fluid)
        else:
            Fluid_data = {'fluid':self.fluid}
        failure_mode['Fluid_data'] = Fluid_data
        pipe = failure_mode.get('Pipe')

        if failure_mode['mode'] == 'gas line':
            Causes = ['small leak', 'large leak', 'rupture', 'weld small leak', 'weld large leak', 'weld rupture']
            for cause in Causes:
                self.calculate_gas_leak(cause, failure_mode)

        elif failure_mode['mode'] == 'fluid line':
            N_lines = failure_mode['N_lines'] #For multiple tranfer lines/U-tubes used the min length of pipe should be used
            for cause in ['fluid line leak', 'fluid line rupture']:
                if 'leak' in cause:
                    Prob = N_lines*5*10**(-7)/(ureg.hr) #Probaility and flow will be recalculated for each cause using the same variable names
                    Area = 10*ureg.mm**2
                else:
                    Prob = N_lines*4*10**(-8)/(ureg.hr) #FESHM chapter uses unconservative approach with ~60% confidence; This value gives 90% confidence
                    Area = pipe.Area
            q_std = self.leak_flow(failure_mode, Area)
            tau = self.volume/q_std
            Prob.ito(1/ureg.hr)
            self.Leaks[cause] = (Prob, q_std, tau.to(ureg.min))

        elif failure_mode['mode'] == 'dewar':
            cause = 'loss of vacuum to air'
            Prob = 4*10**(-6)/(ureg.hr) #FESHM chapter uses unconservative approach with ~60% confidence; This value gives 90% confidence
            q_std = self.limit_flow(failure_mode['q_relief'], failure_mode)
            tau = self.volume/q_std
            Prob.ito(1/ureg.hr)
            self.Leaks[cause] = (Prob, q_std, tau.to(ureg.min))

        elif failure_mode['mode'] == 'other':
            cause = failure_mode['cause']
            Prob = failure_mode['Prob'] #Probability of a single failure
            Prob_total = Prob*failure_mode.get('N', 1) #Total probability of failure
            try:
                Prob.ito(1/ureg.hr)
            except AttributeError:
                if Prob != 1:
                    logger.error('Flat probability not equal to 1 (constant leak) is used: {}'.format(Prob))
            q_std = self.limit_flow(failure_mode['q'], failure_mode)
            tau = self.volume/q_std
            self.Leaks[cause] = (Prob, q_std, tau.to(ureg.min))

        else:
            raise  ValueError ('Mode is not supported: %r. Try "gas line", "fluid line", "dewar" or "other".', failure_mode['mode'])


    def leak_flow(self, failure_mode, Area):
        d = (4*Area/math.pi)**0.5 #diameter for the leak openning
        pipe = failure_mode['Pipe'] #For fluid line whole length is affecting the flow; the length should be to the nearest bayonet from the source
        if failure_mode['mode'] == 'gas line':
            pipe.L = pipe.L/2 #Average path for the flow will be half of piping length for gas piping
        openning = Openning(d)
        piping = Piping(pipe, failure_mode['Fluid_data'])
        piping.add(openning)
        m_dot = piping.m_dot(P_NTP)
        return self.limit_flow(m_dot, failure_mode)

    def limit_flow(self, flow, failure_mode):
        '''
        Leak through the openning sometimes cannot be larger than a certain number, e.g. a compressor throughput. This function limits flow to this value if it exists.
        Max flow should be specified for flowing conditions (Fluid_data).
        Additionally converts the flow to standard conditions (required for comparison).
        '''
        max_flow = failure_mode.get('max_flow')
        #logger.debug('Failure mode: {} Flow: {:.3~}, STD flow: {:.3~}'.format(failure_mode.get('cause', failure_mode['mode']), flow, to_standard_flow(flow, failure_mode['Fluid_data']) ))
        flow = to_standard_flow(flow, failure_mode['Fluid_data'])
        if max_flow:
            max_flow = to_standard_flow(max_flow, failure_mode['Fluid_data']) #converting max_flow to standard volumetric flow  
            return min (flow, max_flow)
        else:
            return flow

    def print_leaks (self):
        for key in sorted(self.Leaks.keys()):
            print ('Failure mode: ' + key + ' pipe')
            print ('Failure rate: {:.2~}'.format(self.Leaks[key][0]))
            print ('Flow rate: {:.2~}'.format(self.Leaks[key][1].to(ureg.ft**3/ureg.min)))
            print ('Event duration: {:.2~}'.format(self.Leaks[key][2]))
            print ()


    

class odh_volume:
    '''
    Volume/building affected by inert gases.
    '''
    show_sens = 1e-7
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
        PFD_prev = 1
        for m in reversed(range(1, N_fans+1)):
            PFD_fan = failure_on_demand(m, N_fans, Test_period, Fail_rate, Mean_repair_time) #Probability of failure if m units is required; Equivalent to 0..m units working
            PFD_m_fan_work = PFD_prev - PFD_fan #Probability exactly m fans work = P(0..m+1) - P(0..m) = PFD_prev-PFD_fan
            Fan_flowrates.insert(0, (PFD_m_fan_work, Q_fan*m))
            PFD_prev = PFD_fan
        Fan_flowrates.insert(0, (PFD_fan, 0*ureg('ft^3/min'))) #last value, probability of system failure with only 1 required to work = probability of exactly 0 fans working
        unity_check = sum([x for x,y in Fan_flowrates]) 
        if unity_check != 1:
            logger.error('Unity check for fan probabilities failed: {}'.format(unity_check))
        self.Fan_flowrates = Fan_flowrates
        self.mdt = Test_period/(N_fans-1)+Mean_repair_time #Mean Down Time for fan system, only 1 fan required to work - used for constant leak calculation

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
            Fi = 10**(6.6956-76.087*O2_conc) #Fi formula, reverse engineered using 8.8% and 18% thresholds; These values are used in FESHM chapter to approximate O2 partial pressure
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
        PFD_power_build = Power*PFD_power+(not Power)*1 #Probability of power failure in the building: PFD_power if no outage, 1 if there is outage
        for source in odh_source.instances: #Sources is a list of odh_source objects
            if source.fluid not in self.Fluids:
                continue
            if len (source.Leaks) == 0:
                if self.source_safe(source):
                    if Power: #to avoid excess output for power outage case
                        print ("{} doesn't contain enough gas to cause fatality in {}\n".format(source.name, self.name))
                    continue
                else:
                    logger.warning('Potentially dangerous for {} source "{}" has no defined leaks!\n'.format(self.name, source.name))
            PFD_sol_source = source.isol_valve*PFD_sol+(not source.isol_valve)*1 #If there is an isolation valve, probability = PFD_sol; if there are no valve the probability to fail = 1 i.e. the flow cannot be stopped
            for cause, leak in source.Leaks.items():
                (P_leak, q_leak, tau) = leak
                if hasattr(P_leak, 'dimensionality'):
                    P_i = P_leak*(PFD_power_build*PFD_sol_source+(1-PFD_power_build)*PFD_odh) #power and solenoid failure or power is on and ODH system fails: situation, when leak continues and fans do not start
                    #The previous equation is equal to P_i = P_leak*self.PFD_system(PFD_power, PFD_odh) probability of power or ODH system failure: fans do not start
                    O2_conc = conc_vent (self.volume, q_leak, 0*ureg('ft^3/min'), tau) #is limited by ammount of inert gas the source has; fans are not operational
                    F_i = self.fatality_prob(O2_conc)
                    self.phi += P_i*F_i
                    self.info(source, self.show_sens, cause, O2_conc, P_leak, P_i, F_i, Power, q_leak)
                    if hasattr(self, 'Fan_flowrates'):
                        for (P_fan, Q_fan) in self.Fan_flowrates:
                            P_i = P_leak*(1-self.PFD_system(PFD_power_build, PFD_odh))*PFD_sol_source*P_fan #Probability of leak occuring, ODH system/power working, solenoid not closing and m number of fans working; closed solenoid cuts off the supply thus phi = 0
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
        if phi >= sens/ureg.hr and Power:
            print ('Major ODH cause with fatality rate > {:.2~}'.format(sens/ureg.hr))
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

        



def failure_on_demand (m, n, T, l, MTTR=0*ureg.hr):
    '''
    Failure on demand probability
    The definition in FESHM chapter (probability of m units out of n starting) is incorrect. Definition in D. Smith's Reliability... is confusing.
    This value represents the probability of system failure if for normal functioning only m units out of n is required.
    1 - PFD gives the probability of at least m units out of n starting. That is usually not what required for ODH analysis. The workaround should be used.
    inputs:
        T - test period, hr
        MTTR - mean repair time (when not negligible compared to test period), hr
        l = lambda = failure rate for fan, 1/hr
    '''
    PFD = math.factorial(n)*(l*T)**(n-m+1)/(math.factorial(m-1)*math.factorial(n-m+2))*(1+(n-m+2)*MTTR/T) #Adapted formula from D. Smith's Reliability for unrevealed failures p. 108 to include repair time effect
    return PFD

def to_standard_flow(flow_rate, Fluid_data):
    '''
    Converting volumetric flow at certain conditions or mass flow to volumetric flow at NTP
    '''
    (x, M, D_NTP) = ht.rp_init({'fluid':Fluid_data['fluid'], 'T':T_NTP, 'P':P_NTP})
    if flow_rate.dimensionality == ureg('kg/s').dimensionality: #mass flow, flow conditions are unnecessary
        q_std = flow_rate/(D_NTP*M)
    elif flow_rate.dimensionality == ureg('m^3/s').dimensionality: #volumetric flow given, converting to standard pressure and temperature
        if 'T' in Fluid_data and 'P' in Fluid_data:
            (fluid, T_fluid, P_fluid) = ht.unpack_fluid(Fluid_data)
            (x, M, D_fluid) = ht.rp_init(Fluid_data)
            q_std = flow_rate*D_fluid/D_NTP
        else:
            logger.warning('Flow conditions for volumetric flow {:.3~} are not set. Assuming standard flow at NTP'.format(flow_rate))
            q_std = flow_rate
    q_std.ito(ureg.ft**3/ureg.min)
    return q_std

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










    
if __name__ == "__main__":
    #Testing

    He_storage_dewar_gas = odh_source('Storage dewar', 'helium', Q_(3390000, ureg.cubic_feet), 'vapor', Q_(0, ureg.psig)) #blowdown from 12 psig to 1 atmosphere, estimated by R. Rabehl, TID-N-3A, p. 7
    Test = odh_source('Test', 'helium', Q_(33900, ureg.cubic_feet), 'vapor', Q_(0, ureg.psig)) #blowdown from 12 psig to 1 atmosphere, estimated by R. Rabehl, TID-N-3A, p. 7
    Test_2 = Test+Test+Test
    odh_source.delete()
    He_storage_dewar_gas.leak({'mode':'gas line', 'Pipe':Tube(6*ureg.inch, 0.15*ureg.inch, 350000*ureg.ft), 'N_welds':10, 'Fluid_data':{'P': 95*ureg.psig, 'T':300*ureg.K}, 'max_flow':108*ureg('g/s')}) #Mid stage; Max flow is taken as midstage return from cryoplant
    Test_pipe = Tube(1*ureg.inch, 0.035*ureg.inch, 20*ureg.ft)
    Test_piping = Piping(Test_pipe, {'fluid':'helium', 'P':12*ureg.psig, 'T':5*ureg.K})
    print(Test_piping.m_dot().to(ureg.g/ureg.s))
    print(Test_piping.dP(607*ureg('g/s')))
    Test_period = 1*ureg('month') #IB1 fan test period 
    Test_period.ito(ureg.hr)
    Mean_repair_time = 3*ureg('days') #IB1 fan/louver average repair time: most delay is caused by response time, it has recently improved from about 1 week to 1 day. The average value of 3 days is assumed
    l_fan = 9e-6/ureg.hr #Fan failure rate
    l_louv = 3e-7/ureg.hr #Louver failure rate
    l_vent = l_fan+l_louv #At IB1 louver and fan are installed in series; failure of either one results in no venting
    Q_fan = 4000*ureg.ft**3/ureg.min #Flowrate of 4 ceiling fans is >= 4000 CFM; positive value refers to blowing into the building: FALSE IN CASE OF IB1!
    N_fans = 4

    A = Q_(15360, ureg.square_feet) #IB1 floor area
    IB1_air = odh_volume('IB1 air', ['helium', 'nitrogen'], A*19.8*ureg.feet) #All IB1 air

    IB1_air.fan_fail(Test_period, l_vent, Q_fan, N_fans, Mean_repair_time)
    IB1_air.odh()
    print_result(IB1_air)

    tau = list(range(math.ceil(266300*60/7750)))*ureg.s
    Q = -16000*ureg('ft^3/min')
    R =7750*ureg('ft^3/min')
    tau = 266300*ureg.cubic_feet/R
    H = 9.8*ureg.ft
    V = H*A
    C = conc_vent(V, R, Q, tau)
    Test_vol = odh_volume ('Test', [],V)
    print ('Volume: {:.0f}'.format(V.to(ureg.ft**3)))
    print ('Oxygen concentration: {:.0%}'.format(C))
    print ('Oxygen pressure: {:.0f}'.format(C/0.21*160)) 
    print ('Fatality factor: {:.0g}'.format(Test_vol.fatality_prob(C)))
