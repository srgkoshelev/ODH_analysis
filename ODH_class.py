#!python3
#Defining base classes for ODH analysis

import math
import sys
sys.path.append('D:/Personal/Python repo/')
from heat_transfer import functions as ht
from heat_transfer import piping as pipe

ureg = ht.ureg
Q_ = ureg.Quantity
#ureg.auto_reduce_dimensions = True



T_NTP = Q_(68, ureg.degF) #Normal Temperature (NIST)
P_NTP = Q_(14.7, ureg.psi) #Normal Pressure (NIST)

T_MSC = Q_(15, ureg.degC) #Metric Standard Conditions (used by Crane TP-410)
P_MSC = Q_(101325, ureg.Pa) #Metric Standard Conditions (used by Crane TP-410)


class odh_source:
    """Define the possible source of inert gas"""
    def __init__ (self,fluid, Volume, phase = 'vapor', pressure = Q_(0, ureg.psig) ):
        self.fluid = fluid
        self.Leaks = {}
        assert phase in ['vapor','liquid', 'gas'], 'Phase can only be liquid or vapor: %r' % phase
        if phase == 'vapor' or phase == 'gas':
            self.volume = Volume*pressure/Q_(14.7,ureg.psi)
        elif phase == 'liquid':
            Fluid_data = ht.pack_fluid(fluid, T_NTP, P_NTP) #standard conditions
            (x, M, D_fluid_std) = ht.rp_init(Fluid_data)
            satur = ht.satp(101325*ureg.Pa, x)
            #D_fluid_sat = Q_(satur ['Dliq'], ureg('mol/L'))
            D_fluid_sat = satur ['Dliq']*ureg('mol/L')
            self.volume = Volume*D_fluid_sat/D_fluid_std
            #print (D_fluid_sat,D_fluid_std)
        self.volume.ito(ureg('ft**3'))

    def safe (self, escape = True):
        A = Q_(1550, ureg.square_feet) #IB1 floor area
        He_height = [3.3, 9.8, 19.8]*ureg.feet
        N2_height = 0.6*ureg.m #height of a sitting person
        if self.fluid == 'helium':
            V_effect = [A*h for h in He_height]
        elif self.fluid == 'nitrogen':
            V_effect = [A*N2_height] #Volume affected by inert gas
        if escape == True: #if mixed air is allowed to escape within considered volume
            O2_conc = [0.21*V/(V+self.volume) for V in V_effect]
        else: #worst case; inert gas is trapped and expells the air outside the considered volume
            O2_conc = [0.21*(1-self.volume/V) for V in V_effect]
        #print ( [self.fatality_prob(conc) == 0 for conc in O2_conc])
        return all([self.fatality_prob(conc) == 0 for conc in O2_conc])

        
    def fatality_prob(self, O2_conc):
        if O2_conc >= 0.18: #Lowest oxygen concentration above 18%
            Fi = 0
        elif O2_conc <= 0.088: #8.8% of oxygen is assumed to be 100% fatal
            Fi = 1
        else:
            Fi = 10**(6.6956-76.087*O2_conc) #Fi formula, reverse engineered using 8.8% and 18% thresholds; These values are used in FESHM chapter to approximate O2 partial pressure
        return Fi

    def leak (self, failure_mode = {'mode': 'gas line', 'Pipe':{'L':10*ureg.m, 'D_nom':3*ureg.inch, 'SCH':5}, 'Fluid_data': {'P':2.33*ureg('bar'), 'T':Q_(40,ureg.degC)}, 'max_flow':0.01*ureg('m^3/s')}):
        '''Calculating leak and probability for different cases of equipment failure
        '''
        (fluid, T_fluid, P_fluid) = ht.unpack_fluid(failure_mode['Fluid_data'])
        fluid = self.fluid #is defined by the source
        Fluid_data = ht.pack_fluid(fluid, T_fluid, P_fluid)
        failure_mode['Fluid_data'] = Fluid_data
        Pipe = failure_mode['Pipe']

        if failure_mode['mode'] == 'gas line':
            Pipe = failure_mode['Pipe']
            cause = 'small leak ' + str(Pipe['D_nom'])
            Prob = 10**(-9)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
            Prob.ito(1/ureg.hr)
            q = leak_flow(Fluid_data,  10*ureg.mm**2)
            self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            if pipe.D_pipe(Pipe) > 2:
                cause = 'large leak ' + str(Pipe['D_nom'])
                Prob = 10**(-10)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
                q = leak_flow(Fluid_data,  1000*ureg.mm**2)
                Prob.ito(1/ureg.hr)
                self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            cause = 'rupture ' + str(Pipe['D_nom'])
            Prob = 3*10**(-11)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
            Prob.ito(1/ureg.hr)
            q = leak_flow(Fluid_data,  pipe.Area(Pipe))
            self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            cause = 'weld small leak ' + str(Pipe['D_nom'])
            N_welds = failure_mode['N_welds']
            Prob = N_welds*2*10**(-11)/(ureg.hr)*pipe.OD(Pipe)/pipe.wall(Pipe) #Probaility and flow will be recalculated for each cause using the same variable names
            Prob.ito(1/ureg.hr)
            q = leak_flow(Fluid_data,  10*ureg.mm**2)
            self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            if pipe.D_pipe(Pipe) > 2:
                cause = 'weld large leak ' + str(Pipe['D_nom'])
                Prob = N_welds*2*10**(-12)/(ureg.hr)*pipe.OD(Pipe)/pipe.wall(Pipe) #Probaility and flow will be recalculated for each cause using the same variable names
                Prob.ito(1/ureg.hr)
                q = leak_flow(Fluid_data,  1000*ureg.mm**2)
                self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            cause = 'weld rupture ' + str(Pipe['D_nom'])
            Prob = N_welds*6*10**(-13)/(ureg.hr)*pipe.OD(Pipe)/pipe.wall(Pipe) #Probaility and flow will be recalculated for each cause using the same variable names
            Prob.ito(1/ureg.hr)
            q = leak_flow(Fluid_data,  pipe.Area(Pipe))
            self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

        elif failure_mode['mode'] == 'fluid line':
            cause = 'fluid line leak'
            N_lines = failure_mode['N_lines']
            Prob = N_lines*5*10**(-7)/(ureg.hr) #Probaility and flow will be recalculated for each cause using the same variable names
            Prob.ito(1/ureg.hr)
            q = leak_flow(Fluid_data,  10*ureg.mm**2)
            self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            cause = 'fluid line rupture'
            Prob = N_lines*5*10**(-8)/(ureg.hr) #FESH chapter uses unconservative approach with ~60% confidence; This value gives 95% confidence and is widely used: see rule of 3
            Prob.ito(1/ureg.hr)
            q = leak_flow(Fluid_data,  pipe.Area(Pipe))
            self.Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

        elif failure_mode['mode'] == 'dewar':
            cause = 'loss of vacuum to air'
            Prob = 1*10**(-6)/(ureg.hr) #Probaility and flow will be recalculated for each cause using the same variable names
            Prob.ito(1/ureg.hr)
            q = failure_mode['q_relief']
            self.Leaks[cause] = (Prob, q)

        elif failure_mode['mode'] == 'other':
            cause = failure_mode['cause']
            Prob = failure_mode['Prob']
            Prob.ito(1/ureg.hr)
            q = failure_mode['q']
            self.Leaks[cause] = (Prob, q)

        else:
            raise  ValueError ('Mode is not supported: %r. Try "gas line", "fluid line", "dewar" or "other".', failure_mode['mode'])

        return self.Leaks

    def limit_flow(self, flow, failure_mode):
        '''
        Leak through the openning sometimes cannot be larger than a certain number, e.g. a compressor throughput. This function limits flow to this value if it exists
        '''
        max_flow = failure_mode.get('max_flow', float('inf')*ureg.ft**3/ureg.min) #If max_flow value is not specified infinity is assumed
        return min (flow, max_flow)


#class odh_failure(odh_source):
#    """Define failure mode for existing source of inert gas"""
#    def __init__ (self, ODH_source, 


def failure_on_demand (m, n, Test_period, failure_rate):
    '''Failure on demand probability
    probability of m units starting out of n
    '''
    F_fod = ((Test_period/2*failure_rate*m)**(n+1-m))/math.factorial(n+1-m)
    return F_fod

def to_standard_flow(flow_rate, Fluid_data):
    (fluid, T_fluid, P_fluid) = ht.unpack_fluid(Fluid_data)
    (x, M, D_fluid) = ht.rp_init(Fluid_data)
    D_NTP = ht.flsh('TP', T_NTP, P_NTP, x)['D'] #Molar density of the fluid at NTP
    if flow_rate.dimensionality == ureg('m^3/s').dimensionality:
        q_std = flow_rate*D_fluid/D_NTP
    elif flow_rate.dimensionality == ureg('kg/s'):
        q_std = flow_rate/D_NTP
    q_std.ito(ureg.m**3/ureg.s)
    return q_std



def leak_flow (Fluid_data = {'fluid':'air', 'P':2.33*ureg('bar'), 'T':Q_(40,ureg.degC)}, A = 10*ureg('mm**2')):
    d = (4*A/math.pi)**0.5
    Y = 1 #conservative value; from Crane TP-410 A-21
    C = 0.7 #conservative value; from Crane TP-410 A-20 
    (fluid, T_fluid, P_fluid) = ht.unpack_fluid(Fluid_data)
    (x, M, D_fluid) = ht.rp_init(Fluid_data)
    rho = D_fluid*M
    k = ht.gamma(Fluid_data) #adiabatic coefficient
    rc = (2/(k+1))**(k/(k-1)) #Critical pressure drop; Note: according to Crane TP-410 is depndent on the hydraulic resistance of the flow path
    if P_MSC > P_fluid*rc: #subsonic flow
        DeltaP = P_fluid - P_MSC
    else: #Sonic flow
        DeltaP = P_fluid*(1-rc) #Crane TP-410, p 2-15
    quality = ht.flsh('TP', T_fluid, P_fluid, x)['q']
    if quality >= 1:
        (x_a, M_a) = ht.rp_init({'fluid':'air'})
        S_g = M/M_a #Specific gravity
        q = orifice_flow_gas(Y, d, C, S_g, DeltaP, rho)
    else:
        q = orifice_flow_liquid(d, C, DeltaP, rho) #Calculating liquid flow through the opening. Assuming no evaporation happens

    q_NTP = to_standard_flow(q, Fluid_data) #convert to NTP standard point
    q_NTP.ito(ureg.cubic_feet/ureg.min)
    return q_NTP

def orifice_flow_gas (Y, d, C, S_g, DeltaP, rho):
    '''Wrapper function for flow through orifice equation.
    Original formula in crane uses a non-dimensionless coefficient 0.0002864 and non-SI units for values in formulas so straightforward approach is complicated.
    '''
    d_1 = d.to(ureg.mm).magnitude #diameter in mm
    Deltap = DeltaP.to(ureg.bar).magnitude #Pressure drop in bar
    rho_1 = rho.to(ureg.kg/ureg.m**3).magnitude #Density of inflow in kg/m^3
    q = 0.0002864*Y*d_1**2*C/(S_g)*(Deltap*rho_1)**0.5 #Flow at MSC in m^3/s; Crane TP-410, Eq.3-22
    q = q*ureg('m^3/s')
    return q


def orifice_flow_liquid (d, C, DeltaP, rho):
    '''Wrapper function for flow through orifice equation.
    Original formula in crane uses a non-dimensionless coefficient 0.0002864 and non-SI units for values in formulas so straightforward approach is complicated.
    '''
    d_1 = d.to(ureg.mm).magnitude #diameter in mm
    Deltap = DeltaP.to(ureg.bar).magnitude #Pressure drop in bar
    rho_1 = rho.to(ureg.kg/ureg.m**3).magnitude #Density of inflow in kg/m^3
    q = 0.0003512*d_1**2*C*(Deltap/rho_1)**0.5 #Flow at MSC in m^3/s; Crane TP-410, Eq.3-21
    q = q*ureg('m^3/s')
    return q






#Adding the inert gas sources

He_storage_dewar_liquid = odh_source('helium', Q_(10000, ureg.L), 'liquid')
He_storage_dewar_gas = odh_source('helium', Q_(33900, ureg.cubic_feet), 'vapor', Q_(0, ureg.psig)) #blowdown from 12 psig to 1 atmosphere, estimated by R. Rabehl, TID-N-3A, p. 7
Buffer_tank_1 = odh_source('helium', Q_(4000,'ft**3'), 'gas', Q_(250, ureg.psig))
Buffer_tank_2 = odh_source('helium', Q_(4000,'ft**3'), 'gas', Q_(250, ureg.psig))
Buffer_tank_3 = odh_source('helium', Q_(4000,'ft**3'), 'gas', Q_(195, ureg.psig))
Buffer_tank_4 = odh_source('helium', Q_(4000,'ft**3'), 'gas', Q_(195, ureg.psig))
Buffer_tank_5 = odh_source('helium', Q_(4000,'ft**3'), 'gas', Q_(195, ureg.psig))
Buffer_tank_6 = odh_source('helium', Q_(4000,'ft**3'), 'gas', Q_(195, ureg.psig))
VTS_1 = odh_source('helium', Q_(1560, ureg.L), 'liquid')
VTS_2 = odh_source('helium', Q_(3256, ureg.L), 'liquid')
VTS_3 = odh_source('helium', Q_(3256, ureg.L), 'liquid')
VMTF = odh_source('helium', Q_(1450, ureg.L), 'liquid') #without displacer
Portable_500L = odh_source('helium', Q_(500, ureg.L), 'liquid')
Stand_4 = odh_source('helium', Q_(450, ureg.L), 'liquid') #LHC magnet and feed can, single phase
Suction_buffer = odh_source('helium', Q_(500, ureg.cubic_feet), 'gas', Q_(80, ureg.psig)) #suction buffer + return pipe
Quench_tank = odh_source('helium', 2*Q_(134, ureg.cubic_feet), 'gas', Q_(100, ureg.psig)) #two quench tanks that behave as one volume
Distribution_box = odh_source('helium', Q_(68, ureg.L), 'liquid') #distribution box He subcooler
MCTF = odh_source('helium', Q_(114, ureg.L), 'liquid') #calibration/research cryostat
LN2_portable_180L = odh_source('nitrogen', Q_(180, ureg.L), 'liquid')
LN2_storage_dewar = odh_source('nitrogen', Q_(10000, ureg.gallon), 'liquid')

Sources = [He_storage_dewar_liquid, He_storage_dewar_gas, Buffer_tank_1, Buffer_tank_2, Buffer_tank_3, Buffer_tank_4, Buffer_tank_5, Buffer_tank_6, VTS_1, VTS_2, VTS_3, VMTF, Portable_500L, Stand_4, Suction_buffer, Quench_tank, Distribution_box, MCTF, LN2_portable_180L, LN2_storage_dewar]

He_storage_dewar_liquid.leak({'mode':'fluid line', 'Pipe':{'D_nom':0.5, 'SCH':5}, 'N_lines':11, 'Fluid_data':{'P': 12*ureg.psig, 'T':4.2*ureg.K}}) #Transfer lines going between Cold Box and Dewar, and Distribution Box
He_storage_dewar_gas.leak({'mode':'gas line', 'Pipe':{'D_nom':4, 'SCH':5, 'L':35*ureg.ft}, 'N_welds':10, 'Fluid_data':{'P': 285*ureg.psig, 'T':300*ureg.K}, 'max_flow':907*ureg('ft^3/min')}) #Supply to Cold Box; Max flow is taken from Mycom room ODH analysis as Mycom compressor capacity
He_storage_dewar_gas.leak({'mode':'gas line', 'Pipe':{'D_nom':6, 'SCH':5, 'L':35*ureg.ft}, 'N_welds':10, 'Fluid_data':{'P': 35*ureg.psig, 'T':300*ureg.K}, 'max_flow':0*ureg('ft^3/min')}) #Mid stage; Max flow is taken from Mycom room ODH analysis as Mycom compressor capacity
He_storage_dewar_gas.leak({'mode':'gas line', 'Pipe':{'D_nom':8, 'SCH':5, 'L':35*ureg.ft}, 'N_welds':10, 'Fluid_data':{'P': 16.5*ureg.psig, 'T':300*ureg.K}, 'max_flow':6044*ureg('ft^3/min')}) #Return (suction); Max flow is taken from Mycom room ODH analysis as Total supply rate to suction piping


def print_leaks (odh_source):
    for key in sorted(odh_source.Leaks.keys()):
        print (key)
        print (odh_source.Leaks[key])
        print ()


for source in Sources:
    #if source.safe():
    #    print (source.volume)
    print (source.volume)
    print()
    print_leaks (source)
    print('\n'*2)














    
if __name__ == "__main__":
    pass
#    Fluid = {'fluid':'air', 'P':20*101325*ureg('Pa'), 'T':Q_(315,ureg.K)}
#    d_leak = 1*ureg.mm
    #print (leak_flow())
#    print (Q_(14.7,ureg('psig')).to(ureg.psi))
#    print (Q_(0,ureg('psig')).to(ureg.psi))

