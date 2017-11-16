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

    def leak (self, failure_mode = {'mode': 'piping', 'Pipe':{'L':10*ureg.m, 'D_nom':3*ureg.inch, 'SCH':5}, 'fluid': {'fluid':'air', 'P':2.33*ureg('bar'), 'T':Q_(40,ureg.degC)}, 'max_flow':0.01*ureg('m^3/s')}):
        '''Calculating leak and probability for different cases of equipment failure
        '''
        Leaks = {}
        if failure_mode['mode'] == 'piping':
            Pipe = failure_mode['Pipe']
            cause = 'small leak'
            Prob = 10**(-9)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
            q = leak_flow(failure_mode['fluid'], 10*ureg.mm**2)
            Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

            if (Pipe.get('ID') or Pipe.get('D_nom') or Pipe.get('OD')) > 2*ureg.inch:
                cause = 'large leak'
                Prob = 10**(-10)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
                q = leak_flow(failure_mode['fluid'], 1000*ureg.mm**2)
                Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

        cause = 'rupture'
        Prob = 3*10**(-11)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
        q = leak_flow(failure_mode['fluid'], pipe.Area(Pipe))
        Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

        cause = 'weld small leak'
        Prob = 3*10**(-11)/(ureg.m*ureg.hr)*Pipe['L'] #Probaility and flow will be recalculated for each cause using the same variable names
        q = leak_flow(failure_mode['fluid'], pipe.Area(Pipe))
        Leaks[cause] = (Prob, self.limit_flow(q, failure_mode))

        return Leaks

    def limit_flow(self, flow, failure_mode):
        '''
        Leak through the openning sometimes cannot be larger than a certain number, e.g. a compressor throughput. This function limits flow to this value if it exists
        '''
        max_flow = failure_mode.get('max_flow', float('inf')*ureg.m**3/ureg.s) #If max_flow value is not specified infinity is assumed
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




def leak_flow (Fluid_data = {'fluid':'air', 'P':2.33*ureg('bar'), 'T':Q_(40,ureg.degC)}, A = 10*ureg('mm**2')):
    d = (4*A/math.pi)**0.5
    Y = 1 #conservative value; from Crane TP-410 A-21
    C = 0.7 #conservative value; from Crane TP-410 A-20 
    (x_a, M_a) = ht.rp_init({'fluid':'air'})
    (fluid, T_fluid, P_fluid) = ht.unpack_fluid(Fluid_data)
    (x, M, D_fluid) = ht.rp_init(Fluid_data)
    S_g = M/M_a #Specific gravity
    rho = D_fluid*M
    k = ht.gamma(Fluid_data) #adiabatic coefficient
    rc = (2/(k+1))**(k/(k-1)) #Critical pressure drop; Note: according to Crane TP-410 is depndent on the hydraulic resistance of the flow path
    if P_MSC > P_fluid*rc: #subsonic flow
        DeltaP = P_fluid - P_MSC
    else: #Sonic flow
        DeltaP = P_fluid*(1-rc) #Crane TP-410, p 2-15
    q = orifice_flow(Y, d, C, S_g, DeltaP, rho)
    return q

def orifice_flow (Y, d, C, S_g, DeltaP, rho):
    '''Wrapper function for flow through orifice equation.
    Original formula in crane uses a non-dimensionless coefficient 0.0002864 and non-SI units for values in formulas so straightforward approach is complicated.
    '''
    d_1 = d.to(ureg.mm).magnitude #diameter in mm
    Deltap = DeltaP.to(ureg.bar).magnitude #Pressure drop in bar
    rho_1 = rho.to(ureg.kg/ureg.m**3).magnitude #Density of inflow in kg/m^3
    q = 0.0002864*Y*d_1**2*C/(S_g)*(Deltap*rho_1)**0.5 #Flow at MSC in m^3/s; Crane TP-410, Eq.3-22
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

for source in Sources:
    if source.safe():
        print (source.volume)
print (VTS_1.leak())








    
if __name__ == "__main__":
    pass
#    Fluid = {'fluid':'air', 'P':20*101325*ureg('Pa'), 'T':Q_(315,ureg.K)}
#    d_leak = 1*ureg.mm
    #print (leak_flow())
#    print (Q_(14.7,ureg('psig')).to(ureg.psi))
#    print (Q_(0,ureg('psig')).to(ureg.psi))

