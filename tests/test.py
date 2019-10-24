import ODH_analysis as odh
import heat_transfer as ht
from pprint import pprint
Test = odh.Source('Storage dewar', ht.ThermState('helium'), ht.Q_(3390000, ht.ureg.cubic_feet)) #blowdown from 12 psig to 1 atmosphere, estimated by R. Rabehl, TID-N-3A, p. 7
Test.Fluid.update('T', ht.T_NTP, 'P', ht.P_NTP)
Test.gas_pipe_failure(ht.piping.Pipe(1))
pprint(Test.Leaks)

#Test = Source('Test', 'helium', Q_(33900, ureg.cubic_feet), 'vapor', Q_(0, ureg.psig)) #blowdown from 12 psig to 1 atmosphere, estimated by R. Rabehl, TID-N-3A, p. 7
#Test_2 = Test+Test+Test
#Source.delete()
#He_storage_dewar_gas.leak({'mode':'gas line', 'Pipe':Tube(6*ureg.inch, 0.15*ureg.inch, 350000*ureg.ft), 'N_welds':10, 'Fluid_data':{'P': 95*ureg.psig, 'T':300*ureg.K}, 'max_flow':108*ureg('g/s')}) #Mid stage; Max flow is taken as midstage return from cryoplant
#Test_pipe = Tube(1*ureg.inch, 0.035*ureg.inch, 20*ureg.ft)
#Test_piping = Piping({'fluid':'helium', 'P':12*ureg.psig, 'T':5*ureg.K}, Test_pipe)
#print(Test_piping.m_dot().to(ureg.g/ureg.s))
#print(Test_piping.dP(607*ureg('g/s')))
#Test_period = 1*ureg('month') #IB1 fan test period 
#Test_period.ito(ureg.hr)
#Mean_repair_time = 3*ureg('days') #IB1 fan/louver average repair time: most delay is caused by response time, it has recently improved from about 1 week to 1 day. The average value of 3 days is assumed
#l_fan = 9e-6/ureg.hr #Fan failure rate
#l_louv = 3e-7/ureg.hr #Louver failure rate
#l_vent = l_fan+l_louv #At IB1 louver and fan are installed in series; failure of either one results in no venting
#Q_fan = 4000*ureg.ft**3/ureg.min #Flowrate of 4 ceiling fans is >= 4000 CFM; positive value refers to blowing into the building: FALSE IN CASE OF IB1!
#N_fans = 4
#
    #A = Q_(15360, ureg.square_feet) #IB1 floor area
    #IB1_air = Volume('IB1 air', ['helium', 'nitrogen'], A*19.8*ureg.feet) #All IB1 air
#
    #IB1_air.fan_fail(Test_period, l_vent, Q_fan, N_fans, Mean_repair_time)
    #IB1_air.odh()
    #print_result(IB1_air)
#
    #Q = -16000*ureg('ft^3/min')
    #R =7750*ureg('ft^3/min')
    #tau = 266300*ureg.cubic_feet/R
    #H = 9.8*ureg.ft
    #V = H*A
    #C = conc_vent(V, R, Q, tau)
    #Test_vol = Volume ('Test', [],V)
    #print ('Volume: {:.0f}'.format(V.to(ureg.ft**3)))
    #print ('Oxygen concentration: {:.0%}'.format(C))
    #print ('Oxygen pressure: {:.0f}'.format(C/0.21*160)) 
    #print ('Fatality factor: {:.0g}'.format(Test_vol.fatality_prob(C)))
    #print (prob_m_of_n(1, 1, 1*ureg.hr, 9e-6/ureg.hr))
