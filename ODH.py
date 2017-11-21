#!python3
#Defining some functions that will be of use soon

import numpy as np
import matplotlib.pyplot as plt
from pint import UnitRegistry

ureg = UnitRegistry()
ureg.auto_reduce_dimensions = True

def conc_vent (V, R, Q, t):
    #V - volume of the confined space (ft3 or m3)
    #R - spill rate into confined space (scfm or m3/s)
    #Q = ventilation rate of fan(s), (cfm or m3/s); positive value corresponds to blowing air into the confined space, negative - drawing contaminated air outside
    #t = time, (minutes or seconds) beginning of release is at t=0
    #C - oxygen concentration in confined space
    #Case B
    if Q > 0*ureg('ft^3/min'):
            C = 0.21/(Q+R)*(Q+R*np.e**-((Q+R)/V*t))
    elif abs(Q) <= abs(R):
        C = 0.21*np.e**-(R/V*t)
    elif abs(Q) > abs(R):
            C = 0.21*(1-R/abs(Q)*(1-np.e**-(abs(Q)*t/V)))
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


V = 1000*ureg('ft^3')
R = 10*ureg('ft^3/min')
Q = -20*ureg('ft^3/min')
t_e = 5*10**3*ureg('sec')
time = np.linspace(0,10**4, 10**3)*ureg('s')
time1 = [t for t in time if (t-t_e)<0*ureg('s')]
time2 = [t for t in time if (t-t_e)>0*ureg('s')]
time2.insert(0,time1[-1])
concentration1 = list(map (lambda t: conc_vent(V, R, Q, t), time1))
C_e = conc_vent(V, R, Q, t_e)
concentration2 = list(map(lambda t: conc_after (V, C_e, Q, t, t_e), time2))

plt.plot([t.magnitude for t in time1], concentration1, [t.magnitude for t in time2], concentration2)
a = 2*ureg("m")/(1*ureg('ft'))
# a.ito_base_units()
print (a)

# plt.plot(time1, concentration1, time2, concentration2)
# plt.show()




# V = 1000*ureg('ft^3')
# R = 10*ureg('ft^3/min')
# Q1 = -20*ureg('ft^3/min')
# Q2 = 0*ureg('ft^3/min')
# Q3 = 20*ureg('ft^3/min')

# time = np.linspace(0,10**5,100)*ureg('s')
# # print(time)
# concentration1 = list(map (lambda t: conc_vent(V, R, Q1, t), time))
# concentration2 = list(map (lambda t: conc_vent(V, R, Q2, t), time))
# concentration3 = list(map (lambda t: conc_vent(V, R, Q3, t), time))
# # print ('{:P}'.format(concentration[3]))
# print (conc_final(V,R,Q3))
# plt.plot(time, concentration1, time, concentration2, time, concentration3)
# plt.show()
# # # print (vent_out_little(1000, 100, 20, 0))
