#!python3
#Defining base classes for ODH analysis

from pint import UnitRegistry

ureg = UnitRegistry()
ureg.auto_reduce_dimensions = True


class odh_source:
    """Define the possible source of inert gas"""
    def __init__ (self, Volume):
        self.volume = Volume

#class odh_failure(odh_source):
#    """Define failure mode for existing source of inert gas"""
#    def __init__ (self, ODH_source, 


def leak_flow (d1, Fluid_data = {'fluid':'air', 'P':101325*ureg('Pa'), 'T':38*ureg('degC')}):
    A = math.pi*d1**2/4
    Y = 1 #conservative value; from Crane TP-410 A-21
    C = 0.7 #conservative value; from Crane TP-410 A-20 
    P_atm = 101325*ureg('Pa')
    DeltaP = 
    q = Y*C*A*(2*DeltaP/rho)**0.5
