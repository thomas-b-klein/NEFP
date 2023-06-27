from utility import *
from ifr import ISENTROPIC as ISEN
import numpy as np
import matplotlib.pyplot as plt

#Normal Shock Relations

class NORMAL():

    def __init__(self, mach, gamma = 1.4):
        self.M = np.asanyarray(mach)
        self.gamma = gamma
    
    @property
    def u2u1(self): #F1 = u2/u1
        u2_u1 = ( 2 + (self.gamma - 1) * self.M**2 ) / ( (self.gamma + 1) * self.M**2 )
        return u2_u1
    F1 = u2u1
    
    @property
    def r2r1(self): #F2 = rho2/rho1
        rho2_rho1 = 1 / self.F1
        return rho2_rho1
    F2 = r2r1

    @property
    def t2t1(self): #F3 = T2/T1
        T2_T1 = 1 + (self.gamma - 1) / 2 * self.M**2 - (self.gamma - 1) / 2 * self.M**2 * self.F1**2
        return T2_T1
    F3 = t2t1
    
    @property
    def p2p1(self): #F4 = p2/p1
        p2_p1 = self.F2 * self.F3
        return p2_p1
    F4 = p2p1

    @property
    def M2(self): #F5 = M2
        M2 = self.F1 * self.M / self.F3
        return M2
    F5 = M2

    @property
    def ds(self): #F6 = (s2-s1)/Cv
        ds = np.log( self.F3 * (1 / self.F2)**(self.gamma - 1) )
        return ds
    F6 = ds

    @property
    def r02r01(self): #F7 = rho02/rho01
        rho02_rho01 = np.exp( self.F6 )**(1 / (1 - self.gamma))
        return rho02_rho01
    F7 = r02r01

    @property
    def p02p01(self): #F8 = p02/p01
        p02_p01 = self.F7
        return p02_p01
    F8 = p02p01

    @property
    def a2sa1s(self): #F9 = A2*/A1*
        A2sA1s = 1 / self.F8
        return A2sA1s
    F9 = a2sa1s

    @property
    def p1p02(self): 
        p1_p01 = ISEN(self.gamma).pp0(self.M)
        p02_p01 = self.p02p01
        p1_p02 = p1_p01 / p02_p01
        return p1_p02