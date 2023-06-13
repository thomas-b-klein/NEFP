from utility import *
from ifr import ISENTROPIC as ISEN
from nsr import NORMAL as NORM
import numpy as np

#Oblique Shock Relations

class OBLIQUE():

    #all angles input and output are in radians unless specifically stated
    rad = True

    def __init__(self, M1, delta, weak = True, gamma=1.4):
        self.M1 = np.array(M1)
        self.delta = delta
        self.weak = weak
        self.gamma = gamma

    @property
    def M2(self):
        M2 = self.M2n / np.sin(self.beta - self.delta)
        return M2

    @property
    def M1n(self):
        M1n = self.M1*np.sin(self.beta)
        return M1n

    @property
    def M2n(self):
        M2n = NORM(self.M1n, self.gamma).M2
        return M2n

    @property    
    def beta(self):
        #!FIX: change function to iterate Newton's method - slightly faster
        #wave angle of the flow due to the shockwave
        weak = 1 if self.weak else 0
        lam = np.sqrt( (self.M1**2-1)**2 - 3*(1+(self.gamma-1)*self.M1**2/2)*(1+(self.gamma+1)*self.M1**2/2)*np.tan(self.delta)**2)
        chi = ( (self.M1**2 - 1)**3 - 9*(1+(self.gamma-1)*self.M1**2/2)*(1+(self.gamma-1)*self.M1**2/2 + (self.gamma+1)*self.M1**4/4) * np.tan(self.delta)**2 ) / lam**3
        beta = np.arctan(((self.M1**2 - 1) + 2*lam*np.cos((4*np.pi*weak + np.arccos(chi))/3)) / (3*(1 + (self.gamma-1)*self.M1**2/2)*np.tan(self.delta)))
        return beta

    @property
    def p2p1(self):
        p2p1 = NORM(self.M1n, self.gamma).p2p1
        return p2p1

    @property    
    def r2r1(self):
        r2r1 = NORM(self.M1n, self.gamma).r2r1
        return r2r1

    @property
    def t2t1(self):
        t2t1 = NORM(self.M1n, self.gamma).t2t1
        return t2t1

    @property
    def p02p01(self):
        p02p01 = NORM(self.M1n, self.gamma).p02p01
        return p02p01

    @classmethod
    def from_M1n(cls, M1, M1n, gamma=1.4):
        beta = np.arcsin( M1n / M1 )
        return cls.from_beta(M1, beta, gamma)

    @classmethod
    def from_beta(cls, M1, beta, gamma=1.4):
        #!FIX: Find out how to determine whether a shock is strong or weak
        delta = np.arctan( 2*(M1**2*np.sin(beta)**2 - 1) / (np.tan(beta) * ( 2 + M1**2*(gamma + np.cos(2*beta)) ) ) )
        weak = None
        return cls(M1, delta, weak, gamma)
    