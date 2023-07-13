from utility import *
import numpy as np
from scipy.optimize import newton

#Functions for Optimization (Newton's Method)

def aas_func(M, gamma):
    return (1/M) * ((2/(gamma+1)) * (1 + ((gamma-1)/2) * M**2)) ** ((gamma+1) / (2*(gamma-1)))

def pm_func(M, gamma):
    return np.sqrt((gamma+1)/(gamma-1)) * np.arctan( np.sqrt( (M**2-1)/((gamma+1)/(gamma-1)) ) ) - np.arctan( np.sqrt( M**2-1 ) )

#Isentropic Flow Class

class ISENTROPIC():

    #all angles input and output are in radians unless specifically stated
    rad = True 

    def __init__(self, M, gamma = 1.4):
        if gamma <= 0: raise ValueError('Gamma must be greater than 1.')
        self.gamma = gamma
        
        if type(M) == int or type(M) == float: self.M = np.array([M])
        else: self.M = np.array(M)
        if any(i <= 0.0 for i in self.M): raise ValueError('Mach must be greater than zero.')
        
    @property    
    def tt0(self): 
        tt0 =  (1 + (self.gamma-1)/2 * self.M**2) ** (-1)
        return tt0

    @property
    def rr0(self):
        rr0 = (1 + (self.gamma-1)/2 * self.M**2) ** (-1/(self.gamma-1))
        return rr0

    @property    
    def pp0(self):
        pp0 = (1 + (self.gamma-1)/2 * self.M**2) ** (-self.gamma/(self.gamma-1))
        return pp0

    @property
    def aas(self):
        aas = aas_func(self.M, self.gamma)
        return aas

    @property
    def mu(self, rad = True):
        mu = np.arcsin( 1/self.M )
        return mu if self.rad else np.degrees(mu)

    @property
    def pm(self, rad = True):
        pm = pm_func(self.M, self.gamma)
        return pm if self.rad else np.degrees(pm)

    @classmethod
    def from_tt0(cls, tt0, gamma = 1.4): #m_tt0
        if type(tt0) == int or type(tt0) == float: tt0 = np.array([tt0])
        else: tt0 = np.array(tt0)
        if any(i <= 0.0 or i >= 1.0 for i in tt0): raise ValueError('T/T0 must be between 0 and 1.') 
        tt0 = np.asanyarray(tt0)
        M = np.sqrt( (tt0**(-1) - 1) * 2/(gamma-1) )
        return cls(M, gamma)

    @classmethod
    def from_rr0(cls, rr0, gamma = 1.4): #m_rr0
        if type(rr0) == int or type(rr0) == float: rr0 = np.array([rr0])
        else: rr0 = np.array(rr0)
        if any(i <= 0.0 or i >= 1.0 for i in rr0): raise ValueError('rho/rho0 must be between 0 and 1.') 
        M = np.sqrt( (rr0**(-(gamma-1)) - 1) * 2/(gamma-1) )
        return cls(M, gamma)

    @classmethod
    def from_pp0(cls, pp0, gamma = 1.4): #m_pp0
        if type(pp0) == int or type(pp0) == float: pp0 = np.array([pp0])
        else: pp0 = np.array(pp0)
        if any(i <= 0.0 or i >= 1.0 for i in pp0): raise ValueError('p/p0 must be between 0 and 1.') 
        M = np.sqrt( (pp0**(-(gamma-1)/gamma) - 1) * 2/(gamma-1) )
        return cls(M, gamma)

    @classmethod
    def from_aas(cls, aas, gamma = 1.4, subsonic = False): #m_aas
        if type(aas) == int or type(aas) == float: aas = np.array([aas])
        else: aas = np.array(aas)
        if any(i <= 0.0 or i < 1.0 for i in aas): raise ValueError('A/A* must be greater or equal to 1.') 
        else: 
            est = 0.1 if subsonic else 2.0
            eq = make_implicit(aas_func)
            M = [newton(eq, est, args= (aas_i, gamma)) for aas_i in aas]   
        return cls(M, gamma)

    @classmethod
    def from_mu(cls, mu, gamma = 1.4, rad = True): #m_mu
        if not rad: mu = np.radians(mu)
        if type(mu) == int or type(mu) == float: mu = np.array([mu])
        else: mu = np.array(mu)
        if any(i < 0.0 or i > np.pi/2 for i in mu):
            raise ValueError('The Mach angle must be between (0, pi/2) radians.')
        M = (np.sin(mu))**(-1)
        return cls(M, gamma)

    @classmethod
    def from_pm(cls, pm, gamma = 1.4, rad = True): #m_pm
        if not rad: pm = np.radians(pm)
        if type(pm) == int or type(pm) == float: pm = np.array([pm])
        else: pm = np.array(pm)
        pm_max = (np.pi/2 * ( np.sqrt( (gamma+1)/(gamma-1) ) - 1 ))
        if any(i <= 0.0 or i > pm_max for i in pm):
            raise ValueError('The Prandtl-Meyer angle must be between (0, {}) radians.'.format(pm_max))
        eq = make_implicit(pm_func)
        M = [newton(eq, 2.0, args = (pm_i, gamma)) for pm_i in pm]
        return cls(M, gamma)