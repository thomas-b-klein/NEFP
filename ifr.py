from utility import *
import numpy as np
from scipy.optimize import newton

#Isentropic Flow Relations

class ISENTROPIC():

    def __init__(self, gamma = 1.4):
        if gamma <= 0: raise ValueError('Gamma must be greater than 1.')
        self.gamma = gamma
    
    def m_tt0(self, tt0):
        if tt0 <= 0.0 or tt0 >= 1.0: raise ValueError('T/T0 must be between 0 and 1.') 
        tt0 = np.asanyarray(tt0)
        M = np.sqrt( (tt0**(-1) - 1) * 2/(self.gamma-1) )
        return M

    def m_rr0(self, rr0):
        if rr0 <= 0.0 or rr0 >= 1.0: raise ValueError('rho/rho0 must be between 0 and 1.') 
        rr0 = np.asanyarray(rr0)
        M = np.sqrt( (rr0**(-(self.gamma-1)) - 1) * 2/(self.gamma-1) )
        return M

    def m_pp0(self, pp0):
        if pp0 <= 0.0 or pp0 >= 1.0: raise ValueError('p/p0 must be between 0 and 1.') 
        pp0 = np.asanyarray(pp0)
        M = np.sqrt( (pp0**(-(self.gamma-1)/self.gamma) - 1) * 2/(self.gamma-1) )
        return M

    def m_aas(self, aas, isen_object = None, subsonic = False):
        if not isen_object: isen_object = ISENTROPIC(gamma=1.4)
        if aas <= 0.0 or aas < 1.0: raise ValueError('A/A* must be between 0 and 1.') 
        elif aas == 1.0: M = 1.0
        else: 
            est = 0.5 if subsonic else 2.0
            aas = np.asanyarray(aas)
            eq = make_implicit(isen_object.aas)
            M = newton(eq, est, args= (aas,))    
        return M

    def m_mu(self, mu, rad = True):
        if not rad: mu = np.radians(mu)
        if mu < 0.0 or mu > np.pi/2:
            raise ValueError('The Mach angle must be between (0, pi/2) radians.')
        mu = np.asanyarray(mu)
        M = (np.sin(mu))**(-1)
        return M

    def m_pm(self, pm, isen_object = None, rad = True):
        if not rad: pm = np.radians(pm)
        if not isen_object: isen_object = ISENTROPIC(gamma=1.4)
        pm = np.asanyarray(pm)
        pm_max = (np.pi/2 * ( np.sqrt( (self.gamma+1)/(self.gamma-1) ) - 1 ))
        if pm <= 0.0 or pm > pm_max:
            raise ValueError('The Prandtl-Meyer angle must be between (0, {}) radians.'.format(pm_max))
        eq = make_implicit(ISENTROPIC().pm)
        M = newton(eq, 2.0, args = (pm,))
        return M

    def tt0(self, M = 1.0):
        if M <= 0.0: raise ValueError('Mach must be greater than zero.')
        M = np.asanyarray(M)
        tt0 =  (1 + (self.gamma-1)/2 * M**2) ** (-1)
        return tt0

    def rr0(self, M = 1.0):
        if M <= 0.0: raise ValueError('Mach must be greater than zero.')
        M = np.asanyarray(M)
        rr0 = (1 + (self.gamma-1)/2 * M**2) ** (-1/(self.gamma-1))
        return rr0
    
    def pp0(self, M = 1.0):
        if M <= 0.0: raise ValueError('Mach must be greater than zero.')
        M = np.asanyarray(M)
        pp0 = (1 + (self.gamma-1)/2 * M**2) ** (-self.gamma/(self.gamma-1))
        return pp0

    def aas(self, M = 1.0):
        if M <= 0.0: raise ValueError('Mach must be greater than zero.')
        M = np.asanyarray(M)
        aas = (1/M) * ((2/(self.gamma+1)) * (1 + ((self.gamma-1)/2) * M**2)) ** ((self.gamma+1) / (2*(self.gamma-1)))
        return aas

    def mu(self, M = 1.0, rad = True):
        if M <= 0.0: raise ValueError('Mach must be greater than zero.')
        M = np.asanyarray(M)
        mu = np.arcsin( 1/M )
        return mu if rad else np.degrees(mu)

    def pm(self, M = 1.0, rad = True):
        if M <= 0.0: raise ValueError('Mach must be greater than zero.')
        M = np.asanyarray(M)
        gamma_ratio = (self.gamma+1)/(self.gamma-1)
        pm = np.sqrt(gamma_ratio) * np.arctan( np.sqrt( (M**2-1)/gamma_ratio ) ) - np.arctan( np.sqrt( M**2-1 ) )
        return pm if rad else np.degrees(pm)