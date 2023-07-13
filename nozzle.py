from utility import *
from plotting import *
from ifr import ISENTROPIC as ISEN
from nsr import NORMAL as NORM
# from osr import OBLIQUE as OBLQ

import numpy as np

class GEOMETRY():
    def __init__(self, aeas, gamma = 1.4):
        if aeas < 1.0: raise ValueError('A/A* must be greater or equal to 1.') 

        self.aeas = aeas
        self.gamma = gamma

        self.exit_sub = ISEN.from_aas(aeas, gamma = gamma, subsonic = True)
        self.exit_sup = ISEN.from_aas(aeas, gamma = gamma, subsonic = False)

    @property
    def subsonic(self):
        return (1,self.choked)

    @property 
    def choked(self):
        return float(self.exit_sub.pp0)

    @property
    def shock_in_nozzle(self):
        return (self.choked, self.shock_at_exit)

    @property
    def shock_at_exit(self):
        p2p1 = NORM(self.exit_sup.M, gamma = self.gamma).p2p1
        return float(p2p1 * self.design)
    
    @property
    def over_expanded(self):
        return (self.shock_at_exit, self.design)

    @property
    def design(self):
        return float(self.exit_sup.pp0)

    @property
    def under_expanded(self):
        return (self.design, 0)

    def ranges(self):
        return [self.subsonic, self.shock_in_nozzle, self.over_expanded, self.under_expanded]
    
    def point(self):
        return [self.choked, self.shock_at_exit, self.design]



class NOZZLE():
    def __init__(self, aeas, pbpc, gamma = 1.4):
        if pbpc <= 0.0: raise ValueError('pb/pc value must be greater than 0.')
        elif pbpc > 1: raise ValueError('Reversed flow: pb/pc value must be less than 1.')
        elif pbpc == 1: raise ValueError('No flow: pb/pc value must be less than 1.')

        self.pbpc = pbpc
        self.aeas = aeas

        pbpc_ranges = GEOMETRY(aeas, gamma=gamma).ranges()
        pbpc_points = GEOMETRY(aeas, gamma=gamma).points()
        functions = ([NOZZLE.choked, NOZZLE.shock_at_exit, NOZZLE.design],
                    [NOZZLE.subsonic, NOZZLE.shock_in_nozzle, NOZZLE.over_expanded, NOZZLE.under_expanded])
        for func, point in zip(functions[0], pbpc_points):
            if pbpc == point: func(self); break
        for func, range in zip(functions[1], pbpc_ranges):
            if within(pbpc, range): func(self); break

    def subsonic(self):
        print('subsonic')

    def shock_in_nozzle(self):
        print('sin')
    
    def over_expanded(self):
        print('over')

    def under_expanded(self):
        print('under')

    def choked(self):
        print('choked')
    
    def shock_at_exit(self):
        print('shock at exit')

    def design(self):
        print('design')

# print(GEOMETRY(5).ranges())
# print(NOZZLE(5,0.995))