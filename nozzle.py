from utility import *
from ifr import ISENTROPIC as ISEN
from nsr import NORMAL as NORM
# from osr import OBLIQUE as OBLQ

import matplotlib.pyplot as plt
import numpy as np

class GEOMETRY():
    def __init__(self, aeas, gamma = 1.4):
        self.aeas = aeas
        self.gamma = gamma
        self.exit_sub = ISEN.from_aas(aeas, gamma = gamma, subsonic = True)
        self.exit_sup = ISEN.from_aas(aeas, gamma = gamma, subsonic = False)

    @property
    def subsonic(self):
        return (1,self.choked)

    @property 
    def choked(self):
        return self.exit_sub.pp0

    @property
    def shock_in_nozzle(self):
        return (self.choked, self.shock_at_exit)

    @property
    def shock_at_exit(self):
        p2p1 = NORM(self.exit_sup.M, gamma = self.gamma).p2p1
        return p2p1 * self.design
    
    @property
    def over_expanded(self):
        return (self.shock_at_exit, self.design)

    @property
    def design(self):
        return self.exit_sup.pp0

    @property
    def under_expanded(self):
        return (self.design, np.inf)

    def find_condition(self, pbpc):
        if pbpc >= 1.0 or pbpc <= 0.0: raise ValueError('pb/pc cannot be greater than 1 or less than 0.')
        if within(pbpc, self.subsonic): return (0,'Subsonic')
        elif pbpc == self.choked: return (1, 'Choked')
        elif within(pbpc, self.shock_in_nozzle): return (2, 'Shock in Nozzle')
        elif pbpc == self.shock_at_exit: return (3, 'Shock at Exit')
        elif within(pbpc, self.over_expanded): return (4, 'Over-Expanded')
        elif pbpc == self.design: return (5, 'Design')
        elif within(pbpc, self.under_expanded): return (6, 'Under-Expanded')
        else: 
            categories = [self.subsonic, self.shock_in_nozzle, self.over_expanded, self.under_expanded]
            raise ValueError('pb/pc does not fit into any of the categories:\n{}'.format(categories))
