from utility import *
from ifr import ISENTROPIC as ISEN
import numpy as np

#Oblique Shock Relations

class OBLIQUE():

    def __init__(self, M1, gamma=1.4):
        self.M1 = np.asanyarray(M1)
        self.gamma = gamma

