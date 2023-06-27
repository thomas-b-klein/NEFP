from ifr import ISENTROPIC as ISEN
from nsr import NORMAL as NORM
from osr import OBLIQUE as OBLQ

import matplotlib.pyplot as plt
import numpy as np

#NORMAL SHOCK RELATIONS
def normal_shock_relations():
    # u2/u1 vs. M
    M = np.arange(0.1,10,0.1)
    fig1 = plt.figure(1)
    plt.plot(M,NORM(M).F1)
    plt.axvline(x=1, color='k', linestyle='--'); plt.axhline(y=1, color='k', linestyle='--')
    plt.xlim(0,10); plt.ylim(0,2)
    plt.xlabel('Mach'); plt.ylabel('u_2/u_1')

    # rho2/rho1, T2/T1, p2/p1 vs. M
    M = np.arange(0.1,10,0.1)
    fig2 = plt.figure(2)
    plt.plot(M, NORM(M).r2r1, color='b')
    plt.plot(M, NORM(M).F3, color='r')
    plt.plot(M, NORM(M).F4, color='g')
    plt.axvline(x=1, color='k', linestyle='--'); plt.axhline(y=1, color='k', linestyle='--')
    plt.xlim(0,10); plt.ylim(0,20)
    plt.xlabel('Mach'); plt.legend(['rho_2/rho_1', 'T_2/T_1', 'p_2/p_1'])

    # M2 vs. M
    M = np.arange(0.1,10,0.1)
    fig3 = plt.figure(3)
    plt.plot(M,NORM(M).F5)
    plt.axvline(x=1, color='k', linestyle='--'); plt.axhline(y=1, color='k', linestyle='--')
    plt.xlim(0,10); plt.ylim(0,2)
    plt.xlabel('M1'); plt.ylabel('M2')

    # ds vs. M
    M = np.arange(0.1,10,0.1)
    fig4 = plt.figure(4)
    plt.plot(M,NORM(M).F6)
    plt.axvline(x=1, color='k', linestyle='--'); plt.axhline(y=0, color='k', linestyle='--')
    plt.xlim(0,10); plt.ylim(-0.5,2.5)
    plt.xlabel('Mach'); plt.ylabel('s_2-s_1')

    # # p02/p01 and rho02/rho01 vs. M
    M = np.arange(0.1,10,0.1)
    fig5 = plt.figure(5)
    plt.plot(M,NORM(M).F8)
    plt.axvline(x=1, color='k', linestyle='--'); plt.axhline(y=1, color='k', linestyle='--')
    plt.xlim(0,10); plt.ylim(0,1.5)
    plt.xlabel('Mach'); plt.ylabel('p0_2/p0_1 and rhp0_2/rho0_1')

    # A2*/A1* vs. M
    M = np.arange(0.1,10,0.1)
    fig6 = plt.figure(6)
    plt.plot(M,NORM(M).F9)
    plt.axvline(x=1, color='k', linestyle='--'); plt.axhline(y=1, color='k', linestyle='--')
    plt.xlim(0,10); plt.ylim(0,20)
    plt.xlabel('Mach'); plt.ylabel('A_2*/A_1*')

    plt.show()

def oblique_shock_relations():
    ...