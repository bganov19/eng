from math import exp, sqrt, pi
import scipy.optimize as sp
import numpy as np
from constants import *


class PPV:
    N_T = 5*10 ** 17
    N_Sh = 6*10 ** 17
    u_ph = 5 * 10 ** 12
    sigma_h = 10 ** (-13)
    sigma_e = 10 ** (-16)
    E_HC = 0.35
    E_HE = 0.85
    E_EC = 0.1

    def __init__(self, temp=300, light=0):
        self.T = temp
        self.n = light

    @property
    def v_th_e(self):
        return 100 * sqrt(3 * kb * self.T / m_e)  # cm/s

    @property
    def v_th_h(self):
        return 100 * sqrt(3 * kb * self.T / m_h)  # cm/s

    @property
    def N_V(self):
        return 2 * (2 * pi * m_h * kb * self.T / (h ** 2)) ** (3 / 2) / 10 ** 6

    def p(self, f):
        return self.N_Sh + self.N_T * (2 * f - 1) + self.n

    def occupation(self, f, vte, vth, nv):
        return (-self.u_ph**(-1)*(vth*self.sigma_h*self.p(f))**2 * exp(-self.E_HC*eV/(kb*self.T)) * f +\
                                self.u_ph**(-1)*(vth*self.sigma_h*nv)**2 * exp(-self.E_HE*eV/(kb*self.T)) * (1-f) + \
                                self.u_ph**(-1)*vth*self.sigma_h*self.n*vte*self.sigma_e*nv * exp(-self.E_EC*eV/(kb*self.T)) * (1-f))/10**12

    def simple_occupation(self, f): return (1+((self.p(f)/self.N_V)**2)*exp((self.E_HE-self.E_HC)*eV/(kb*self.T)))**(-1) - f

    def equilibrium(self):
        buff = self.n
        self.n = 0
        oc = float(sp.fsolve(self.simple_occupation, 0))
        self.n = buff
        return oc

    def rk4_step(self, f, dt, vte, vth, nv):
        k1 = self.occupation(f, vte, vth, nv)
        k2 = self.occupation(f+k1*dt/2, vte, vth, nv)
        k3 = self.occupation(f+k2*dt/2, vte, vth, nv)
        k4 = self.occupation(f+k3, vte, vth, nv)
        return (k1 + 2*k2 + 2*k3 + k4)/6*10**12

    def ls_simulation(self, dt, lstime, T0=None):

        if T0 is not None:
            buff = self.T
            self.T = T0
            c0 = self.equilibrium()
            self.T = buff
        else:
            c0 = self.equilibrium()

        print(self.equilibrium())
        time = np.linspace(0, lstime, int(lstime/dt))
        conc = [c0]
        h0 = self.p(c0) - self.n
        holes = [h0]

        vte = self.v_th_e
        vth = self.v_th_h
        nv = self.N_V

        for t in time:
            o = conc[-1] + self.rk4_step(conc[-1], dt, vte, vth, nv) * dt
            conc.append(o)
            holes.append(self.p(o) - self.n)

        return holes

    def ls_wide_simulation(self, dt, T0, Ti, Tf, m, lstime):
        temperatures = np.linspace(Ti, Tf, m)
        holes = []

        for Te in temperatures:
            self.T = Te
            hol = self.ls_simulation(dt, lstime, T0)
            holes.append(hol)
            # print(self.n, Te, " {0:1.2}".format(hol))

        return holes

    def tHE(self):
        return self.u_ph**(-1)*(self.sigma_h*self.v_th_h*self.N_V)**2*exp(-self.E_HE*eV/(kb*self.T))

    def tEC(self):
        return self.u_ph**(-1)*self.sigma_e*self.v_th_e*self.n*self.sigma_h*self.v_th_h*self.N_V*exp(-self.E_EC*eV/(kb*self.T))

    def tHC(self, f):
        return self.u_ph**(-1)*(self.sigma_h*self.v_th_h*self.p(f))**2*exp(-self.E_HC*eV/(kb*self.T))

    def occupation2(self, f):
        return (self.tHE() + self.tEC()) / (self.tHE() + self.tEC() + self.tHC(f)) - f

    def eq_state(self, temp, light):
        self.T = temp
        self.n = light
        eq = float(sp.fsolve(self.occupation2, 0))
        return eq

