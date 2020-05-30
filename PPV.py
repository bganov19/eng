import scipy.optimize as sp
from constants import *
from dist_functions import *

np.random.seed(1)
bins = 1000


class PPV:
    N_T = 5*10 ** 15
    N_Sh = 5*10 ** 15
    u_ph = 5 * 10 ** 12
    sigma_h = 10 ** (-15)
    sigma_e = 10 ** (-16)
    E_HC = 0.35
    E_HE = 0.85
    E_EC = 0.1

    def __init__(self, temp=300, light=0):
        self.T = temp
        self.n = light
        # self.e_ec = np.random.normal(self.E_EC, 0.05, bins)
        # self.e_ec = np.array([self.E_EC]*bins)
        self.e_ec = inv_cum(np.random.uniform(0, 1, bins))
        self.f = np.zeros(bins)


    @property
    def v_th_e(self):
        return 100 * sqrt(3 * kb * self.T / m_e)  # cm/s

    @property
    def v_th_h(self):
        return 100 * sqrt(3 * kb * self.T / m_h)  # cm/s

    @property
    def N_V(self):
        return 2 * (2 * pi * m_h * kb * self.T / (h**2)) ** (3/2) / 10**6

    @property
    def p(self):
        return self.N_Sh + self.N_T*np.average(2*self.f - 1)

    def p0(self, f0):
        return self.N_Sh + self.N_T*(2*f0 - 1)

    def occupation(self, f, vte, vth, nv):
        return (-self.u_ph**(-1)*(vth*self.sigma_h*self.p)**2 * exp(-self.E_HC*eV/(kb*self.T)) * f +\
                                self.u_ph**(-1)*(vth*self.sigma_h*nv)**2 * exp(-self.E_HE*eV/(kb*self.T)) * (1-f) + \
                                self.u_ph**(-1)*vth*self.sigma_h*self.n*vte*self.sigma_e*nv * exp(-self.E_EC*eV/(kb*self.T)) * (1-f))/10**12

    def simple_occupation(self, f): return (1+((self.p0(f)/self.N_V)**2)*exp((self.E_HE-self.E_HC)*eV/(kb*self.T)))**(-1) - f

    def matrix_occupation(self, f, vte, vth, nv):
        diff = []
        for fi, ei in zip(f, self.e_ec):
            self.E_EC = ei
            diff.append(self.occupation(fi, vte, vth, nv))
        return np.array(diff)

    def equilibrium(self):
        buff = self.n
        self.n = 0
        diff = [float(sp.fsolve(self.simple_occupation, 0))] * bins
        self.n = buff
        return np.array(diff)

    def rk4_step(self, dt, vte, vth, nv):
        k1 = self.matrix_occupation(self.f, vte, vth, nv)
        k2 = self.matrix_occupation(self.f+k1*dt/2, vte, vth, nv)
        k3 = self.matrix_occupation(self.f+k2*dt/2, vte, vth, nv)
        k4 = self.matrix_occupation(self.f+k3, vte, vth, nv)
        return (k1 + 2*k2 + 2*k3 + k4)/6*10**12

    def ls_simulation(self, lstime, dt=None, T0=None, cut=True):

        if dt is None:
            # dt = min(10**(floor(-log10(self.tEC0(np.min(self.e_ec))))+1), lstime/100)
            dt = min((self.tEC0(np.min(self.e_ec))**-1), lstime / 100)
            print("dt = ", dt)

        if T0 is not None:
            buff = self.T
            self.T = T0
            self.f = self.equilibrium()
            self.T = buff
        else:
            self.f = self.equilibrium()

        conc = [self.f]
        h0 = self.p - self.n
        holes = [h0]

        vte = self.v_th_e
        vth = self.v_th_h
        nv = self.N_V

        p_ss = self.eq_state()

        for t in range(int(lstime/dt)):
            if holes[-1] < 0.999*p_ss or not cut:
                self.f = conc[-1] + self.rk4_step(dt, vte, vth, nv) * dt
                conc.append(self.f)
                holes.append(self.p - self.n)
            else:
                holes.append(0.999*p_ss)


        print (self.T, " K Done!")
        return holes

    def short_ls_simulation(self, lstime, dt=None, T0=None, cut=True):

        if dt is None:
            # dt = min(10**(floor(-log10(self.tEC0(np.min(self.e_ec))))+1), lstime/100)
            dt = min((self.tEC0(np.min(self.e_ec))**-1), lstime / 100)
            print("dt = ", dt)

        if T0 is not None:
            buff = self.T
            self.T = T0
            self.f = self.equilibrium()
            self.T = buff
        else:
            self.f = self.equilibrium()

        current_holes = self.p - self.n

        vte = self.v_th_e
        vth = self.v_th_h
        nv = self.N_V

        p_ss = self.eq_state()

        for t in range(int(lstime/dt)):
            if current_holes < 0.999*p_ss or not cut:
                self.f += self.rk4_step(dt, vte, vth, nv) * dt
                current_holes = self.p - self.n
            else:
                break

        print (self.T, " K Done!")
        return current_holes

    def ls_wide_simulation(self, T0, Ti, Tf, m, lstime, dt=None):
        temperatures = np.linspace(Ti, Tf, m)
        holes = []

        for Te in temperatures:
            self.T = Te
            hol = self.ls_simulation(lstime, dt, T0)
            holes.append(hol)
            # print(self.n, Te, " {0:1.2}".format(hol))

        return holes

    def short_ls_wide_simulation(self, T0, Ti, Tf, m, lstime, dt=None):
        temperatures = np.linspace(Ti, Tf, m)
        holes = []

        for Te in temperatures:
            self.T = Te
            holes.append(self.short_ls_simulation(lstime, dt, T0))

        return holes

    def tHE(self):
        return self.u_ph**(-1)*(self.sigma_h*self.v_th_h*self.N_V)**2*exp(-self.E_HE*eV/(kb*self.T))

    def tEC(self):
        return self.u_ph**(-1)*self.sigma_e*self.v_th_e*self.n*self.sigma_h*self.v_th_h*self.N_V*exp(-self.E_EC*eV/(kb*self.T))

    def tEC0(self, e_ec0):
        return self.u_ph**(-1)*self.sigma_e*self.v_th_e*self.n*self.sigma_h*self.v_th_h*self.N_V*exp(-e_ec0*eV/(kb*self.T))

    def tHC(self, p1):
        return self.u_ph**(-1)*(self.sigma_h*self.v_th_h*p1)**2*exp(-self.E_HC*eV/(kb*self.T))

    # def occupation2(self, fi):
    #     return (self.tHE() + self.tEC()) / (self.tHE() + self.tEC() + self.tHC(999)) - fi
    #
    # def occupation3(self, p1):
    #     return (self.tHE() + self.tEC()) / (self.tHE() + self.tEC() + self.tHC(p1)) + (self.N_Sh-self.N_T-p1)/(2*self.N_T)

    def occupation4(self, p1):
        t0 = 0
        for ei in self.e_ec:
            t0 += (self.tHE() + self.tEC0(ei)) / (self.tHE() + self.tEC0(ei) + self.tHC(p1))

        return self.N_Sh + self.N_T*(2*t0/bins - 1) - p1

    def eq_state(self):
        return float(sp.fsolve(self.occupation4, 10**17)) - self.n
