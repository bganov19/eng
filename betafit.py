from PPV import *
import matplotlib.pyplot as plt

ppv = PPV(temp=100, light=10**11)
p_ss = ppv.p0(ppv.eq_state(temp=100, light=10**11))
holes = ppv.ls_simulation(dt=10**0, lstime=3600, T0=300)
time = np.linspace(0, 3600, 3600)

p_dark = holes[0]
dp = p_ss - p_dark

ppv.E_EC = 0.1
# print(p_ss)


def F(t, beta): return [p_ss - dp*exp(-(elem*ppv.tEC())**beta) for elem in t]


popt, pcov = sp.curve_fit(F, time, holes[:-1])
print(popt)
teo_curve = F(time, *popt)

plt.plot(time, holes[:-1])
plt.plot(time, teo_curve, 'r')
plt.yscale('log')
plt.show()

