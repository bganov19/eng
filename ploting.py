from PPV import *
import matplotlib.pyplot as plt

ppv = PPV(temp=300, light=10**14)

ss = []

holes = ppv.short_ls_wide_simulation(T0=300, Ti=100, Tf=300, m=9, lstime=3600)

temps = range(100, 325, 25)
#
for t in temps:
    ppv.T = t
    ss.append(ppv.eq_state())

plt.scatter(temps, ss, label="steady state")
plt.scatter(temps, holes, label='after 3600s LS')
plt.yscale("log")
plt.xlabel("T (K)")
plt.ylabel("p (cm-3)")
plt.legend()
plt.show()
