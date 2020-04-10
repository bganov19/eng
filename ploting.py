from PPV import *
import matplotlib.pyplot as plt

ppv = PPV()
ppv.eq_state(temp=300, light=10**16)
# intensities = [10**11]
# results = []
# temperatures = np.linspace(100, 300, 9)
#
# for n in intensities:
#     ppv.n = n
#     results.append(ppv.ls_wide_simulation(dt=10**(-2), Ti=100, Tf=300, m=9, T0=300, lstime=3600))
#
# plt.plot(temperatures, results[0], 'b')
# plt.yscale("log")
# plt.ylim(10**15, 2*10**16)
# plt.xlabel("T (K)")
# plt.ylabel("p (cm-3)")
# plt.show()
