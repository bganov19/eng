from PPV import *
import matplotlib.pyplot as plt

ppv = PPV()
ppv.eq_state(temp=300, light=10**12)

# light = [10**14, 10**15]
sigmas = [10**(-19), 10**(-18), 10**(-17), 10**(-16), 10**(-15), 10**(-14), 10**(-13)]
light5p = [1.4*10**14, 1.4*10**13, 1.4*10**12, 2.1*10**11, 2.2*10**12, 2.2*10**13, 2.2*10**14]
results_n = []
results_10n = []

# i = 6
# ppv.sigma_h = sigmas[i]
# ppv.n = light5p[i]
# print(ppv.ls_simulation(dt=10**(-8), lstime=10**(-3), T0=300)[-1])
# ppv.n = light5p[i]*10
# print(ppv.ls_simulation(dt=10**(-8), lstime=10**(-3), T0=300)[-1])

p5 = [5.08, 5.08, 5.07, 5.04, 5.1, 5.1, 5.1]
p10 = [10.95, 10.95, 10.85, 8.2, 8.26, 8.26, 8.25]
res_p = [elem2/elem1 for elem1, elem2 in zip(p5, p10)]
rel_sig = [elem/ppv.sigma_e for elem in sigmas]
print(rel_sig)
print(res_p)

plt.scatter(rel_sig, res_p)
plt.plot(rel_sig, res_p)
plt.xscale('log')
plt.xlabel('sigma_h/sigma_e')
plt.ylabel('p(10n)/p(n)')
plt.show()

# for l, s in zip(light5p, sigmas):
#     ppv.sigma_h = s
#     ppv.n = l
#     results_n.append(ppv.ls_simulation(dt=10**(-8), lstime=10**(-2), T0=300)[-1])
    # ppv.n = 10*l
    # results_10n.append(ppv.ls_simulation(dt=10**(-8), lstime=10**(-2), T0=300)[-1])

# results = [p10n/pn for pn, p10n in zip(results_n, results_10n)]
# sigmas_rel = [elem/ppv.sigma_e for elem in sigmas]

# print(results_n)
# print(results)
# plt.plot(sigmas_rel, results)
# plt.xscale('log')

# plt.show()

# x = np.linspace(-2., 2.5, 1000)
# plt.plot(x, (x/1.25)**2)
# plt.plot(x, (x/1.25-1.35)**2 + 0.7)
# plt.plot(x, (x/1.25)**2 + 1.1)
# plt.ylim(-0.1, 2)
# plt.xlim(-2, 2.5)
# plt.xticks([])
# plt.yticks([])
# plt.xlabel("Q")
# plt.ylabel("E")
# plt.show()
