from PPV import *
import matplotlib.pyplot as plt

ppv = PPV(temp=150, light=10**16)
p_ss = ppv.eq_state()

holes = ppv.ls_simulation(dt=0.01, lstime=1, T0=300)
time = np.linspace(0, 1, np.size(holes))

p_dark = holes[0]
dp = p_ss - p_dark


def kww(x, beta, tau):
    return [p_ss - dp*exp(-(elem*tau)**beta) for elem in x]


sigma_holes = [(p_ss - elem)/dp for elem in holes]

# popt = sp.curve_fit(kww, time[1:], holes[1:], p0=(0.5, 1000), bounds=((0, 0), (1., np.inf)))
# print(popt)
print(holes[-1])
print(p_ss)

plt.plot(time[1:], holes[1:])
# plt.plot(time[1:], sigma_holes[1:])
# plt.plot(time, (p_ss - kww(time, popt[0][0], popt[0][1]))/dp, 'r')
plt.yscale('log')
plt.xlabel("t")
plt.ylabel('$\Phi$(t)')
plt.show()