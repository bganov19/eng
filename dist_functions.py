from math import exp, sqrt, pi, log, log10, floor
import scipy.integrate as integrate
import numpy as np

C = integrate.quad(lambda e: 1/sqrt(e), 0.1, 0.5)[0]**(-1)


def grain(E):
    if 0.1 <= E <= 0.25:
        return C/sqrt(E)
    else:
        return 0


def grain_cum(E):
    if E < 0.1:
        return 0
    elif E <= 0.25:
        return 2*C*(sqrt(E)-sqrt(0.1))
    else:
        return 1


def inv_cum(E):
    return np.square(E/(2*C) + sqrt(0.1))


def grain_dist(bins):
    d = grain(0.1)
    u = np.random.uniform(0.1, 0.2, bins)
    y = np.random.uniform(0, d, bins)
    ok = [elem1 < grain(elem2) for elem1, elem2 in zip(y, u)]
    u = u[ok]
    return u
