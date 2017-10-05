import matplotlib.pyplot as plt
import numpy as np
from math import pow, atanh
from numpy import tanh, sqrt

def enb(r):
    return C - K * tanh(a * r + b)

alphaMin = 0.2

C = ( 1 + alphaMin ) / 2
K = ( 1 - alphaMin ) / 2

# N = 4
# M = 8
# a = 4. / ( 10**M - 10**N )
# b = 2. - 4. / ( 1 - 10**(N-M) )
r_half_target = 10**4
delta_1 = 0.02
a = atanh ((delta_1 - 1 + C) / K) / (1 - r_half_target)
b = - r_half_target * a

r_half = - b / a
alpha_1 = enb(1)

rs = 10**(np.linspace(-1,10,200))
fs = enb(rs)

# plt.plot(rs, fs)
# ax = plt.gca()
fig, ax = plt.subplots()
ax.plot(rs, fs, label="Error norm based alpha")
ax.set_xscale("log")
ax.set_ylim(0,1)
# ax.annotate("Alpha_max", xy = (1, alpha_1), xytext=(0.1, alpha_1))
ax.axhline(alpha_1, ls='-.', color='g', label="Alpha_max")
# ax.annotate("Alpha_min", xy = (1, alphaMin), xytext=(0.1, alphaMin))
ax.axhline(alphaMin, ls='-.', color='r', label="Alpha_min")
# ax.annotate("r_h", xy = (r_half, 1), xytext=(r_half, 1.02))
ax.axvline(r_half, ls='-.', color='k', label="r_h")
plt.legend(loc="lower left")
plt.show()
