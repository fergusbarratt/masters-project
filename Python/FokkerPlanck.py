'''1d linear fokker planck solutions'''

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from scipy.linalg import expm
import matplotlib.pyplot as plt


def diffuse(A, D, initial_conds=[1, 1, 0]):
    A = np.array([[A]])
    D = np.array([[D]])
    def x(t):
        return initial_conds[0]*expm(A*t)[0]

    def sig(t):
        return (initial_conds[1]*expm(2*A*t) - (D/(2*A))*(1-expm(2*A*t)))[0]

    def broadening(x, t):
        return ((1/np.sqrt(2*np.pi*(D/2*A)*(expm(2*A*t)-1)))*expm(-0.5*((x-initial_conds[2]*expm(A*t))**2)/(D/(2*A))*(expm(2*A*t - 1))))[0]
    return [x, sig, broadening]


funcs = diffuse(0.1, 0.1)
x = funcs[0]
sig = funcs[1]
broadening = funcs[2]
T = np.linspace(-10, 10, 100)
X = [x(t) for t in T]
Y = [sig(t) for t in T]
# Z = [broadening(x, 5) for x in X]
q, = plt.plot(T, X)
p, = plt.plot(T, Y)
# o, = plt.plot(X, Z)
plt.legend([q, p], ["$x$", "$sig^2$"])
plt.xlabel("Time")
plt.ylabel("Value")
plt.show()
