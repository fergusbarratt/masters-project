'1d linear fokker planck solutions'''

import matplotlib
matplotlib.use('TkAgg')
import numpy as np
from numpy import Inf
from scipy.linalg import expm
from scipy.integrate import quad as integr
import matplotlib.pyplot as plt


class diffuse:

    def __init__(self, A, D, initial_conds=[1, 1, 0]):
        self.A = np.array([[A]])
        self.D = np.array([[D]])
        self.funcs = [self.x, self.sig, self.broadening]
        self.initial_conds = [1, 1, 0]
        self.normval = None

    def x(self, t):
        return self.initial_conds[0]*expm(self.A*t)[0]

    def sig(self, t):
        return (self.initial_conds[1]*expm(2*self.A*t) - (self.D/(2*self.A))*(1-expm(2*self.A*t)))[0]

    def broadening(self, x, t):
        def unnormalised(x, t):
            return ((1/np.sqrt(2*np.pi*(self.D/2*self.A)*(expm(2*self.A*t)-1)))*expm(-0.5*((x-self.initial_conds[2]*expm(self.A*t))**2)/(self.D/(2*self.A))*(expm(2*self.A*t - 1))))[0]
        if not self.normval: self.normval = integr(unnormalised, -Inf, Inf, args=(t))
        return unnormalised(x, t)/self.normval[0]

funcs = diffuse(0.1, 0.1).funcs
x = funcs[0]
sig = funcs[1]
broadening = funcs[2]
X = np.linspace(-3, 3, 100)
# X = [x(t) for t in T]
# Y = [sig(t) for t in T]
Z = [broadening(x, 5) for x in X]
# q, = plt.plot(T, X)
# p, = plt.plot(T, Y)
o, = plt.plot(X, Z)
# plt.legend([q, p], ["$x$", "$sig^2$"])
plt.xlabel("Time")
plt.ylabel("Value")
plt.show()
#
