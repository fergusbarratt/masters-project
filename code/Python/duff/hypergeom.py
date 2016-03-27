import numpy as np
import mpmath as mp
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()


def c(kappa, delta, chi):
    return (kappa-1j*delta)/chi

def d(kappa, delta, chi):
    return (kappa+1j*delta)/chi

def z(E0, chi):
    return 2*(E0/chi)**2

def f02(a, b, z):
    return mp.hyper([], [a, b], z)

def fs_drive(chi, kappa, delta, E0s):
     return np.asarray([(chi/(kappa-1j*delta))*np.abs(E0/chi)*f02(c(kappa, delta, chi)+1, d(kappa, delta, chi), z(E0, chi))/f02(c(kappa, delta, chi), d(kappa, delta, chi), z(E0, chi)) for E0 in E0s])

kappa = 1
delta = -10
E0s = np.linspace(0, 40, 200)
chi = 10 

plt.plot(E0s, np.abs(fs_drive(chi, kappa, delta, E0s)))
plt.show()
