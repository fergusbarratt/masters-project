from qutip import *
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np

def jc_hamiltonian(omega_qubit = 1.0, omega_cavity = 1.25, g = 0.05,  N = 2):
	# operator definitions
	sz = tensor(sigmaz(), identity(N))
	sm = tensor(destroy(2), identity(N))
	a = tensor(identity(2), destroy(N))
	return 0.5*omega_qubit*sz+omega_cavity*a.dag()*a+g*(a.dag()*sm+sm.dag()*a)

def kcollapse(kappa = 0.001, N=2):
	return sqrt(2*kappa)*tensor(create(2), identity(N))

def jc_liouvillian(omega_qubit=1.0, omega_cavity = 1.25, g = 0.05, kappa = 0.01, N = 2, *args):
	H = jc_hamiltonian(omega_qubit, omega_cavity, g, N)
	C  = tensor(kcollapse(kappa, 2), identity(N))
	L = liouvillian(H, [C])
	return L

psi0 = tensor(basis(2, 0), fock(2, 0))
times = np.linspace(0.0, 10.0, 100)
result = mesolve(jc_hamiltonian(), psi0, times, [kcollapse()], [tensor(sigmaz(), identity(2))])

# t = np.arange(0., 5., 0.2)
plt.plot(times, result.expect[0])
plt.show()


print("\n")