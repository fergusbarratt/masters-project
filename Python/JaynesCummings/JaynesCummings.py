# from qutip import *
# from math import sqrt, cos

import matplotlib.pyplot as plt
import numpy as np
from qutip import *

# Example from manual
g = 0.5 # coupling strength
kappa = 0.01 # Cavity decay rate
t = np.linspace(0, 100, 200) # Define time vector
N = 2 # Set where to truncate Fock state for cavity
M = 2 # Number of qubit levels.
omega_qubit = 1.0
omega_cavity = 1.0
omega_driving = 1.0
driving_strength = 10
# kappa = 0.01
# g = 0.05
# N = 2

#States
# ustate = basis(M, 0)
excited = basis(M, 1)
ground = basis(M, 0)

# Operators
sz = tensor(sigmaz(), qeye(N))
sm = tensor(destroy(2), qeye(N))
a = tensor(destroy(N), qeye(M))
ada = tensor(num(N), qeye(M))
c_ops = [] # Build collapse operators
c_ops.append(np.sqrt(kappa) * a)

# Initial state and expectation projectors
psi0 = tensor(basis(N, 1), ground) # Define initial state

state_GG = tensor(basis(N, 1), ground) # Define states onto which to project
sigma_GG = state_GG * state_GG.dag()
# state_UU = tensor(basis(N, 0), ustate)
# sigma_UU = state_UU * state_UU.dag()

# Hamiltonian Components
H0 = 0.5 * omega_qubit*sz + omega_cavity*a.dag()*a + g*(a.dag()*sm + sm.dag()*a)
H1 = ((driving_strength)/np.sqrt(2))*(a+a.dag())
def H1_coeff(t, args, omega_driving=1.0):
	return np.cos(omega_driving*t)
H = [H0, [H1, H1_coeff]]

# Solve, Plot
output = mesolve(H, psi0, t, c_ops, [ada, sigma_GG])

plt.plot(t, output.expect[0])
plt.show()


print("\n")
