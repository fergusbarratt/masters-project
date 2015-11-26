import qutip as qt
from qutip import wigner, qfunc, steadystate
import numpy as np
import math as m
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt

# Cavity Resolution
N = 40

def jc_hamiltonian(E, det, g = 10):
    # Identities
    idcavity = qt.qeye(N)
    idqubit = qt.qeye(2)

    # Cavity field and atomic operators
    a = qt.tensor(qt.destroy(N), idqubit)
    sm = qt.tensor(idcavity, qt.sigmam())

    #bare and int
    H0 = -det * (sm.dag() * sm + a.dag() * a) + g * (sm.dag() * a + sm * a.dag())
    # Drive
    H1 = E * (a + a.dag())

    H = H0 + H1
    return H

def c_ops(kappa=1, gamma=0):
    # Identities
    idcavity = qt.qeye(N)
    idqubit = qt.qeye(2)

    # Cavity field and atomic operators
    a = qt.tensor(qt.destroy(N), idqubit)
    sm = qt.tensor(idcavity, qt.sigmam())

    c_ops = []
    c1 = m.sqrt(2 * kappa) * a
    c2 = m.sqrt(2*gamma) * sm
    c_ops.append(c1)
    c_ops.append(c2)
    return c_ops

# Q Functions
# xvec = np.linspace(-8, 8, 200)
# yvec = np.linspace(-8, 8, 200)
# qf = qt.wigner(rho_ss, xvec, yvec)

xvec = np.linspace(-8,8,200)
def jc_qp(E, D, fl = 'q'):
    if fl != 'q':
        return wigner(steadystate(jc_hamiltonian(E, D), c_ops()).ptrace(0), xvec, xvec)
    else:
        return qfunc(steadystate(jc_hamiltonian(E, D), c_ops()).ptrace(0), xvec, xvec)
W_1 = jc_qp(4, 0)
W_2 = jc_qp(5, 0)
W_3 = jc_qp(6, 0)
# plot the results
fig, axes = plt.subplots(1, 3, figsize=(12,3))
cont0 = axes[0].contourf(xvec, xvec, W_1, 100)
# bl0 = axes[0].set_title("Coherent state")
cont1 = axes[1].contourf(xvec, xvec, W_2, 100)
# lbl1 = axes[1].set_title("Thermal state")
cont0 = axes[2].contourf(xvec, xvec, W_3, 100)
# lbl2 = axes[2].set_title("Fock state")
plt.show()

