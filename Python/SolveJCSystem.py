import qutip as qt
import numpy as np
import math as m
import matplotlib as mpl
from matplotlib import cm
import matplotlib.pyplot as plt
from timeit import default_timer as timer
import numbers

class Steady_State_JC_System:
    def __init__(self, drive_range, det_range, g=10, kappa=1, gamma=1, N=60):
        if isinstance(drive_range, numbers.Number):
            self.drive_range = [drive_range]
        else:
            self.drive_range = drive_range
        if isinstance(det_range, numbers.Number):
            self.det_range = [det_range]
        else:
            self.det_range = det_range
        self.g = g
        self.kappa = kappa
        self.gamma = gamma
        self.N = N
        self.idcavity = qt.qeye(N)
        self.idqubit = qt.qeye(2)
        self.a = qt.tensor(qt.destroy(N), self.idqubit)
        self.sm = qt.tensor(self.idcavity, qt.sigmam())
    def hamiltonian(self, drive_index=0, det_index=0):
        #bare and int
        H0 = -self.det_range[det_index] * (self.sm.dag() * self.sm + self.a.dag() * self.a) + self.g * (self.sm.dag() * self.a + self.sm * self.a.dag())
        # Drive
        H1 = self.drive_range[drive_index] * (self.a + self.a.dag())

        H = H0 + H1
        return H
    def c_ops(self, *arg):
        c_ops = []
        c1 = m.sqrt(2 * self.kappa) * self.a
        c2 = m.sqrt(2 * self.gamma) * self.sm
        c_ops.append(c1)
        c_ops.append(c2)
        c_ops+=[op for op in arg]
        return c_ops
    def rho_ss(self, drive_index, det_index):
        return qt.steadystate(self.hamiltonian(drive_index, det_index), self.c_ops())
    def qps(self, xvec, yvec, type='Q'):
        qps = []
        for drive_index, E in enumerate(self.drive_range):
            for det_index, D in enumerate(self.det_range):
                if type != 'Q':
                    qps.append(qt.wigner(self.rho_ss(drive_index, det_index).ptrace(0), xvec, yvec))
                else:
                    qps.append(qt.qfunc(self.rho_ss(drive_index, det_index).ptrace(0), xvec, yvec))
        return  qps
    def intracavity_photon_numbers(self):
        intracavity_photon_numbers = []
        for drive_index, E in enumerate(self.drive_range):
            for det_index, D in enumerate(self.det_range):
                intracavity_photon_number = qt.expect(self.rho_ss(drive_index, det_index), self.a.dag()*self.a)
                intracavity_photon_numbers.append(intracavity_photon_number)
        return intracavity_photon_numbers

    def purities(self):
        purities = []
        return [(self.rho_ss(drive[0], det[0])**2).tr() for drive in enumerate(self.drive_range) for det in enumerate(self.det_range)]

vec = np.linspace(-8, 8, 200)
carmichael = Steady_State_JC_System(5, [0, 1, 2])
for i in carmichael.purities():
    for j in carmichael.intracavity_photon_numbers():
        print i, j


# Q_1 =  carmichael.qps(vec, vec)[0]
# Q_2 = carmichael.qps(vec, vec)[1]
# Q_3 = carmichael.qps(vec, vec)[2]
# # plot the results
# fig, axes = plt.subplots(1, 3, figsize=(12,3))
# cont0 = axes[0].contourf(vec, vec, Q_1, 100)
# # bl0 = axes[0].set_title("Coherent state")
# cont1 = axes[1].contourf(vec, vec, Q_2, 100)
# # lbl1 = axes[1].set_title("Thermal state")
# cont0 = axes[2].contourf(vec, vec, Q_3, 100)
# # lbl2 = axes[2].set_title("Fock state")
# plt.show()
