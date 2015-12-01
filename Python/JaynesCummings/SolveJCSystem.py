# IF qps is given an array of dets and drives it will go along dets
# Reparametrised system doesn't have i in interaction hamiltonian but
# unchanged sys does?
import qutip as qt
import numpy as np
import math as m
from matplotlib import pyplot as plt
from matplotlib import animation
import numbers
import warnings

class JC_System:
    def __init__(self, drive_range, omega_qubit_range, omega_cavity_range, omega_drive_range, c_op_params=[1, 1], g=10, N=40):
        # Turn all ranges into numpy arrays
        if isinstance(drive_range, numbers.Number):
            self.drive_range = np.array([drive_range])
        else:
            self.drive_range = np.array(drive_range)
        if isinstance(omega_drive_range, numbers.Number):
            self.omega_drive_range = np.array([omega_drive_range])
        else:
            self.omega_drive_range = np.array(omega_drive_range)
        if isinstance(omega_cavity_range, numbers.Number):
            self.omega_cavity_range = np.array([omega_cavity_range])
        else:
            self.omega_cavity_range = np.array(omega_cavity_range)
        if isinstance(omega_qubit_range, numbers.Number):
            self.omega_qubit_range = np.array([omega_qubit_range])
        else:
            self.omega_qubit_range = np.array(omega_qubit_range)

        if not (len(self.omega_drive_range) == len(self.omega_qubit_range) == len(self.omega_cavity_range)):
            raise ArithmeticError('Frequency ranges are not of equal length')

        # Instantiate cavity and qubit parameters
        self.g = g
        self.kappa = c_op_params[0]
        self.gamma = c_op_params[1]
        self.N = N
        self.idcavity = qt.qeye(N)
        self.idqubit = qt.qeye(2)
        self.a = qt.tensor(qt.destroy(N), self.idqubit)
        self.sm = qt.tensor(self.idcavity, qt.sigmam())

        # reparametrise - set flag if sys has not been reparametrised
        self.reparam = False
        if not (self.omega_cavity_range-self.omega_qubit_range).any():
            #any non zero cavity-qubit detunings
            self.reparam = True
            self.det_range = self.omega_cavity_range-self.omega_drive_range

class Steady_State_JC_System(JC_System):
    def __init__(self, drive_range, omega_qubit_range, omega_cavity_range, omega_drive_range, c_op_params, g, N):
        JC_System.__init__(self, drive_range, omega_qubit_range, omega_cavity_range, omega_drive_range, c_op_params, g, N)

    def hamiltonian(self, drive_index=0, det_index=0):
        #bare and int
        if self.reparam:
            hamiltonian_bare_int = -self.det_range[det_index] * (self.sm.dag() * self.sm + self.a.dag() * self.a) + self.g * (self.sm.dag() * self.a + self.sm * self.a.dag())
            # Drive
            hamiltonian_drive = self.drive_range[drive_index] * (self.a + self.a.dag())

            hamiltonian = hamiltonian_bare_int+hamiltonian_drive
            return hamiltonian
        else:
            self.q_d_det = (self.omega_qubit_range-self.omega_drive_range)[det_index]
            self.c_d_det = (self.omega_cavity_range-self.omega_drive_range)[det_index]
            self.c_q_det = -(self.omega_cavity_range-self.omega_qubit_range)[det_index]

            print '\n', self.drive_range[drive_index], '\ncavity-drive detuning: ', self.c_d_det, '\nqubit-drive detuning: ', self.q_d_det, '\ncavity-qubit detuning: ', self.c_q_det
            hamiltonian_bare = self.q_d_det*self.sm.dag()*self.sm + self.c_d_det*self.a.dag()*self.a
            hamiltonian_int = 1j*self.g*(self.a.dag()*self.sm - self.sm.dag()*self.a)
            hamiltonian_drive = self.drive_range[drive_index]*(self.a.dag()+self.a)
            hamiltonian = hamiltonian_bare + hamiltonian_int + hamiltonian_drive
            return hamiltonian

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
            for det_index, D in enumerate(self.omega_qubit_range):
                if type != 'Q':
                    qps.append(qt.wigner(self.rho_ss(drive_index, det_index).ptrace(0), xvec, yvec))
                else:
                    qps.append(qt.qfunc(self.rho_ss(drive_index, det_index).ptrace(0), xvec, yvec))
        return  qps
    def intracavity_photon_numbers(self):
        return [qt.expect(self.rho_ss(drive[0], det[0]), self.a.dag()*self.a) for drive in enumerate(self.drive_range) for det in enumerate(self.omega_qubit_range)]

    def purities(self):
        return [(self.rho_ss(drive[0], det[0])**2).tr() for drive in enumerate(self.drive_range) for det in enumerate(self.omega_qubit_range)]

class Sys_Build:
    def __init__(self, drives=[4, 5, 6], omega_qubits=1, omega_cavities=1, omega_drives=1, c_op_params=[1, 1], g=10, N=40):
        # defaults move through critical on zero det
        self.drives = drives
        self.omega_qubits = omega_qubits
        self.omega_cavities = omega_cavities
        self.omega_drives = omega_drives
        self.c_op_params = c_op_params
        self.g = g
        self.N = N

    def params(self):
        return (self.drives, self.omega_qubits, self.omega_cavities, self.omega_drives, self.c_op_params, self.g, self.N)

    def re_params(self, drives=[4, 5, 6], dets=0):
        # Completely untested
        # allows setting just drives and detunings.
        if isinstance(dets, numbers.Number):
            self.dets = [dets]
        self.omega_qubits = np.ones(len(self.dets))
        self.omega_cavities = omega_qubits
        self.omega_drives = omega_qubits+np.array(self.dets)
        return (self.drives, self.omega_qubits, self.omega_cavities, self.omega_drives, self.c_op_params, self.g, self.N)

def plot_qps(sys, type='W', infigsize=(12, 3), vec=np.linspace(-10, 10, 200)):
    W =  sys.qps(vec, vec, type)
    plots = []

    if len(W) == 1:
        fig, axes = plt.subplots(1, len(W), figsize=(6, 6))
        plots.append(axes.contourf(vec, vec, W[0], 100))
    else:
        fig, axes = plt.subplots(1, len(W), figsize=infigsize)
        for i in range(len(W)):
            plots.append(axes[i].contourf(vec, vec, W[i], 100))
    # plot the results - suppress comparison warning
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plt.show()

# drives, Q, C, D, c_op_params, g, N
# build empty system and invoke reparams(drives, dets) for old way
carmichael_system_params = Sys_Build(5, 1, 1, 2, [1, 1], 10, 40).params()

plot_qps(Steady_State_JC_System(*carmichael_system_params), 'W')

fig, axes = plt.subplots(1, len(W), figsize=(6, 6))

