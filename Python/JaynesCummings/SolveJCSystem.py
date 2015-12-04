""" IF qps is given an array of dets and drives it will go along dets
Reparametrised system doesn't have i in interaction hamiltonian but
unchanged sys does? Is this not a docstring?"""

import math as m
import qutip as qt
import numpy as np
from matplotlib import pyplot as plt
plt.rcParams['image.cmap'] = 'seismic_r'
from matplotlib import animation
import numbers
import warnings


class QuantumOpticsSystem:
    '''superclass for all qoptics systems'''

    def __init__(self, c_op_params, g, N_cavity_modes, N_qubits=1):

        # Instantiate cavity and qubit parameters
        self.g = g
        self.kappa = c_op_params[0]
        self.gamma = c_op_params[1]
        self.N_cavity_modes = N_cavity_modes
        self.N_qubits = N_qubits
        self.idcavity = qt.qeye(self.N_cavity_modes)
        self.idqubit = qt.qeye(2)
        self.a_bare = qt.destroy(self.N_cavity_modes)
        self.sm_bare = qt.sigmam()
        self.sz_bare = qt.sigmaz()

        # for back compatibility
        self.N = N_cavity_modes


class QuantumOpticsModel:
    '''superclass for all models'''

    def __init__(self):
        pass


class DickeSystem(QuantumOpticsSystem):
    '''http://rsta.royalsocietypublishing.org/content/roypta/369/19
    39/1137.full.pdf'''

    def __init__(
            self,
            omega_cavity_range,
            omega_qubit_range,
            N_qubits,
            N_cavity_modes,
            c_op_params,
            g):
        QuantumOpticsSystem.__init__(
            self, c_op_params, g, N_cavity_modes, N_qubits)
        # sort out omega ranges
        if isinstance(omega_cavity_range, numbers.Number):
            self.omega_cavity_range = np.array([omega_cavity_range])
        else:
            self.omega_cavity_range = np.array(omega_cavity_range)
        if isinstance(omega_qubit_range, numbers.Number):
            self.omega_qubit_range = np.array([omega_qubit_range])
        else:
            self.omega_qubit_range = np.array(omega_qubit_range)

        # tensor together N identities or sigmams or sigmazs
        self.idqubits = qt.tensor(
            *tuple(self.idqubit for qubit in range(N_qubits)))
        self.sz = qt.tensor(self.idcavity, qt.tensor(
            *tuple(self.sz_bare for qubit in range(N_qubits))))
        self.sm = qt.tensor(self.idcavity, qt.tensor(
            *tuple(self.sm_bare for qubit in range(N_qubits))))
        self.a = qt.tensor(self.a_bare, self.idqubits)
        if self.a.dims != self.sm.dims != self.sz.dims:
            raise ArithmeticError('dimensions issues')

    def hamiltonian(self, det_index=0):
        H_bare = self.omega_cavity_range[det_index] * (self.a.dag() * self.a)
        + self.omega_qubit_range[det_index] * self.sz
        H_int = self.g * (self.a.dag() * self.sm + self.a * self.sm.dag())
        H = H_bare + H_int
        return H

    def c_ops(self):
        c_ops = []
        c1 = m.sqrt(2 * self.kappa) * self.a
        c2 = m.sqrt(2 * self.gamma) * self.sm
        c_ops.append(c1)
        c_ops.append(c2)
        return c_ops


class JaynesCummingsSystem(QuantumOpticsSystem):
    '''Jaynes Cummings System. Takes parameter ranges, builds hamiltonians at
    each value'''

    def __init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            g,
            N):
        QuantumOpticsSystem.__init__(self, c_op_params, g, N)

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
        # set the important range data and name - for plotting against etc.
        # drive_range by default, changes in block below
        self.important_range = self.drive_range
        self.important_range_name = 'Drive'
        if not (len(self.omega_drive_range) == len(
                self.omega_qubit_range) == len(self.omega_cavity_range)):
            # test if lengths of freq vectors are the same. Correct if two
            # scalars and one vector, else raise
            if len(
                    self.omega_cavity_range) == 1 and len(
                    self.omega_drive_range) == 1:
                self.omega_cavity_range = np.ones_like(
                    omega_qubit_range) * self.omega_cavity_range[0]
                self.omega_drive_range = np.ones_like(
                    omega_qubit_range) * self.omega_drive_range[0]
                self.important_range = self.omega_qubit_range
                self.important_range_name = 'Omega Qubit'
            elif len(self.omega_qubit_range) == 1 and (
                    len(self.omega_cavity_range)) == 1:
                self.omega_qubit_range = np.ones_like(
                    omega_drive_range) * self.omega_qubit_range[0]
                self.omega_cavity_range = np.ones_like(
                    omega_drive_range) * self.omega_cavity_range[0]
                self.important_range = self.omega_drive_range
                self.important_range_name = 'Omega Drive'
            elif (len(self.omega_qubit_range)) == 1 and (
                    len(self.omega_drive_range)) == 1:
                self.omega_drive_range = np.ones_like(
                    omega_cavity_range) * self.omega_drive_range[0]
                self.omega_qubit_range = np.ones_like(
                    omega_cavity_range) * self.omega_qubit_range[0]
                self.important_range = self.omega_cavity_range
                self.important_range_name = 'Omega Cavity'
            else:
                raise ArithmeticError(
                    'Frequency ranges are not of equal length')
        # normalise parameters to kappa
        self.drive_range = self.drive_range / self.kappa
        self.omega_cavity_range = self.omega_cavity_range / self.kappa
        self.omega_qubit_range = self.omega_qubit_range / self.kappa
        self.omega_drive_range = self.omega_drive_range / self.kappa

        # reparametrise - set flag if sys has not been reparametrised
        self.reparam = False
        if not (self.omega_cavity_range - self.omega_qubit_range).any():
            # any non zero cavity-qubit detunings
            self.reparam = True
            self.det_range = self.omega_cavity_range - self.omega_drive_range

        # 1 atom 1 cavity operators
        self.a = qt.tensor(self.a_bare, self.idqubit)
        self.sm = qt.tensor(self.idcavity, self.sm_bare)

    def hamiltonian(self, drive_index=0, det_index=0):
        # bare and int
        if self.reparam:
            hamiltonian_bare_int = -self.det_range[det_index] * (
                self.sm.dag() * self.sm + self.a.dag() * self.a) + self.g * (
                self.sm.dag() * self.a + self.sm * self.a.dag())

            hamiltonian_drive = self.drive_range[drive_index] * (
                self.a + self.a.dag())

            hamiltonian = hamiltonian_bare_int + hamiltonian_drive
            return hamiltonian
        else:
            self.q_d_det = (
                self.omega_qubit_range -
                self.omega_drive_range)[det_index]
            self.c_d_det = (
                self.omega_cavity_range -
                self.omega_drive_range)[det_index]
            self.c_q_det = (
                self.omega_cavity_range -
                self.omega_qubit_range)[det_index]

            print '\ndrive strength: ',
            self.drive_range[drive_index], '\ncavity-drive detuning: ',
            self.c_d_det, '\nqubit-drive detuning: ', self.q_d_det,
            '\ncavity-qubit detuning: ', self.c_q_det

            hamiltonian_bare = self.q_d_det * self.sm.dag() * self.sm + (
                self.c_d_det * self.a.dag() * self.a)

            hamiltonian_int = 1j * self.g * (self.a.dag() * self.sm -
                                             self.sm.dag() * self.a)

            hamiltonian_drive = self.drive_range[drive_index] * (
                self.a.dag() + self.a)

            hamiltonian = hamiltonian_bare + hamiltonian_int + hamiltonian_drive
            return hamiltonian

    def c_ops(self):
        c_ops = []
        c1 = m.sqrt(2 * self.kappa) * self.a
        c2 = m.sqrt(2 * self.gamma) * self.sm
        c_ops.append(c1)
        c_ops.append(c2)
        return c_ops


class SteadyStateDickeModel(DickeSystem):
    '''does steady state things to Dicke model'''

    def __init__(
            self,
            omega_cavity_range,
            omega_qubit_range,
            N_qubits,
            N_cavity_modes,
            c_op_params,
            g):
        DickeSystem.__init__(
            self,
            omega_cavity_range,
            omega_qubit_range,
            N_qubits,
            N_cavity_modes,
            c_op_params,
            g)

    def rho(self, det_index):
        return qt.steadystate(self.hamiltonian(det_index), self.c_ops())

    def qps(self, xvec, yvec, type='Q'):
        qps = []
        for det_index in range(len(self.omega_qubit_range)):
            if type != 'Q':
                qps.append(
                    qt.wigner(
                        self.rho(det_index).ptrace(0),
                        xvec,
                        yvec))
            else:
                qps.append(qt.qfunc(self.rho(det_index).ptrace(0), xvec, yvec))
        return qps

    def intracavity_photon_numbers(self):
        return [qt.expect(self.rho(det[0]), self.a.dag() * self.a)
                for det in enumerate(self.omega_qubit_range)]

    def purities(self):
        return [(self.rho(det[0]) ** 2).tr()
                for det in enumerate(self.omega_qubit_range)]


class SteadyStateJaynesCummingsModel(JaynesCummingsSystem):
    '''does steadystate things to the super jaynes cummings system'''

    def __init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            g,
            N):
        JaynesCummingsSystem.__init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            g,
            N)

    def rho(self, drive_index, det_index):
        return qt.steadystate(
            self.hamiltonian(
                drive_index,
                det_index),
            self.c_ops())

    def qps(self, xvec, yvec, type='Q'):
        qps = []
        for drive_index in range(len(self.drive_range)):
            for det_index in range(len(self.omega_qubit_range)):
                if type != 'Q':
                    qps.append(
                        qt.wigner(
                            self.rho(
                                drive_index,
                                det_index).ptrace(0),
                            xvec,
                            yvec))
                else:
                    qps.append(
                        qt.qfunc(
                            self.rho(
                                drive_index,
                                det_index).ptrace(0),
                            xvec,
                            yvec))
        return qps

    def intracavity_photon_numbers(self):
        return [
            qt.expect(
                self.rho(
                    drive[0],
                    det[0]),
                self.a.dag() *
                self.a) for drive in enumerate(
                self.drive_range) for det in enumerate(
                    self.omega_qubit_range)]

    def purities(self):
        return [(self.rho(drive[0], det[0]) ** 2).tr()
                for drive in enumerate(self.drive_range)
                for det in enumerate(self.omega_qubit_range)]


class JaynesCummingsParameters:
    ''' interface to ssjcm class for unpacking parameters and
    reparametrising'''

    def __init__(
            self,
            drives=[
                4,
                5,
                6],
            omega_qubits=1,
            omega_cavities=1,
            omega_drives=1,
            c_op_params=[
                1,
                1],
            g=10,
            N=40):
        # defaults move through critical on zero det
        self.drives = drives
        self.omega_qubits = omega_qubits
        self.omega_cavities = omega_cavities
        self.omega_drives = omega_drives
        self.c_op_params = c_op_params
        self.g = g
        self.N = N

    def params(self):
        return (
            self.drives,
            self.omega_qubits,
            self.omega_cavities,
            self.omega_drives,
            self.c_op_params,
            self.g,
            self.N)

    def re_params(self, drives=[4, 5, 6], dets=0):
        # Completely untested
        # allows setting just drives and detunings.
        if isinstance(dets, numbers.Number):
            self.dets = [dets]
        self.omega_qubits = np.ones(len(self.dets))
        self.omega_cavities = self.omega_qubits
        self.omega_drives = self.omega_qubits + np.array(self.dets)
        return (
            self.drives,
            self.omega_qubits,
            self.omega_cavities,
            self.omega_drives,
            self.c_op_params,
            self.g,
            self.N)


def plot_qps(sys, type='Q', infigsize=(12, 3), vec=np.linspace(-10, 10, 200)):
    W = sys.qps(vec, vec, type)
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


def draw_qps(sys, type='Q', ininterval=50, contno=100, save=False, infigsize=(
        6, 6), xvec=np.linspace(-8, 7, 50), yvec=np.linspace(-10, 4, 50),
        form='mp4'):

    W = zip(sys.important_range, sys.qps(xvec, yvec, type))  # SYSTEM SPECIFIC

    fig, axes = plt.subplots(1, 1, figsize=infigsize)

    def init():
        cont = axes.contourf(xvec, yvec, W[0][1], contno)
        return cont

    def animate(i):
        axes.cla()
        plt.cla()
        plt.title(sys.important_range_name +
                  ': %d' % W[i][0])  # SYSTEM SPECIFIC
        cont = axes.contourf(xvec, yvec, W[i][1], contno)
        return cont

    if len(W) != 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            anim = animation.FuncAnimation(
                fig,
                animate,
                init_func=init,
                frames=len(W),
                interval=ininterval)
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            cont = axes.contourf(xvec, yvec, W[0], contno)
    if save and len(W) != 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            if form == 'mp4':
                anim.save(
                    'qp_anim.mp4', fps=30, extra_args=[
                        '-vcodec', 'libx264'])
            if form == 'gif':
                anim.save('qp_anim.gif', writer='imagemagick', fps=4)
    plt.show()


def plot_exp(sys, expval='iphnum'):
    xs = sys.important_range
    if expval == 'iphnum':
        ys = sys.intracavity_photon_numbers()
    elif expval == 'purities':
        ys = sys.purities()
    else:
        raise NameError('no such expectation value')
    plt.plot(xs, ys)
    plt.show()

# omega_cavity_range, omega_qubit_range, N_qubits, N_cavity_modes,
# c_op_params, g
# dicke = SteadyStateDickeModel(2, 2, 3, 10, [1, 1], 10)
# plot_qps(dicke)
# print dicke.hamiltonian()

# drives, Q, C, D, c_op_params, g, N
# build empty system and invoke reparams(drives, dets) for old way
carmichael_system_params = JaynesCummingsParameters(
    5, 100, 100, np.linspace(95, 105, 50), [1, 1], 10, 40).params()

carmichael_sys = SteadyStateJaynesCummingsModel(*carmichael_system_params)

draw_qps(carmichael_sys, 'W', 50, 100, False)
# plot_qps(carmichael_sys, 'Q')
# plot_exp(carmichael_sys)


def draw_figs(estart, estep, estop, detstart, detstep, detstop):
    '''broken'''
    fig = plt.figure()
    ax = plt.axes()
    line, = ax.plot([], [])

    def init():
        line.set_data([], [])
        return line,

    def animate(frame):
        system = SteadyStateJaynesCummingsModel(*JaynesCummingsParameters(
            estart+frame/(1/estep), 100, 100,
            np.linspace(detstart, detstop, 1 / detstep),
            [1, 1], 10, 10).params())
        # line.set_data(sys.important_range, sys.intracavity_photon_numbers())
        xs, ys = system.important_range, system.intracavity_photon_numbers()
        return plt.plot(xs, ys)

    # def animate(i):
    #     x = np.linspace(0, 2, 1000)
    #     y = np.sin(2 * np.pi * (x - 0.01 * i))
    #     line.set_data(x, y)
    #     return line,

    anim = animation.FuncAnimation(
        fig,
        animate,
        init_func=init,
        frames=int(estop),
        interval=20,
        blit=False,
        repeat=False)

    plt.show()

# INTERESTING
# carmichael_system_params = JaynesCummingsParameters(5, np.linspace(60, 140,
# 20), 100, 100, [1, 0.8], 10, 40).params()
