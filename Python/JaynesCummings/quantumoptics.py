""" a module for defining quantum optical systems. Wrap up interesting parameters using the wrapper class and feed them to the interesting model, and use available methods on that model object"""

from __future__ import print_function
import math as m
import qutip as qt
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numbers
import warnings
from matplotlib import cm
plt.rcParams['image.cmap'] = 'viridis'


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
        self.sx_bare = qt.sigmax()
        self.sy_bare = qt.sigmay()

        # for back compatibility TOFIX (N cavity modes)
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
        self.omega_qubit_range = np.atleast_1d(np.asarray(omega_qubit_range))
        self.omega_cavity_range = np.atleast_1d(np.asarray(omega_cavity_range))

        # replicate non important range to length of important one
        if len(
               np.atleast_1d(self.omega_cavity_range)) == len(
               np.atleast_1d(self.omega_qubit_range)):
            self.important_range = self.omega_cavity_range
        elif len(np.atleast_1d(omega_cavity_range)) == 1:
            self.important_range = self.omega_qubit_range
            self.omega_cavity_range = self.omega_cavity_range[
                0]*np.ones_like(self.omega_qubit_range)
        elif len(np.atleast_1d(self.omega_qubit_range)) == 1:
            self.important_range = self.omega_cavity_range
            self.omega_qubit_range = self.omega_qubit_range * \
                np.ones_like(self.omega_cavity_range)
        else:
            raise ArithmeticError('wrong input form')

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
        print(self.omega_cavity_range[0], self.omega_qubit_range[0])
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
    each value. JaynesCummingsParameters for wrapping'''

    def __init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            g,
            N,
            noisy=False):
        QuantumOpticsSystem.__init__(self, c_op_params, g, N)
        self.noisy = noisy
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
            # scalars and one vector, else raise - bad to mix important
            # range name business and parameters TOFIX
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
        self.sx = qt.tensor(self.idcavity, self.sx_bare)
        self.sy = qt.tensor(self.idcavity, self.sy_bare)
        self.sz = qt.tensor(self.idcavity, self.sz_bare)

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
            if self.noisy:
                print(
                    '\ndrive strength: ',
                    self.drive_range[drive_index],
                    '\ncavity-drive detuning: ',
                    self.c_d_det,
                    '\nqubit-drive detuning: ',
                    self.q_d_det,
                    '\ncavity-qubit detuning: ',
                    self.c_q_det)

            hamiltonian_bare = self.q_d_det * self.sm.dag() * self.sm + (
                self.c_d_det * self.a.dag() * self.a)

            hamiltonian_int = 1j * self.g * (self.a.dag() * self.sm -
                                             self.sm.dag() * self.a)

            hamiltonian_drive = self.drive_range[drive_index] * (
                self.a.dag() + self.a)

            hamiltonian = hamiltonian_bare + \
                hamiltonian_int + hamiltonian_drive
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
        for det_index in range(len(self.important_range)):
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
                for det in enumerate(self.important_range)]

    def purities(self):
        return [(self.rho(det[0]) ** 2).tr()
                for det in enumerate(self.important_range)]


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

    def plot_qps(self, type='Q', infigsize=(12, 3),
                 vec=np.linspace(-10, 10, 200)):
        W = self.qps(vec, vec, type)
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

    def draw_qps(self, type='Q', plottype='c', ininterval=50, contno=100,
                 save=False, form='mp4', infigsize=(6, 6),
                 xvec=np.linspace(-8, 7, 70), yvec=np.linspace(-10, 4, 70)):

        W = list(zip(self.important_range, self.qps(xvec, yvec, type)))
        if plottype == 'c':
            fig, axes = plt.subplots(1, 1, figsize=infigsize)
        elif plottype == 's':
            fig = plt.figure()
            axes = fig.gca(projection='3d')
            X, Y = np.meshgrid(xvec, yvec)

        def init():
            if plottype == 'c':
                plot = axes.contourf(xvec, yvec, W[0][1], contno)
            elif plottype == 's':
                Z = W[0][1]
                plot = axes.plot_surface(
                    X, Y, Z, rstride=1, cstride=1, linewidth=0,
                    antialiased=True, shade=True, cmap=cm.coolwarm)
                axes.set_zlim(0.0, 0.1)
            if plottype == 'c':
                plt.colorbar(plot)
            return plot

        def animate(i):
            axes.cla()
            plt.cla()
            plt.title(self.important_range_name +
                      ': %d' % W[i][0])
            if plottype == 'c':
                plot = axes.contourf(xvec, yvec, W[i][1], contno)
            elif plottype == 's':
                Z = W[i][1]
                plot = axes.plot_surface(
                    X, Y, Z, rstride=1, cstride=1, linewidth=0,
                    antialiased=False, shade=True, cmap=cm.coolwarm)
                axes.set_zlim(0.0, 0.4)
            return plot

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
                cont = axes.contourf(xvec, yvec, W[0][1], contno)
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

    def plot_exp(self, expval='iphnum'):
        xs = self.important_range
        if expval == 'iphnum':
            ys = self.intracavity_photon_numbers()
        elif expval == 'purities':
            ys = self.purities()
        else:
            raise NameError('no such expectation value')
        plt.plot(xs, ys)
        plt.show()


class TimeDependentJaynesCummingsModel(SteadyStateJaynesCummingsModel):

    def __init__(self,
                 drive_range,
                 omega_qubit_range,
                 omega_cavity_range,
                 omega_drive_range,
                 c_op_params,
                 g,
                 N,
                 tlist,
                 initial_state=None):

        SteadyStateJaynesCummingsModel.__init__(self,
                                                drive_range,
                                                omega_qubit_range,
                                                omega_cavity_range,
                                                omega_drive_range,
                                                c_op_params,
                                                g,
                                                N)
        self.tlist = tlist
        if initial_state is None:
            self.initial_state = qt.tensor(self.idcavity, self.idqubit)
        else:
            self.initial_state = initial_state

    def solve(self, exps=[]):
        return qt.mesolve(self.hamiltonian(), self.initial_state,
                          self.tlist, self.c_ops(), exps)

    def show_time_dependence(self, exps):
        plt.plot(self.tlist, self.solve(exps).expect[0])
        plt.show()

    def draw_bloch_sphere(self, ininterval=1, contno=100,
                          save=False, form='mp4', infigsize=(6, 6)):
        soln = self.solve([self.sx, self.sy, self.sz])
        sx, sy, sz = soln.expect[0], soln.expect[1], soln.expect[2]
        vecs = [np.real(qt.Qobj([[sx[i], sy[i], sz[i]]]).unit().full())[0]
                for i in range(1, len(sx))]
        fig = plt.figure()
        ax = Axes3D(fig)
        sphere = qt.Bloch(axes=ax, fig=fig)
        sphere.add_vectors(vecs[3])
        sphere.show()


class JaynesCummingsParameters:

    ''' interface to ssjcm class for unpacking parameters and
    reparametrising'''

    def __init__(
            self,
            g=10,
            N=40):
        self.g = g
        self.N = N

    def params(self,
               drives,
               omega_cavities,
               omega_drives,
               omega_qubits,
               c_op_params,
               ):
        return (drives,
        omega_qubits,
        omega_cavities,
        omega_drives,
        c_op_params,
        self.g,
        self.N)

    def t_d_params(self,
               drives,
               omega_cavities,
               omega_drives,
               omega_qubits,
               c_op_params,
               tlist):
        return (drives,
        omega_qubits,
        omega_cavities,
        omega_drives,
        c_op_params,
        self.g,
        self.N,
        tlist)

    def det_params(self,
                  drive_strengths,
                  drive_cavity_detunings,
                  qubit_cavity_detunings,
                  c_op_params,
                  omega_cavity=1000000):
        # Completely untested
        # allows setting just drives and detunings.
        self.drives = drive_strengths
        self.omega_qubits = np.asarray(
                                       [omega_cavity + \
                                        qcd for qcd in np.atleast_1d(qubit_cavity_detunings)])
        self.omega_drives = np.asarray([omega_cavity + \
                                        dcd for dcd in np.atleast_1d(drive_cavity_detunings)])
        self.omega_cavities = np.asarray([omega_cavity])
        self.c_op_params = c_op_params
        return (
            self.drives,
            self.omega_qubits,
            self.omega_cavities,
            self.omega_drives,
            self.c_op_params,
            self.g,
            self.N)

class TestSystem(object):
    def __init__(self, system):
        self.system = system

    def test_jaynes_cummings_system(sys_type='normal', draw_bloch_sphere=False, time_dependent=False, draw_qps=False, plot_exp=False, ham=False):

        # time dependent and non time dependent systems
        def test_sys(sys_type):
            '''produces a system, either 'normal' (parameters given explicitly)
            time_dependent, or det_param (parameters given as default middle and
            two detunings)'''
            if sys_type=='normal':
                carmichael_system_params = \
                    JaynesCummingsParameters(10, 40).params(5, 100, 100, 100, [1, 1])
                return SteadyStateJaynesCummingsModel(*carmichael_system_params)
            elif sys_type=='time_dependent':
                time_dependent_carmichael_sys_params = \
                    JaynesCummingsParameters(10, 40).t_d_params(5, 100, 100, 100, [1, 1], np.linspace(0, 4, 100))
                return TimeDependentJaynesCummingsModel(*time_dependent_carmichael_sys_params)
            elif sys_type=='det_param':
                det_sys_params = \
                    JaynesCummingsParameters(10, 40).det_params(5, np.linspace(-3, 3, 20), 1000, [1, 1])
                return SteadyStateJaynesCummingsModel(*det_sys_params)

        # time_dependent_carmichael_sys.show_time_dependence(num)
        if draw_bloch_sphere and sys_type=='time_dependent':
            test_sys(sys_type).draw_bloch_sphere()
        if draw_qps:
            test_sys(sys_type).draw_qps()
        if plot_exp:
            test_sys(sys_type).plot_exp()
        if ham:
            print(carmichael_sys.hamiltonian())
