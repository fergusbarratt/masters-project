""" a module for defining quantum optical systems.
Wrap up interesting parameters using the wrapper class and feed them
to the model, and use available methods on that model object"""

import qutip as qt
import math as m
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

    def __init__(self,
                 c_op_params,
                 N_field_modes,
                 N_qubits=1,
                 coupling=0):

        # Instantiate cavity and qubit parameters
        self.g = coupling
        self.kappa = c_op_params[0]
        if len(c_op_params)>1:
            self.gamma = c_op_params[1]
        self.N_field_modes = N_field_modes
        self.N_qubits = N_qubits
        self.idcavity = qt.qeye(self.N_field_modes)
        self.idqubit = qt.qeye(2)
        self.a_bare = qt.destroy(self.N_field_modes)
        self.sm_bare = qt.sigmam()
        self.sz_bare = qt.sigmaz()
        self.sx_bare = qt.sigmax()
        self.sy_bare = qt.sigmay()

        # for back compatibility TOFIX (N cavity modes)
        self.N = N_field_modes

        # Figure, axes
        self.fig = plt.figure()
        self.ax = plt.axes3D(fig, azim=-40, elev=30)

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
            N_field_modes,
            c_op_params,
            g):
        QuantumOpticsSystem.__init__(
            self, c_op_params, g, N_field_modes, N_qubits)
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
        if hasattr(self, 'gamma'):
            c2 = m.sqrt(2 * self.gamma) * self.sm
        c_ops.append(c1)
        if c2 in locals():
            c_ops.append(c2)
        return c_ops

class QuantumDuffingOscillator(QuantumOpticsSystem):

    ''' Walls, Drummond, Quantum Theory of Optical Bistability I Model '''

    def __init__(self, cavity_freqs, anharmonicity_parameter):
        self.omega_cavity_range = cavity_freq
        self.anharmonicity_parameter=anharmonicity_parameter

    def hamiltonian(self, drive_index=0, det_index=0):
        return self.omega_c * self.a.dag() * self.a + self.anharmonicity_parameter * self.a.dag() ** 2 * self.a ** 2


class JaynesCummingsSystem(QuantumOpticsSystem):

    '''Jaynes Cummings System. Takes parameter ranges, builds hamiltonians
    at each value. JaynesCummingsParameters for wrapping'''

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

        QuantumOpticsSystem.__init__(self, c_op_params, N, coupling=g)

        self.noisy = noisy

        def to_arrays(arrs):
            ret = []
            for arr in arrs:
                if isinstance(arr, numbers.Number):
                    ret.append(np.asarray([arr]))
                else:
                    ret.append(np.asarray(arr))
            return ret

        (self.drive_range,
         self.omega_drive_range,
         self.omega_cavity_range,
         self.omega_qubit_range) = to_arrays([drive_range,
                                             omega_drive_range,
                                             omega_cavity_range,
                                             omega_qubit_range])

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
        """hamiltonian
        Build the steadystate hamiltonian.

        :param drive_index: chooses system drive for hamiltonian. Defaults to zero
        :param det_index: chooses system detuning for hamiltonian. Defaults to zero
        """
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

    def rho(self, drive_index=0, det_index=0):
        """rho
        Solve for the steady state system density matrix using the drive-detuning
        parameters at drive_index, det_index
        :param drive_index: index of Drive in internal list of drives. Defaults to
        zero
        :param det_index: index of Detuning in internal detunings lists. Defaults
        to zero
        """
        return qt.steadystate(
            self.hamiltonian(
                drive_index,
                det_index),
            self.c_ops())

    def c_ops(self):
        """c_ops
        Build list of collapse operators
        """
        c_ops = []
        c1 = m.sqrt(2 * self.kappa) * self.a
        if hasattr(self, 'gamma'):
            c2 = m.sqrt(2 * self.gamma) * self.sm
        c_ops.append(c1)
        if 'c2' in locals():
            c_ops.append(c2)
        return c_ops

    def qps(self, xvec, yvec, type='Q'):
        """qps
        returns an array of Q or W functions for all system parameters.
        :param xvec: X vector over which function is evaluated F(X+iY)
        :param yvec: Y vector over which function is evaluated F(X+iY)
        :param type: 'Q' or 'W' : Husimi Q or Wigner
        """
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

    def draw_qps(self,
                 type='Q',
                 plottype='c',
                 ininterval=50,
                 contno=100,
                 save=False,
                 form='mp4',
                 infigsize=(6, 6),
                 xvec=np.linspace(-8, 7, 70),
                 yvec=np.linspace(-10, 4, 70)):
        '''draw_qps
        Animate the system quasiprobability function list using matplotlib
        builtins. kwargs are pretty similar to matplotlib options.
        frame rate gets set by a range length vs the ininterval
        parameter'''
        W = list(zip(self.important_range, self.qps(xvec, yvec, type)))
        if plottype == 'c' or plottype == 'cf':
            fig, axes = plt.subplots(1, 1, figsize=infigsize)
        elif plottype == 's':
            fig = plt.figure()
            axes = fig.gca(projection='3d')
            X, Y = np.meshgrid(xvec, yvec)

        def init():
            if plottype == 'c':
                plot = axes.contour(xvec, yvec, W[0][1], contno)
            elif plottype == 'cf':
                plot = axes.contourf(xvec, yvec, W[0][1], contno)
            elif plottype == 's':
                Z = W[0][1]
                plot = axes.plot_surface(
                    X, Y, Z, rstride=1, cstride=1, linewidth=0,
                    antialiased=True, shade=True, cmap=cm.coolwarm)
                axes.set_zlim(0.0, 0.1)
            return plot

        def animate(i):
            axes.cla()
            plt.cla()
            plt.title(self.important_range_name +
                      ': %d' % W[i][0])
            if plottype == 'c':
                plot = axes.contour(xvec, yvec, W[i][1], contno)
            elif plottype == 'cf':
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
                cont = axes.contour(xvec, yvec, W[0][1], contno)
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


class SteadyStateJaynesCummingsModel(JaynesCummingsSystem):
    """SteadyStateJaynesCummingsModel
    Steady state modelling of the Jaynes Cummings system. Contains
    convenience functions for easy calculation of: the absolute
    cavity field, the purity"""

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
        JaynesCummingsSystem.__init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            g,
            N)
        '''precalculates steady state density matrices. set
	noisy flag for progress bar'''
        self.rhos_ss = []
        for drive in enumerate(self.drive_range):
            if noisy:
                print(
 '\ndrive={0:.4f}, len(drives)={1}, len(detunings)={2}, N={3}'.format(
     drive[1],
     len(self.drive_range),
     len(self.omega_drive_range),
     self.N))
            for det in enumerate(self.omega_drive_range):
                if noisy:
                    print('|', end='', flush=True)
                self.rhos_ss.append(self.rho(drive[0], det[0]))

    class qp_list(object):
        '''calculate on demand list'''
        def __init__(self, xvec, yvec, density_matrix_list, functype='Q'):
            self.xvec = xvec
            self.yvec = yvec
            self.functype = functype
            self.density_matrix_list = density_matrix_list

        def _qps(self, new_dens_mat_list=None):
            if new_dens_mat_list == None:
                new_dens_mat_list = self.density_matrix_list
            self.qps = []
            for rho in new_dens_mat_list:
                if self.functype != 'Q':
                    self.qps.append(
                        qt.wigner(
                            rho,
                            self.xvec,
                            self.yvec))
                else:
                    qps.append(
                        qt.qfunc(
                            rho,
                            self.xvec,
                            self.yvec))
            return self.qps

        def __getitem__(self, slice):
            '''calculates the q functions only along the 
            requested slice'''
            return self._qps(self.xvec, self.yvec, 
                    self.density_matrix_list[slice], 
                    self.functype)

    def qps(self, xvec, yvec, functype='Q', slice = None):
        """qps
        overrides jcm qps to use precalculated density matrices if available
        pass slice kwarg to generate only a subset
        :param xvec: list of xvalues to calculate QP function at
        :param yvec: list of yvalues to calculate QP function at
        """
        self.rhos_cav_ss = [rho.ptrace(0) for rho in self.rhos_ss]
        if slice is None:
            qps = self.qp_list(xvec, yvec, self.rhos_cav_ss, functype)[:]
        else:
            qps = self.qp_list(xvec, yvec, self.rhos_cav_ss, functype)[slice]
        return qps 

    def draw_bloch_sphere(self):
        """draw the qubit bloch sphere for the system steady states"""
        self.rhos_qb_ss = [rho.ptrace(1) for rho in self.rhos_ss]
        self.b_sphere = qt.Bloch()
        self.b_sphere.add_states(self.rhos_qb_ss)
        self.b_sphere.show()

    def correlator(self):
        """correlator
        Measure of quantum vs semiclassical"""
        return np.abs(np.asarray([qt.expect(self.a*self.sm,
                rho)
                for rho in self.rhos_ss])-\
                np.asarray([qt.expect(self.a,
                rho)
                for rho in self.rhos_ss])*\
                np.asarray([qt.expect(self.sm,
                rho)
                for rho in self.rhos_ss]))

    def abs_cavity_field(self):
        """abs_cavity_field
        Convenience function, calculates abs(expect(op(a)))"""
        return np.absolute([qt.expect(self.a,
                rho)
            for rho in self.rhos_ss])

    def purities(self):
        """purities
        Convenience function, calculates Tr(rho^2)"""
        return np.asarray(
                [(rho** 2).tr() for rho in self.rhos_ss])

class TimeDependentJaynesCummingsModel(JaynesCummingsSystem):
    """TimeDependentJaynesCummingsModel
    Time dependent modelling of the Jaynes-Cummings System. Takes
    two additional parameters, a list of times: tlist and an
    initial state: initial_state, default None"""

    def __init__(self,
                 drive_range,
                 omega_qubit_range,
                 omega_cavity_range,
                 omega_drive_range,
                 c_op_params,
                 g,
                 N,
                 tlist,
                 initial_state=None,
                 noisy=False):

        JaynesCummingsSystem.__init__(self,
                                                drive_range,
                                                omega_qubit_range,
                                                omega_cavity_range,
                                                omega_drive_range,
                                                c_op_params,
                                                g,
                                                N,
                                                noisy)
        self.tlist = tlist
        if initial_state is None:
            self.initial_state = qt.tensor(self.idcavity, self.idqubit)
        else:
            self.initial_state = initial_state

    def mesolve(self, exps=[]):
        """solve
        Interface to qutip mesolve for the system.
        :param exps: List of expectation values to calculate at
        each timestep.
        Defaults to empty.
        """
        return qt.mesolve(self.hamiltonian(), self.initial_state,
                          self.tlist, self.c_ops(), exps)

    def mcsolve(self, ntrajs=500, exps=[], initial_state=None):
        """mcsolve
        Interface to qutip mcsolve for the system
        :param ntrajs: number of quantum trajectories to average.
        Default is QuTiP
        default of 500
        :param exps: List of expectation values to calculate at
        each timestep
        """
        if initial_state is None:
            initial_state = qt.tensor(
                    qt.basis(self.N, 0), qt.basis(2, 0))
        return qt.mcsolve(
                self.hamiltonian(), initial_state,
                self.tlist, self.c_ops(), exps, ntraj=ntrajs)

    def trajectory(self, exps=None, initial_state=None, draw=False):
        '''for convenience. Calculates the trajectory of an
        observable for one montecarlo run. Default expectation is
        cavity amplitude, default initial state is bipartite
        vacuum. todo: draw: draw trajectory on bloch sphere.
        Write in terms of mcsolve??'''
        if exps is None or draw is True:
            exps = []
        if initial_state is None:
            initial_state = qt.tensor(
                    qt.basis(self.N, 0), qt.basis(2, 0))

        self.one_traj_soln = qt.mcsolve(
                self.hamiltonian(), initial_state,
                self.tlist, self.c_ops(), exps, ntraj=1)
        if self.noisy:
            print(self.one_traj_soln.states[0][2].ptrace(1))

        if not draw:
            return self.one_traj_soln
        else:
            self.b_sphere = qt.Bloch()
            self.b_sphere.add_states(
             [state.ptrace(1) for state in self.one_traj_soln.states[0]],
             'point')
            self.b_sphere.point_markers=['o']
            self.b_sphere.size = (10, 10)
            self.b_sphere.show()

class JaynesCummingsParameters:

    ''' interface to ssjcm class for unpacking parameters and
    reparametrising'''

    def __init__(
            self,
            g,
            N):
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
                  omega_cavity):
        # Completely untested
        # allows setting just drives and detunings.
        self.drives = drive_strengths
        self.omega_qubits = np.asarray(
     [omega_cavity + qcd for qcd in np.atleast_1d(qubit_cavity_detunings)])
        self.omega_drives = np.asarray(
     [omega_cavity + dcd for dcd in np.atleast_1d(drive_cavity_detunings)])
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

if __name__ == '__main__':
    main()
