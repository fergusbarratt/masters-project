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

class QuantumOpticsSystem(object):

    '''superclass for all qoptics systems'''

    def __init__(self,
                 N_field_levels,
                 c_op_params,
                 N_qubits=1,
                 coupling=None):

        # basic parameters
        self.N = N_field_levels
        self.N_qubits = N_qubits

        if coupling is None:
            self.g = 0
        else:
            self.g = coupling

        self.kappa = c_op_params[0]
        if len(c_op_params)>1:
            self.gamma = c_op_params[1]

        # bare operators
        self.idcavity = qt.qeye(self.N_field_modes)
        self.idqubit = qt.qeye(2)
        self.a_bare = qt.destroy(self.N_field_modes)
        self.sm_bare = qt.sigmam()
        self.sz_bare = qt.sigmaz()
        self.sx_bare = qt.sigmax()
        self.sy_bare = qt.sigmay()

        # 1 atom 1 cavity operators
        self.a = qt.tensor(self.a_bare, self.idqubit)
        self.sm = qt.tensor(self.idcavity, self.sm_bare)
        self.sx = qt.tensor(self.idcavity, self.sx_bare)
        self.sy = qt.tensor(self.idcavity, self.sy_bare)
        self.sz = qt.tensor(self.idcavity, self.sz_bare)

        # Figure, axes
        self.fig = plt.figure()
        self.ax = plt.axes3D(fig, azim=-40, elev=30)

    def _to_even_arrays(self, arrays):
        ''' Takes a list of arrays and pads them all with zeros to
        the length of the longest'''
        def to_arrays(arrs):
            ret = []
            for arr in arrs:
                if isinstance(arr, numbers.Number):
                    ret.append(np.asarray([arr]))
                else:
                    ret.append(np.asarray(arr))
            return ret

        def pad_arr(arr, new_len):
            if len(arr) > new_len:
                return np.asarray(arr[:new_len])
            else:
                return np.append(arr, np.zeros(new_len-len(arr)))

        arrs = to_arrays(arrays)

        ret_arrs = []
        max_len = max(map(len, arrs))

        for arr in enumerate(arrs):
            if len(arr) < max_len:
                arrs[arr[0]] = pad_arr(arr[1], max_len)
        return arrs

class SteadyStateSystem(QuantumOpticsSystem):

    def __init__(self,
                 param_indices,
                 N_field_levels,
                 c_op_params,
                 N_qubits=1,
                 coupling=None):
        QuantumOpticsSystem.__init__(self,
                                     N_field_levels,
                                     c_op_params,
                                     N_qubits=1,
                                     coupling=None):

        self.rhos_ss = []
        for index in enumerate(param_indices):
            self.rhos_ss.append(self.rho(*index))

    def rho(self, index):
        '''rho
        return steadystate density matrices at param_indices'''
        return qt.steadystate(
            self.hamiltonian(
                *index),
            self._c_ops())

    def _c_ops(self):
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


    def qps(self, xvec, yvec, rhos=None, functype='Q', slice = None):
        """qps
        overrides ssm qps to use precalculated density matrices if
        available
        pass slice or rhos kwarg to generate only a subset
        :param xvec: list of xvalues to calculate QP function at
        :param yvec: list of yvalues to calculate QP function at
        """

        class qp_list(object):
            """qps
            returns an array of Q or W functions for all 
            system parameters.
            :param xvec: X vector over which function is evaluated 
            F(X+iY)
            :param yvec: Y vector over which function is evaluated 
            F(X+iY)
            :param type: 'Q' or 'W' : Husimi Q or Wigner
            """
            def __init__(self,
                         xvec,
                         yvec,
                         density_matrix_list,
                         functype='Q'):

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

        if rhos=None:
            rhos = self.rhos_ss

        if slice is None:
            qps = self.qp_list(xvec, yvec, rhos, functype)[:]
        else:
            qps = self.qp_list(xvec, yvec, rhos, functype)[slice]
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
        Animate the system quasiprobability function list using
        matplotlib
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

class QuantumDuffingOscillator(SteadyStateSystem):
    ''' Walls, Drummond, Quantum Theory Optical Bistability I Model '''

    def __init__(self,
                 cavity_freqs,
                 anharmonicity_parameter,
                 N_field_levels,
                 c_op_params,
                 N_qubits=1,
                 coupling=None):

        SteadyStateSystem.__init__(N_field_levels,
                 c_op_params,
                 N_qubits=1,
                 coupling=None)

        self.params = self._to_even_arrays([cavity_freqs,
                                            anharmonicity_parameter])
        self.param_indices = ((a, b) 
                for a in range(len(self.params[0]))
                for b in range(len(self.params[1])))

    def hamiltonian(self):

        hamiltonians = np.asarray([self.params[0][omega_index] * self.a.dag() * self.a + \
        self.params[1][anh_index] * self.a.dag() ** 2 * self.a ** 2 
        for omega_index, anh_index in self.param_indices])

        return hamiltonians

class JaynesCummingsSystem(QuantumOpticsSystem):

    def __init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            g,
            N):

        QuantumOpticsSystem.__init__(self, c_op_params, N, coupling=g)

         self.params = self._to_even_arrays([drive_range,
                                             omega_drive_range,
                                             omega_cavity_range,
                                             omega_qubit_range])
         (self.drive_range,
          self.omega_drive_range,
          self.omega_cavity_range,
          self.omega_qubit_range) = self.params

        self.param_indices = ((a, b) 
                for a in range(len(self.params[0])) 
                for b in range(len(self.params[1])))

    def hamiltonian(self):
        """hamiltonian
        Build the steadystate hamiltonians.
        """
        self.q_d_det = (
            self.omega_qubit_range -
            self.omega_drive_range)
        self.c_d_det = (
            self.omega_cavity_range -
            self.omega_drive_range)
        self.c_q_det = (
            self.omega_cavity_range -
            self.omega_qubit_range)

        hamiltonian_bare = np.asarray(
                [q_d_det * self.sm.dag() * self.sm + (
            c_d_det * self.a.dag() * self.a) 
            for q_d_det in self.q_d_det
            for c_d_det in self.c_d_det]

        hamiltonian_int = np.ones_like(hamiltonian_bare)*\
             1j * self.g * (self.a.dag() * self.sm -
                            self.sm.dag() * self.a)

        hamiltonian_drive = np.asarray([
            drive * (
            self.a.dag() + self.a)
            for drive in self.drive_range]) 

        hamiltonians = hamiltonian_bare + \
            hamiltonian_int + hamiltonian_drive

        return hamiltonians


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
