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
                 coupling=None,
                 N_qubits=1):

        # basic parameters
        self.N_field_levels = N_field_levels
        self.N_qubits = N_qubits

        if coupling is None:
            self.g = 0
        else:
            self.g = coupling

        if len(c_op_params)>0:
            self.kappa = c_op_params[0]
        if len(c_op_params)>1:
            self.gamma = c_op_params[1]

        # bare operators
        self.idcavity = qt.qeye(self.N_field_levels)
        self.idqubit = qt.qeye(2)
        self.a_bare = qt.destroy(self.N_field_levels)
        self.sm_bare = qt.sigmam()
        self.sz_bare = qt.sigmaz()
        self.sx_bare = qt.sigmax()
        self.sy_bare = qt.sigmay()

        # 1 atom 1 cavity operators
        self.jc_a = qt.tensor(self.a_bare, self.idqubit)
        self.jc_sm = qt.tensor(self.idcavity, self.sm_bare)
        self.jc_sx = qt.tensor(self.idcavity, self.sx_bare)
        self.jc_sy = qt.tensor(self.idcavity, self.sy_bare)
        self.jc_sz = qt.tensor(self.idcavity, self.sz_bare)

    def _to_even_arrays(self, arrays):
        ''' Takes a list of arrays and pads them all with
        last element to the length of the longest'''
        def to_arrays(arrs):
            ''' convert list of numbers and arrays to 1
            and many element arrays'''
            ret = []
            for arr in arrs:
                if isinstance(arr, numbers.Number):
                    ret.append(np.asarray([arr]))
                else:
                    ret.append(np.asarray(arr))
            return ret

        def pad_arr(arr, new_len):
            '''pad arrays with final element of array'''
            if len(arr) > new_len:
                return np.asarray(arr[:new_len])
            else:
                return np.append(arr,
                                 arr[-1]*np.ones(new_len-len(arr)))

        arrs = to_arrays(arrays)

        ret_arrs = []
        max_len = max(map(len, arrs))

        for arr in enumerate(arrs):
            arrs[arr[0]] = pad_arr(arr[1], max_len)
        return arrs

class SteadyStateSystem(QuantumOpticsSystem):

    def __init__(self,
                 N_field_levels,
                 c_op_params,
                 coupling=None,
                 N_qubits=1,
                 precalc=True):

        super().__init__(N_field_levels,
                         c_op_params,
                         coupling,
                         N_qubits)

        if precalc:
            self._calculate()
            self.precalc = precalc

    def _calculate(self):
        self.rhos_ss = self.rhos()

    def _c_ops(self):
        """c_ops
        Build list of collapse operators
        """
        c_ops = []
        if hasattr(self, 'kappa'):
            c1 = m.sqrt(2 * self.kappa) * self.a
        if hasattr(self, 'gamma'):
            c2 = m.sqrt(2 * self.gamma) * self.sm
        if 'c1' in locals():
            c_ops.append(c1)
        if 'c2' in locals():
            c_ops.append(c2)
        return c_ops

    def rhos(self, nslice=None):
        '''rho
        return steadystate density matrices'''
        precalc=True
        if nslice is not None:
            return np.asarray([qt.steadystate(ham,
            self._c_ops())
            for ham in list(self.hamiltonian())[nslice]])
        return np.asarray([qt.steadystate(ham,
            self._c_ops()) for ham in self.hamiltonian()])

    def qps(self, xvec, yvec, nslice=None, tr=None, functype='Q'):

        class qp_list(object):
            """qps
            lazy calculate qps
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

            def qps(self):
                self.qps = []
                for rho in self.density_matrix_list:
                    if self.functype != 'Q':
                        self.qps.append(
                            qt.wigner(
                                rho,
                                self.xvec,
                                self.yvec))
                    else:
                        self.qps.append(
                            qt.qfunc(
                                rho,
                                self.xvec,
                                self.yvec))
                return self.qps

        if not self.precalc:
            self._calculate_rhos_ss(nslice)
        if nslice is not None:
            if tr=='cavity':
                rhos = [rho.ptrace(0) for rho in self.rhos_ss][nslice]
            elif tr=='qubit':
                rhos = [rho.ptrace(1) for rho in self.rhos_ss][nslice]
            else:
                rhos = self.rhos_ss[nslice]
        else:
            if tr=='cavity':
                rhos = [rho.ptrace(0) for rho in self.rhos_ss]
            elif tr=='qubit':
                rhos = [rho.ptrace(1) for rho in self.rhos_ss]
            else:
                rhos = self.rhos_ss

        qps = qp_list(xvec, yvec,
                      rhos,
                      functype).qps()

        return qps

    def correlator(self):
        """correlator
        Measure of quantum vs semiclassical"""
        if not self.precalc:
            self._calculate_rhos_ss()
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
        if not self.precalc:
            self._calculate_rhos_ss()
        return np.absolute([qt.expect(self.a,
                rho)
            for rho in self.rhos_ss])

    def purities(self):
        """purities
        Convenience function, calculates Tr(rho^2)"""
        if not self.precalc:
            self._calculate_rhos_ss()
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
        if not self.precalc:
            self._calculate_rhos_ss()
        W = enumerate(self.qps(xvec, yvec, type))
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
                plot = axes.plot_surface(X, Y, Z,
                                         rstride=1,
                                         cstride=1,
                                         linewidth=0,
                                         antialiased=True,
                                         shade=True,
                                         cmap=cm.coolwarm)
                axes.set_zlim(0.0, 0.1)
            return plot

        def animate(i):
            axes.cla()
            plt.cla()
            plt.title(': %d' % W[i][0])
            if plottype == 'c':
                plot = axes.contour(xvec, yvec, W[i][1], contno)
            elif plottype == 'cf':
                plot = axes.contourf(xvec, yvec, W[i][1], contno)
            elif plottype == 's':
                Z = W[i][1]
                plot = axes.plot_surface(X, Y, Z,
                                         rstride=1,
                                         cstride=1,
                                         linewidth=0,
                                         antialiased=False,
                                         shade=True,
                                         cmap=cm.coolwarm)
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
                    anim.save('qp_anim.mp4',
                              fps=30,
                              extra_args=['-vcodec', 'libx264'])
                if form == 'gif':
                    anim.save('qp_anim.gif',
                              writer='imagemagick',
                              fps=4)
        plt.show()

class QuantumDuffingOscillator(SteadyStateSystem):
    ''' Walls, Drummond, Quantum Theory Optical Bistability I Model '''

    def __init__(self,
                 drive_strengths,
                 cavity_freqs,
                 drive_freqs,
                 anharmonicity_parameters,
                 N_field_levels,
                 c_op_params,
                 N_qubits=1,
                 coupling=None):

        self.params = np.asarray(self._to_even_arrays([cavity_freqs,
                                        drive_freqs,
                                        drive_strengths,
                                        anharmonicity_parameters])).T

        self.length = len(self.params)

        super().__init__(N_field_levels,
                         c_op_params,
                         N_qubits,
                         coupling)

    def __def_ops(self):

        self.a = self.a_bare
        self.sm = self.sm_bare
        self.sx = self.sx_bare
        self.sy = self.sy_bare
        self.sz = self.sz_bare

    def hamiltonian(self):

        self.__def_ops()

        hamiltonians_bare = np.asarray(
                        [(omega_c-omega_d) * self.a.dag() * self.a + \
                         anh * self.a.dag() ** 2 * self.a ** 2
                            for omega_c,
                                omega_d,
                                _,
                                anh in self.params])

        hamiltonians_drive = np.asarray(
                        [dr_str * (self.a.dag() + self.a)
                            for _,
                                _,
                                dr_str,
                                _,  in self.params])

        hamiltonians = hamiltonians_bare + hamiltonians_drive

        return hamiltonians

class JaynesCummingsSystem(SteadyStateSystem):

    def __init__(
            self,
            drive_range,
            omega_qubit_range,
            omega_cavity_range,
            omega_drive_range,
            c_op_params,
            coupling,
            N_field_levels):

        self.params = self._to_even_arrays([drive_range,
                                            omega_drive_range,
                                            omega_cavity_range,
                                            omega_qubit_range])
        (self.drive_range,
         self.omega_drive_range,
         self.omega_cavity_range,
         self.omega_qubit_range) = self.params

        self.length = len(self.params)

        super().__init__(N_field_levels,
                         c_op_params,
                         coupling)

    def __def_ops(self):

        self.a = self.jc_a
        self.sm = self.jc_sm
        self.sx = self.jc_sx
        self.sy = self.jc_sy
        self.sz = self.jc_sz

    def hamiltonian(self):

        self.__def_ops()

        self.q_d_det = (
            self.omega_qubit_range -
            self.omega_drive_range)
        self.c_d_det = (
            self.omega_cavity_range -
            self.omega_drive_range)
        self.c_q_det = (
            self.omega_cavity_range -
            self.omega_qubit_range)

        self.hamiltonian_bare = np.asarray(
                [q_d_det * self.sm.dag() * self.sm + (
                 c_d_det * self.a.dag() * self.a)
            for q_d_det, c_d_det in zip(self.q_d_det, self.c_d_det)])

        self.hamiltonian_int = np.ones_like(self.hamiltonian_bare)*\
             1j * self.g * (self.a.dag() * self.sm -
                            self.sm.dag() * self.a)

        self.hamiltonian_drive = np.asarray([
            drive * (
            self.a.dag() + self.a)
            for drive in self.drive_range])

        hamiltonians = self.hamiltonian_bare + \
            self.hamiltonian_int + self.hamiltonian_drive

        return hamiltonians

    def draw_bloch_sphere(self):
        """draw the qubit bloch sphere for the system steady states"""
        self.rhos_qb_ss = [rho.ptrace(1) for rho in self.rhos_ss]
        self.b_sphere = qt.Bloch()
        self.b_sphere.add_states(self.rhos_qb_ss)
        self.b_sphere.show()

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

        super().__init__(drive_range,
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
                    qt.basis(self.N_field_levels, 0), qt.basis(2, 0))
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
                    qt.basis(self.N_field_levels, 0), qt.basis(2, 0))

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
