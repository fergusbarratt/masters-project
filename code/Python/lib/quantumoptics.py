""" a module for defining quantum optical systems."""

import qutip as qt
import math as m
import numpy as np
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numbers
import warnings
import time
from matplotlib import cm
import seaborn as sns
sns.set(context='poster')
plt.rcParams['image.cmap'] = 'viridis'

class QuantumOpticsSystem(object):

    '''superclass for all qoptics systems'''

    def __init__(self,
                 N_field_levels,
                 coupling=None,
                 N_qubits=1):

        # basic parameters
        self.N_field_levels = N_field_levels
        self.N_qubits = N_qubits

        if coupling is None:
            self.g = 0
        else:
            self.g = coupling

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

        arrs = to_arrays(arrays)
        
        max_len = max(map(len, arrs))

        def pad_arrs(arrs):
            '''pad a group of arrays with their last elemenst 
               to the maximum length'''
            ret = []
            for arr in arrs:
                if len(arr) >= max_len:
                    ret.append(np.asarray(arr[:max_len]))
                else:
                    ret.append(np.append(arr,
                                         arr[-1]*np.ones(max_len-len(arr))))
            return ret

        return ([arr for arr in arrs if len(arr)==max_len][0], 
                 pad_arrs(arrs))


class SteadyStateSystem(QuantumOpticsSystem):

    def __init__(self,
                 N_field_levels,
                 coupling=None,
                 N_qubits=1,
                 precalc=False
                 ):

        super().__init__(N_field_levels,
                         coupling,
                         N_qubits)

        self.precalc=precalc
        if precalc:
            self._calculate()

    def _calculate(self):
        self.rhos_ss = self.rhos()

    def rhos(self, nslice=None):
        '''rho
        return steadystate density matrices'''
        self.precalc=True

        if self.noisy:
            print('N = {}, len() = {}'.format(self.N_field_levels, 
                                              len(self.long_range)))

        def progress(*args):
            if self.noisy:
                print('|', sep='', end='', flush=True)
            return args

        if nslice is not None:
            return np.asarray([progress(qt.steadystate(ham,
                                                       self._c_ops()))
            for ham in list(self.hamiltonian())[nslice]]).T
        else:
            return np.asarray([progress(qt.steadystate(ham,
                                                       self._c_ops()))
            for ham in self.hamiltonian()]).T

    def qps(self, 
            xvec, 
            yvec, 
            start=None, 
            stop=None, 
            tr=None, 
            functype='Q'):

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
            self._calculate()
        if start is not None and stop is not None:
            if tr=='cavity':
                rhos = [rho.ptrace(0) for rho in self.rhos_ss][start:stop]
            elif tr=='qubit':
                rhos = [rho.ptrace(1) for rho in self.rhos_ss][start:stop]
            else:
                rhos = self.rhos_ss[start:stop]
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
            self._calculate()
        return np.abs(np.asarray([qt.expect(self.a*self.sm,
                rho)
                for rho in self.rhos_ss])-\
                np.asarray([qt.expect(self.a,
                rho)
                for rho in self.rhos_ss])*\
                np.asarray([qt.expect(self.sm,
                rho)
                for rho in self.rhos_ss])).T

    def g2(self):
        if not self.precalc:
            self._calculate()
        return np.abs(np.asarray([qt.expect(self.a.dag()*self.a.dag()*\
                                            self.a*self.a, 
                                            rho)/\
                                            qt.expect(self.a.dag()*self.a, 
                                                      rho)**2
                        for rho in self.rhos_ss])).T

    def abs_cavity_field(self):
        """abs_cavity_field
        Convenience function, calculates abs(expect(op(a)))"""
        if not self.precalc:
            self._calculate()
        return np.absolute([qt.expect(self.a,
                rho)
            for rho in self.rhos_ss]).T

    def purities(self):
        """purities
        Convenience function, calculates Tr(rho^2)"""
        if not self.precalc:
            self._calculate()
        return np.asarray(
                    [(rho** 2).tr() for rho in self.rhos_ss]).T

    def draw_qps(self,
                 animate=False,
                 tr='cavity',
                 colormap='inferno',
                 type='Q',
                 plottype='cf',
                 ininterval=50,
                 contno=40,
                 save=False,
                 form='mp4',
                 infigsize=(3, 3),
                 xvec=np.linspace(-8, 7, 70),
                 yvec=np.linspace(-7, 4, 70),
                 suptitle="",
                 fontdict={}):
        '''draw_qps
        Animate or plots the system quasiprobability function list using
        matplotlib
        builtins. kwargs are pretty similar to matplotlib options.
        frame rate gets set by a range length vs the ininterval
        parameter,
        infigsize gets tiled horizonally if not animating'''
        if not self.precalc:
            self._calculate()
        W = self.qps(xvec, yvec, functype=type, tr='cavity')
        if plottype == 'c' or plottype == 'cf':
            if animate:
                fig, axes = plt.subplots(1, 1, figsize=infigsize)
            else:
                fig, axes = plt.subplots(1, 
                                         len(W), 
                                         figsize=(infigsize[0]*len(W), 
                                                  infigsize[1]))
                fig.suptitle(suptitle, fontdict=fontdict)
        elif plottype == 's':
            if not animate:
                raise Exception("surface subplots not implemented")
            fig = plt.figure()
            axes = fig.gca(projection='3d')
            X, Y = np.meshgrid(xvec, yvec)
        if not animate:    
            for ax in enumerate(axes):
                if plottype=='c':
                    mble=ax[1].contour(xvec, 
                                       yvec, W[ax[0]], contno, cmap=colormap)
                elif plottype=='cf':
                    mble=ax[1].contourf(xvec, 
                                        yvec, W[ax[0]], contno, cmap=colormap)
                ax[1].set_title('{:0.2f}'.format(self.long_range[ax[0]]), 
                                loc='right')
                ax[1].set_xlabel('$\\mathbb{R}e($' + type + ')')
                ax[1].set_ylabel('$\\mathbb{I}m($' + type + ')')

                plt.colorbar(mble, ax=ax[-1])
        else:
            def init():
                if plottype == 'c':
                    plot = axes.contour(xvec, yvec, W[0], contno, cmap=colormap)
                elif plottype == 'cf':
                    plot = axes.contourf(xvec, yvec, W[0], contno, cmap=colormap)
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
                if plottype == 'c':
                    plot = axes.contour(xvec, yvec, W[i], contno, cmap=colormap)
                elif plottype == 'cf':
                    plot = axes.contourf(xvec, yvec, W[i], contno, cmap=colormap)
                elif plottype == 's':
                    Z = W[i]
                    plot = axes.plot_surface(X, Y, Z,
                                             rstride=1,
                                             cstride=1,
                                             linewidth=0,
                                             antialiased=False,
                                             shade=True,
                                             cmap=cm.coolwarm)
                    axes.set_zlim(0.0, 0.4)
                return plot

            if len(list(W)) != 1:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    anim = animation.FuncAnimation(
                        fig,
                        animate,
                        init_func=init,
                        frames=len(list(W)),
                        interval=ininterval)
            else:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    cont = axes.contour(xvec, yvec, W[0], contno)
            if save and len(list(W)) != 1:
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
        if not animate:
            fig.tight_layout()
            if save:
                plt.savefig('t={}, qp_fig.pdf'.format(time.time()))
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
                 coupling=None,
                 precalc=True):

        (self.long_range,
        self.params) = np.asarray(self._to_even_arrays([cavity_freqs,
                                  drive_freqs,
                                  drive_strengths,
                                  anharmonicity_parameters])).T

        self.length = len(self.params)

        super().__init__(N_field_levels,
                         N_qubits,
                         coupling,
                         precalc)

        if len(c_op_params)>0:
            self.kappa = c_op_params[0]
        if len(c_op_params)>1:
            self.gamma = c_op_params[1]

    def _c_ops(self):
        """_c_ops
        Build list of collapse operators
        """
        self.__def_ops()
        _c_ops = []
        if hasattr(self, 'kappa'):
            c1 = m.sqrt(2 * self.kappa) * self.a
        if hasattr(self, 'gamma'):
            c2 = m.sqrt(2 * self.gamma) * self.sm
        if 'c1' in locals():
            _c_ops.append(c1)
        if 'c2' in locals():
            _c_ops.append(c2)
        return _c_ops

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
            N_field_levels,
            noisy=False,
            precalc=True):

        self.noisy=noisy
        self.precalc=precalc

        (self.long_range,
        self.params) = self._to_even_arrays([drive_range,
                                             omega_drive_range,
                                             omega_cavity_range,
                                             omega_qubit_range])
        (self.drive_range,
         self.omega_drive_range,
         self.omega_cavity_range,
         self.omega_qubit_range) = self.params

        self.length = len(self.params)

        super().__init__(
                         N_field_levels,
                         coupling)

        if len(c_op_params)>0:
            self.kappa = c_op_params[0]
        if len(c_op_params)>1:
            self.gamma = c_op_params[1]

    def __def_ops(self):

        self.a = self.jc_a
        self.sm = self.jc_sm
        self.sx = self.jc_sx
        self.sy = self.jc_sy
        self.sz = self.jc_sz
        self.num = self.a.dag()*self.a

    def _c_ops(self):
        """_c_ops
        Build list of collapse operators
        """
        self.__def_ops()
        _c_ops = []
        if hasattr(self, 'kappa'):
            c1 = m.sqrt(2 * self.kappa) * self.a
        if hasattr(self, 'gamma'):
            c2 = m.sqrt(2 * self.gamma) * self.sm
        if 'c1' in locals():
            _c_ops.append(c1)
        if 'c2' in locals():
            _c_ops.append(c2)
        return _c_ops

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
                 initial_state=None):

        super().__init__(drive_range,
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
        self.__def_ops()

    def __def_ops(self):

        self.a = self.jc_a
        self.sm = self.jc_sm
        self.sx = self.jc_sx
        self.sy = self.jc_sy
        self.sz = self.jc_sz
        self.num = self.a.dag()*self.a

    def mesolve(self, exps=[]):
        """solve
        Interface to qutip mesolve for the system.
        :param exps: List of expectation values to calculate at
        each timestep.
        Defaults to empty.
        """
        return qt.mesolve(self.hamiltonian()[0], 
                          self.initial_state,
                          self.tlist, 
                          self._c_ops(), 
                          exps)

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
                self.hamiltonian()[0], initial_state,
                self.tlist, self._c_ops(), exps, ntraj=ntrajs)

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
                self.tlist, self._c_ops(), exps, ntraj=1)
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

    def params(self,
               drives,
               omega_cavities,
               omega_drives,
               omega_qubits,
               c_op_params,
               g,
               N):
        return (drives,
        omega_qubits,
        omega_cavities,
        omega_drives,
        c_op_params,
        g,
        N)

    def t_d_params(self,
               drives,
               omega_cavities,
               omega_drives,
               omega_qubits,
               c_op_params,
               g,
               N,
               tlist):
        return (drives,
        omega_qubits,
        omega_cavities,
        omega_drives,
        c_op_params,
        g,
        N,
        tlist)

    def det_params(self,
                  drive_strengths,
                  drive_cavity_detunings,
                  qubit_cavity_detunings,
                  c_op_params,
                  omega_cavity,
                  g,
                  N):

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
            g,
            N)
