# Builds and saves an animation of the bloch sphere for a two level,
# solutions of problems defined in functions in the first part
from qutip import *
from pylab import *
from scipy import *
import numpy as np
import math as m
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

# PROBLEMS, return SX, SY, SZ expectations


def qubit_integrate(w, theta, gamma1, gamma2, psi0, tlist):
    # operators and the hamiltonian
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()
    sm = sigmam()
    H = w * (m.cos(theta) * sz + m.sin(theta) * sx)
    # collapse operators
    c_op_list = []
    n_th = 0.5  # temperature
    rate = gamma1 * (n_th + 1)
    if rate > 0.0:
        c_op_list.append(sqrt(rate) * sm)
        rate = gamma1 * n_th
        if rate > 0.0:
            c_op_list.append(sqrt(rate) * sm.dag())
            rate = gamma2
            if rate > 0.0:
                c_op_list.append(sqrt(rate) * sz)
    # evolve and calculate expectation values
    output = mesolve(H, psi0, tlist, c_op_list, [sx, sy, sz])
    return output.expect[0], output.expect[1], output.expect[2]


def solve_jc_system(E, det, tlist, psi0, g=0.5, kappa=0.7):

    # Identities
    idcavity = qeye(2)
    idqubit = qeye(2)

    # Cavity field and atomic operators
    a = tensor(destroy(2), idqubit)
    sm = tensor(idcavity, sigmam())

    # Hamiltonian components
    # Bare + int

    H0 = -det * (sm.dag() * sm + a.dag() * a) + \
    g * (sm.dag() * a + sm * a.dag())
    # Drive
    H1 = E * (a + a.dag())

    H = H0 + H1

    # Collapse operators
    c_ops = []
    c1 = m.sqrt(2 * kappa) * a
    # c2 = m.sqrt(gamma) * sm
    # more operators go here - remember to append them to c_ops
    c_ops.append(c1)
    # c_ops.append(c2)

    # Expectation operators
    e_ops = []
    sx = tensor(idcavity, sigmax())
    sy = tensor(idcavity, sigmay())
    sz = tensor(idcavity, sigmaz())
    e_ops.append(sx)
    e_ops.append(sy)
    e_ops.append(sz)
    # ground = tensor(basis(2, 0), basis(2, 1))
    # excited = tensor(basis(2, 0), basis(2, 0))
    # project_ground = ground*ground.dag()
    # project_excited = excited*excited.dag()

    out = mesolve(H, psi0, tlist, c_ops, [sx, sy, sz])
    return out.expect[0], out.expect[1], out.expect[2]


def solve_jc_system_with_detunings(E, w0, wl, wc, tlist, psi0, g=25j, kappa=2):

    # Identities
    idcavity = qeye(2)
    idqubit = qeye(2)

    # Cavity field and atomic operators
    a = tensor(destroy(2), idqubit)
    sm = tensor(idcavity, sigmam())

    H = (w0 - wl) * sm * sm.dag() + (wc - wl) * a.dag() * a + \
    g * (a.dag() * sm - sm.dag() * a) + E * (a.dag() + a)

    # Collapse operators
    c_ops = []
    c1 = m.sqrt(2 * kappa) * a
    # c2 = np.sqrt(gamma) * sm
    # more operators go here - remember to append them to c_ops
    c_ops.append(c1)
    # c_ops.append(c2)

    # Expectation operators
    e_ops = []
    sx = tensor(idcavity, sigmax())
    sy = tensor(idcavity, sigmay())
    sz = tensor(idcavity, sigmaz())
    e_ops.append(sx)
    e_ops.append(sy)
    e_ops.append(sz)
    # ground = tensor(basis(2, 0), basis(2, 1))
    # excited = tensor(basis(2, 0), basis(2, 0))
    # project_ground = ground*ground.dag()
    # project_excited = excited*excited.dag()

    out = mesolve(H, psi0, tlist, c_ops, e_ops)
    return out.expect[0], out.expect[1], out.expect[2]
# ##INITIALIZATIONS##

# Qubit Decay
# w = 1.0 * 2 * pi # qubit angular frequency
# theta = 0.2 * pi
# gamma1 = 0.5
# gamma2 = 0.2
# initial state
# a = 1.0
# psi0 = (a* basis(2,0) + (1-a)*basis(2,1))/(sqrt(a**2 + (1-a)**2))

# Jaynes-Cummings
E = 100
det = 0.0
# initial state
ground = tensor(basis(2, 0), basis(2, 1))
excited = tensor(basis(2, 0), basis(2, 0))
psi0 = ((20*excited+ground).unit()).ptrace(1)


# lists of t for solutions
tlist = linspace(0,3,300)

# Expectation values of states for ploting
sx, sy, sz = solve_jc_system(E, det, tlist, psi0, 25, 1)
# sx, sy, sz = qubit_integrate(w, theta, gamma1, gamma2, psi0, tlist)
# sx, sy, sz = solve_jc_system_with_detunings(E, 1, 1, 1, tlist, psi0, 25, 4)

# ANIMATION CODE
fig = figure()
ax = Axes3D(fig,azim=-40,elev=30)
sphere = Bloch(axes=ax)

def animate(i):
    sphere.clear()
    # sphere.add_vectors([-1,0,0])
    sphere.add_vectors([sx[i], sy[i], sz[i]])
    sphere.make_sphere()
    return ax

def init():
    sphere.vector_color = ['r']
    return ax

ani = animation.FuncAnimation(fig, animate, np.arange(len(sz)),
init_func=init, blit=True, repeat=False)
ani.save('bloch_sphere.mp4', fps=20)
