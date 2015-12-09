from qutip import *
from pylab import *
from scipy import *
import math as m
import numpy as np
import matplotlib.animation as animation
from mpl_toolkits.mplot3d import Axes3D

def qubit_integrate(w, theta, gamma1, gamma2, psi0, tlist):
    # operators and the hamiltonian
    sx = sigmax(); sy = sigmay(); sz = sigmaz(); sm = sigmam()
    H = w * (m.cos(theta) * sz + m.sin(theta) * sx)
    # collapse operators
    c_op_list = []
    n_th = 0.5 # temperature
    rate = gamma1 * (n_th + 1)
    if rate > 0.0: c_op_list.append(m.sqrt(rate) * sm)
    rate = gamma1 * n_th
    if rate > 0.0: c_op_list.append(m.sqrt(rate) * sm.dag())
    rate = gamma2
    if rate > 0.0: c_op_list.append(m.sqrt(rate) * sz)
    # evolve and calculate expectation values
    output = mesolve(H, psi0, tlist, c_op_list, [sx, sy, sz])
    return output.expect[0], output.expect[1], output.expect[2]

## calculate the dynamics
w = 1.0 * 2 * pi # qubit angular frequency
theta = 0.2 * pi
gamma1 = 0.5
gamma2 = 0.2
# initial state
a = 1.0
psi_0 = (a* basis(2,0) + (1-a)*basis(2, 1))/(m.sqrt(a**2 + (1-a)**2))
tlist = linspace(0,4,250)
#expectation values for ploting
sx, sy, sz = qubit_integrate(w, theta, gamma1, gamma2, psi_0, tlist)
fig = figure()
ax = Axes3D(fig,azim=-40,elev=30)
sphere = Bloch(axes=ax)

## Animation
fig = figure()
ax = Axes3D(fig)
sphere = Bloch(axes=ax)
def animate(i):
    sphere.clear()
    sphere.add_vectors([np.sin(theta),0,np.cos(theta)])
    sphere.add_points([sx[:i+1],sy[:i+1],sz[:i+1]])
    sphere.make_sphere()
    return ax

def init():
    sphere.vector_color = ['r']
    return ax
ani = animation.FuncAnimation(fig, animate, np.arange(len(sx)),
                            init_func=init, blit=True, repeat=False)
ani.save('bloch_sphere.mp4', fps=20)
