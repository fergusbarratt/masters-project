from qutip import *
from pylab import *
from scipy import *
import numpy as np
import matplotlib.animation as animation 
from mpl_toolkits.mplot3d import Axes3D

def qubit_integrate(w, theta, gamma1, gamma2, psi0, tlist):
    # operators and the hamiltonian
	sx = sigmax(); sy = sigmay(); sz = sigmaz(); sm = sigmam() 
	H = w * (cos(theta) * sz + sin(theta) * sx)
	# collapse operators
	c_op_list = []
	n_th = 0.5 # temperature
	rate = gamma1 * (n_th + 1)
	if rate > 0.0: c_op_list.append(sqrt(rate) * sm)
	rate = gamma1 * n_th
	if rate > 0.0: c_op_list.append(sqrt(rate) * sm.dag()) 
	rate = gamma2
	if rate > 0.0: c_op_list.append(sqrt(rate) * sz)
    # evolve and calculate expectation values
	output = mesolve(H, psi0, tlist, c_op_list, [sx, sy, sz])
	return output.expect[0], output.expect[1], output.expect[2]

def solve_jc_system(E, det, tlist, g = 0.5, kappa = 0.1):
	# Initial State
	psi0=tensor(basis(2, 0), basis(2, 0))

	# Identities
	idcavity = qeye(2)
	idqubit = qeye(2)

	# Cavity field and atomic operators
	a = tensor(destroy(2), idqubit)
	sm = tensor(idcavity, sigmam())

	# Hamiltonian components
	# Bare + int

	H0 = (sm.dag() * sm + a.dag() * a ) + g * ( sm.dag() * a + sm*a.dag() )
	# Drive
	H1 = E * ( a + a.dag() )

	H = H0+H1

	# Collapse operators
	c_ops = []
	c1 = np.sqrt(2*kappa)*a 
	# more operators go here
	c_ops.append(c1)

	# Expectation operators
	sx = sigmax()
	sy = sigmay()
	sz = sigmaz()

	out = mesolve(H, psi0, tlist, c_ops, [sx, sy, sz])
	return out.expect[0], out.expect[1], out.expect[3]
tlist = np.linspace(10, 20)
