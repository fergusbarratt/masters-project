'''Produce 6 pane resonant plot'''
'''Setup''''''Setup''''''Setup'''

import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='talk', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 10, 10 

import numpy as np
import qutip as qt

import time as T
import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

### Parameters
matrix_size = 20 

# Dispersive
c_q_det = (8.1831-10.5665)
kappa_disp = 0.0024
gamma = 0
coupling_strength = 0.3347
cavity_freq = 10.5665
sigmaz=-1

# Normalising Parameters
xiC1= (abs(c_q_det)*kappa_disp)**(3/2)/(3**(3/4)*coupling_strength**2)
chi0=coupling_strength**2/abs(c_q_det)

def disp_drive(A, omega_cavity, omega_drive, sigmaz, det, g, kappa):
  '''self-consistently determine xi from A'''
  def chi(A):
    return sigmaz*(g**2)/np.sqrt(2*g**2*(A**2+sigmaz)+det**2)
  return np.sqrt(A**2 * (1/omega_cavity**2)*(
            (omega_drive**2-(omega_cavity-chi(A))**2)**2 + kappa_disp**2*omega_drive**2))

# Resonant
kappa=1
g=10

def res_drive(alphasq, detuning, kappa, g):
    return abs(-1j*np.sqrt(alphasq)*(
            kappa - 1j*(detuning + (detuning/abs(detuning))*g**2/np.sqrt(detuning**2 + 4 * g**2 * alphasq))))


# Q Function Method
def abs_cavity_field_Q(sys, xvec=np.linspace(-5, 5, 50), yvec=np.linspace(-5, 5, 50)):
    Qs = sys.qps(xvec, yvec)
    return np.array([abs(np.sum(xvec.T.dot(Q))+1j*np.sum(Q.dot(yvec))) for Q in Qs])/34


### Resonant Plot

fig = plt.figure()

gs = gridspec.GridSpec(4, 3)

ax = [None for _ in range(4*3)] 
ax[0] = fig.add_subplot(gs[0:2, :-1])
ax[1] = fig.add_subplot(gs[2:4, :-1], sharex = ax[0])
for i in enumerate(range(2, 6)):
    ax[i[1]] = fig.add_subplot(gs[i[0], -1])

q_res_drives = np.linspace(4, 5.5, 8)

### Contour plot in top right
## Semiclassical
width = 200
res_alpha = np.linspace(0, 0.15, 100)
res_detrange = np.linspace(-20, 20, width)
res_drives = np.array([[res_drive(a, det, kappa, g) for det in res_detrange] for a in res_alpha])

mble = ax[0].contour(res_detrange, res_alpha, res_drives, levels=q_res_drives, cmap='jet')
ax[0].set_title('Semiclassical', loc='right', fontdict={'fontsize': 12, 'verticalalignment': 'bottom'})
ax[0].set_xlabel('$\omega_c-\omega_d$')
ax[0].set_ylabel('$\left|A\\right|^2$')
res_cbar = plt.colorbar(mble, ax=ax[1], 
                 label='Drive Strength',
                 orientation='horizontal')

# Wider colorbar lines
res_childs = res_cbar.ax.get_children()
res_childs[0].set_linewidths(20)

### Line Plots in bottom left
## Quantum 
q_res_detrange = res_detrange

carmichael_parameters = [qo.JaynesCummingsParameters(10, matrix_size).det_params(
              drive_strengths=q_res_drive,
              drive_cavity_detunings=q_res_detrange,
              qubit_cavity_detunings=0,
              c_op_params=[1, 1],
              omega_cavity=10) for q_res_drive in q_res_drives]

carmichael_systems = [qo.SteadyStateJaynesCummingsModel(*c_params) for c_params in carmichael_parameters]

for sys in enumerate(carmichael_systems):
    ax[1].plot(q_res_detrange, sys[1].abs_cavity_field()**2, c= res_childs[0].get_colors()[sys[0]])
ax[1].set_title('Quantum', loc='right', fontdict={'fontsize': 12, 'verticalalignment': 'bottom'})
ax[1].set_xlabel('$\omega_c-\omega_d$')
ax[1].set_ylabel('$|\langle a \\rangle | ^2 $')

### Q function plots along right hand side
xvec = np.linspace(-8, 8, 100)
yvec = np.linspace(-8, 8, 100)
def plot_drive_QPs(syss, ind, xvec, yvec):
    return [sys.qps(xvec, yvec)[ind] for sys in syss]

for qp in enumerate(plot_drive_QPs(carmichael_systems[::2][0:4], width//2, xvec, yvec)):
    '''for every other system, truncated to four. get colors for qp backgrounds from 
    lines'''
    ax[qp[0]+2].contour(xvec, yvec, qp[1], 40, linewidths=0.8, cmap='inferno')
    ax[qp[0]+2].set_xlabel('Im(Q)')
    ax[qp[0]+2].set_ylabel('Re(Q)')
    ax[qp[0]+2].set_title('Q', loc='right', fontdict={'fontsize':12, 'verticalalignment':'bottom', 'color': res_childs[0].get_colors()[::2][0:4][qp[0]], 'weight':'bold'})
    
gs.update(wspace=0.5, hspace=0.6)
plt.suptitle('Spontaneous Dressed State Polarisation')
plt.show()
