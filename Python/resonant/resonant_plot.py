'''Produce 6 pane resonant plot'''

print('starting setup ...', end='')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='talk', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 12, 8 

import numpy as np
import qutip as qt

import time as T
import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

### Parameters
matrix_size = 85 

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

q_res_drives = np.linspace(4, 5.5, 8)

print('done')
print('Plot setup ... ', end='')

### Resonant Plot
fig = plt.figure()
gs = gridspec.GridSpec(4, 3)
ax = [None for _ in range(4*3)] 
ax[0] = fig.add_subplot(gs[0:2, :-1])
ax[1] = fig.add_subplot(gs[2:4, :-1], sharex = ax[0])
for i in enumerate(range(2, 6)):
    ax[i[1]] = fig.add_subplot(gs[i[0], -1])
print('done')

### Contour plot in top right
## Semiclassical
print('starting semiclassical ... ', end='')
width = 200
res_alpha = np.linspace(0, 0.15, 100)
res_detrange = np.linspace(-20, 20, width)
res_drives = np.array(
        [[res_drive(a, det, kappa, g) for det in res_detrange] 
            for a in res_alpha])

# Plot Contours at q_res_drive on res_drives
mble = ax[0].contour(res_detrange, 
                     res_alpha, 
                     res_drives, 
                     levels=q_res_drives, 
                     cmap='jet')

ax[0].set_title('Semiclassical', 
                loc='right', 
                fontdict={'fontsize': 12, 'verticalalignment': 'bottom'})
ax[0].set_xlabel('$\omega_c-\omega_d$')
ax[0].set_ylabel('$\left|A\\right|^2$')

# Colorbar
res_cbar = plt.colorbar(mble, ax=ax[1], 
                 label='Drive Strength',
                 orientation='horizontal',
                 pad=0.2)

print('done')

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
          omega_cavity=10) 
          for q_res_drive in q_res_drives]

print('starting field setup ...', end='', flush=True)

carmichael_systems = [
        qo.SteadyStateJaynesCummingsModel(
                                          *c_params, 
                                          noisy=True) 
        for c_params in carmichael_parameters]

print('done')

print('starting field data ... ', end='', flush=True)

# preallocate data array
data = np.empty((len(carmichael_systems), len(q_res_detrange)))

for sys in enumerate(carmichael_systems):

    # calculate cavity field (squared)
    cav_field = sys[1].abs_cavity_field()**2

    # store cavity field (squared)
    data[sys[0]] = cav_field

    # plot cavity field (squared) -
    # get colors from semiclassical contours

    ax[1].plot(q_res_detrange, 
               cav_field, 
               c= res_childs[0].get_colors()[sys[0]])

print('done')

# Set plot parameters    
ax[1].set_title('Quantum', 
                loc='right', 
                fontdict={'fontsize': 12, 'verticalalignment': 'bottom'})
ax[1].set_xlabel('$\omega_c-\omega_d$')
ax[1].set_ylabel('$|\langle a \\rangle | ^2 $')

print('starting qps ', end='', flush=True)

### Q function plots along right hand side

xvec = np.linspace(-8, 8, 100)
yvec = np.linspace(-8, 8, 100)

def plot_drive_QPs(syss, ind, xvec, yvec):
    return [sys.qps(xvec, yvec)[ind] for sys in syss]

for qp in enumerate(plot_drive_QPs(carmichael_systems[::2][0:4], 
                                   width//2, 
                                   xvec, 
                                   yvec)):
    '''for every other system, truncated to four. 
    get colors for qp backgrounds from lines'''
    print('.', end='')
    ax[qp[0]+2].contour(xvec, yvec, qp[1], 40, 
            linewidths=0.8, 
            cmap='inferno')

    # Set plot parameters
    ax[qp[0]+2].set_xlabel('Im(Q)')
    ax[qp[0]+2].set_ylabel('Re(Q)')
    ax[qp[0]+2].set_title('Q', 
                      loc='right', 
                      fontdict={'fontsize':12, 'verticalalignment':'bottom', 'color': res_childs[0].get_colors()[::2][0:4][qp[0]], 
                      'weight':'bold'})

print('done') 

print('saving data and displaying plot ... ', end='', flush=True)

# update grid spacings
gs.update(wspace=0.35, hspace=0.53)

# total plot title
fig.suptitle('Spontaneous Dressed State Polarisation', 
             x=0.2, 
             y=0.93,
             verticalalignment='bottom', 
             fontsize=16)

# save figure and data
fig.savefig('resonant.pdf')
np.save('resonant.npy', data)

print('done')

plt.show()
