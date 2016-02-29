'''Produce dispersive plot'''


import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='talk', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 15, 15

import numpy as np
import qutip as qt

import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

### Parameters
matrix_size = 60 

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

## Semiclassical 
print("Setup done")
# Quantum drives (declared here to call for levels in semiclassical)
q_disp_drives = np.linspace(0.015, 0.025, 3)
disp_fig = plt.figure()

gs = gridspec.GridSpec(6, 3)
disp_ax = [None, None, None]
disp_ax[0] = disp_fig.add_subplot(gs[0:2, :])
disp_ax[1] = disp_fig.add_subplot(gs[2:4, :], sharex=disp_ax[0])
disp_ax[2] = disp_fig.add_subplot(gs[4:6, :], sharex=disp_ax[1])
# disp_fig, disp_ax = plt.subplots(2, sharex=True, figsize=(15, 15))

disp_alphas = np.linspace(0, 13, 200)
disp_detrange = np.linspace(0.01, 0.045, 280)
disp_drives = np.array([[disp_drive(A, 
                  cavity_freq, 
                  cavity_freq+d, 
                  -1, 
                  c_q_det, 
                  coupling_strength, 
                  kappa_disp) for d in disp_detrange] for A in disp_alphas])

mble_disp = disp_ax[0].contour(disp_detrange, disp_alphas, disp_drives, 
                               levels=q_disp_drives, 
                               linewidths=1.0,
                               cmap='viridis')
disp_ax[0].set_title('Semiclassical Dispersive', loc='right', 
                fontdict={'fontsize': 16, 'verticalalignment': 'bottom'})
disp_ax[0].set_xlabel('$\omega_c - \omega_d$')
disp_ax[0].set_ylabel('$\left|A\\right|^2$')
if len(q_disp_drives) > 1:
    disp_cbar = plt.colorbar(mble_disp, ax=disp_ax[2], 
                 label='Drive Strength', 
                 orientation='horizontal')
    # Puts chosen set of level lines on the colorbar
    Mark_Colorbar = False
    if Mark_Colorbar:
        cont = ax[0].contour(mble_disp, levels=[q_disp_drives_1, q_disp_drives_2], colors=['r', 'b'])
        cbar.add_lines(cont, erase=False)

    # Wider colorbar lines
    disp_childs = disp_cbar.ax.get_children()
    disp_childs[0].set_linewidths(20)
print("Semiclassical done")
## Quantum
# set parameters for each q_disp_drive
q_disp_detrange = disp_detrange
bishop_parameters = [qo.JaynesCummingsParameters(coupling_strength, matrix_size).det_params(
                              drive_strengths=q_disp_drive, # should have xiC1 for dip @ 6.3 drive
                              drive_cavity_detunings=q_disp_detrange, # Should have * chi0 for dip @ 6.3 drive
                              qubit_cavity_detunings=c_q_det,
                              c_op_params=[kappa_disp],
                              omega_cavity=cavity_freq) for q_disp_drive in q_disp_drives]

# Build system for each q_disp_drive
bishop_systems = [qo.SteadyStateJaynesCummingsModel(*b_params, noisy=True) for b_params in bishop_parameters]
print("\nquantum setup done")

# Plot all absolute cavity fields with colors from semiclassical contours
for sys in enumerate(bishop_systems):
    print('{} / {}'.format(sys[0]+1, len(bishop_systems)))
    disp_ax[1].plot(q_disp_detrange, sys[1].abs_cavity_field(), linewidth=1.0, c=disp_childs[0].get_colors()[sys[0]])

# disp_ax[1].plot(q_disp_detrange, abs_cavity_field_Q(bishop_system)**2, linewidth=1.0)
disp_ax[1].set_title('Quantum Dispersive', loc='right', fontdict={'fontsize': 16, 'verticalalignment': 'bottom'})
disp_ax[1].set_xlabel('$\omega_c - \omega_d$')
disp_ax[1].set_ylabel('$\left | \langle a \\rangle \\right|$')
gs.update(wspace=0.5, hspace=1)
print("field done")

for sys in enumerate(bishop_systems):
    print('{} / {}'.format(sys[0]+1, len(bishop_systems)))
    disp_ax[2].plot(q_disp_detrange, sys[1].correlator(), linewidth=1.0, c=disp_childs[0].get_colors()[sys[0]])
disp_ax[2].set_xlabel('$\omega_c-\omega_d$')
disp_ax[2].set_ylabel('$\langle a \sigma_- \\rangle - \langle a  \\rangle \langle \sigma_- \\rangle$')
disp_ax[2].set_title('Correlation', loc='right',
fontdict={'fontsize': 16, 'verticalalignment':'bottom'})
print("correlator done")
plt.show()
