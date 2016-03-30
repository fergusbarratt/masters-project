'''Produce dispersive plot'''

import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='talk', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 12, 8

import numpy as np
import qutip as qt

import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

### Parameters
matrix_size = 80 

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

  return np.sqrt(
          A**2 * (1/omega_cavity**2)*(
            (omega_drive**2-(omega_cavity-chi(A))**2)**2 + kappa_disp**2*omega_drive**2))

## Semiclassical 

print("Setup done")

# Quantum drives (declared here to call for levels in semiclassical)
q_disp_drives = np.linspace(0.015, 0.025, 5)

# Make all the axes
disp_fig = plt.figure()
gs = gridspec.GridSpec(8, 12)
disp_ax = [None, None, None, None]
disp_ax[0] = disp_fig.add_subplot(gs[0:2, :-1])
disp_ax[1] = disp_fig.add_subplot(gs[2:4, :-1], sharex=disp_ax[0])
disp_ax[2] = disp_fig.add_subplot(gs[4:6, :-1], sharex=disp_ax[1])
disp_ax[3] = disp_fig.add_subplot(gs[6:8, :-1], sharex=disp_ax[2])
cbar_ax = disp_fig.add_subplot(gs[0:8, -1:])

# determine drives from A
disp_alphas = np.linspace(0, 13, 200)
disp_detrange = np.linspace(0.01, 0.045, 280)
disp_drives = np.array([[disp_drive(A, 
                  cavity_freq, 
                  cavity_freq+d, 
                  -1, 
                  c_q_det, 
                  coupling_strength, 
                  kappa_disp) for d in disp_detrange]
                              for A in disp_alphas])

# Plot contours of constant drive at chosen levels
mble_disp = disp_ax[0].contour(disp_detrange, disp_alphas, disp_drives, 
                               levels=q_disp_drives, 
                               linewidths=1.0,
                               cmap='jet')

# Set plot parameters
disp_ax[0].set_title('Semiclassical', 
                     loc='right', 
                     fontdict={'fontsize': 12, 
                               'verticalalignment': 'bottom'})
disp_ax[0].set_xlabel('$\omega_c - \omega_d$')
disp_ax[0].set_ylabel('$\left|A\\right|$')


if len(q_disp_drives) > 1:
    disp_cbar = plt.colorbar(mble_disp, cax=cbar_ax, 
                 label='Drive Strength')
                 

    # Wider colorbar lines
    disp_childs = disp_cbar.ax.get_children()
    disp_childs[0].set_linewidths(20)
    
print("Semiclassical done")

## Quantum
# set parameters for each q_disp_drive
q_disp_detrange = disp_detrange
bishop_parameters = [qo.JaynesCummingsParameters().det_params(
                            drive_strengths=q_disp_drive, 
                            drive_cavity_detunings=q_disp_detrange,
                            qubit_cavity_detunings=c_q_det,
                            c_op_params=[kappa_disp],
                            omega_cavity=cavity_freq,
                            g=coupling_strength,
                            N=matrix_size) 
                     for q_disp_drive in q_disp_drives]

# Build system for each q_disp_drive
bishop_systems = [qo.JaynesCummingsSystem(*b_params)
                  for b_params in bishop_parameters]

print("\nquantum setup done")
field_data = np.empty((len(bishop_systems), len(q_disp_detrange)))
# Plot all absolute cavity fields with colors from semiclassical contours
for sys in enumerate(bishop_systems):
    cav_field = sys[1].abs_cavity_field()
    field_data[sys[0]] = cav_field
    print('{} / {}'.format(sys[0]+1, len(bishop_systems)))
    disp_ax[1].plot(q_disp_detrange, cav_field, 
                linewidth=1.0, 
                c=disp_childs[0].get_colors()[sys[0]])

# disp_ax[1].plot(q_disp_detrange, abs_cavity_field_Q(bishop_system)**2, linewidth=1.0)
disp_ax[1].set_title('Quantum', 
                     loc='right', 
                     fontdict={'fontsize': 12, 
                     'verticalalignment': 'bottom'})
disp_ax[1].set_xlabel('$\omega_c - \omega_d$')
disp_ax[1].set_ylabel('$\left | \langle a \\rangle \\right|$')

print("field done")

correlator_data = np.empty((len(bishop_systems), len(q_disp_detrange)))

# Plot correlators for all systems
for sys in enumerate(bishop_systems):
    corr = sys[1].correlator()
    correlator_data[sys[0]] = corr
    print('{} / {}'.format(sys[0]+1, len(bishop_systems)))
    disp_ax[2].plot(q_disp_detrange, corr, 
            linewidth=1.0, 
            c=disp_childs[0].get_colors()[sys[0]])

disp_ax[2].set_xlabel('$\omega_c-\omega_d$')

disp_ax[2].set_ylabel('$\langle a \sigma_- \\rangle - \langle a  \\rangle \langle \sigma_- \\rangle$')

disp_ax[2].set_title('Quantum', 
                     loc='right',
                     fontdict={'fontsize': 12,
                     'verticalalignment':'bottom'})

print("correlator done")
print("starting g2")

for sys in enumerate(bishop_systems):
    g2 = sys[1].g2()
    print('{} / {}'.format(sys[0]+1, len(bishop_systems)))
    disp_ax[3].plot(q_disp_detrange, g2,
            linewidth=1.0,
            c=disp_childs[0].get_colors()[sys[0]])

disp_ax[3].set_xlabel('$\omega_c-\omega_d$')
disp_ax[3].set_ylabel('$g^{(2)}(0)$')
disp_ax[3].set_title('Quantum',
                     loc='right',
                     fontdict={'fontsize' : 12,
                    'verticalalignment':'bottom'})

print("g2 done")

# Update grid spacings
gs.update(wspace=1, hspace=1)
np.save('correlator_data.npy', correlator_data)
np.save('field_data.npy', field_data)
disp_fig.suptitle('Dispersive Bistability', x=0.2, y=0.93,
        verticalalignment='bottom', fontsize=16)
disp_fig.savefig('dispersive.pdf')
plt.show()
