import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(context='notebook', palette='Paired')
import matplotlib as mpl
import numpy as np
import quantumoptics as qo
import sys

# Dispersive
matrix_size = int(sys.argv[1]) 

c_q_det = -(8.1831-10.5665)
kappa_disp = 0.0024
gamma = 0
coupling_strength = 0.3347
cavity_freq = 10.5665
sigmaz=-1

q_drives = np.append(np.linspace(0, 0.2, 150), np.linspace(0.05, 0.4, 50))

disp_alphas = np.linspace(0, 300, 100)
disp_dets = [0.002] # np.linspace(0.002, 0.0002, 1)

def disp_drive(A, omega_qubit, omega_cavity, omega_drive, sigmaz, g, kappa):
  '''self-consistently determine xi from A'''
  det = omega_cavity-omega_qubit
  def chi(A):
    return sigmaz*(g**2)/np.sqrt(2*g**2*(A**2+sigmaz)+det**2)

  return np.sqrt(
          A**2 * (1/omega_cavity**2)*(
            (omega_drive**2-(omega_cavity-chi(A))**2)**2 +\
                    kappa_disp**2*omega_drive**2))

# Make all the axes
disp_fig = plt.figure()
gs = gridspec.GridSpec(12, 12)
disp_ax = [None, None, None, None]
disp_ax[0] = disp_fig.add_subplot(gs[0:3, :])
disp_ax[1] = disp_fig.add_subplot(gs[3:6, :], sharex=disp_ax[0])
disp_ax[2] = disp_fig.add_subplot(gs[6:9, :], sharex=disp_ax[1])
disp_ax[3] = disp_fig.add_subplot(gs[9:12, :], sharex=disp_ax[2])

## Semiclassical 
print('semiclassical', flush=True)

print('{} curves, {} xvals'.format(len(disp_dets), len(disp_alphas), flush=True))

disp_drives = np.asarray([np.array([disp_drive(A=A, 
                  omega_qubit=(cavity_freq+c_q_det), 
                  omega_cavity=cavity_freq, 
                  omega_drive=cavity_freq + dd, 
                  sigmaz=sigmaz, 
                  g=coupling_strength, 
                  kappa=kappa_disp) for A in disp_alphas]) for dd in disp_dets])

for da in enumerate(disp_drives):
    print('|', sep='', end='', flush=True)
    disp_ax[0].plot(da[1], 
                    disp_alphas, 
                    label="{:0.4f}".format(disp_dets[da[0]]), 
                    linewidth=1)

print('', flush=True)


###QUANTUM

print('quantum')

carms = [qo.JaynesCummingsSystem(drive_range=q_drives,
                                omega_qubit_range=(cavity_freq+c_q_det),
                                omega_cavity_range=cavity_freq,
                                omega_drive_range=cavity_freq + dd,
                                c_op_params=[kappa_disp, gamma],
                                coupling=coupling_strength,
                                N_field_levels=matrix_size,
                                noisy=True) for dd in disp_dets] 

for carm in enumerate(carms):
    disp_ax[1].plot(q_drives, 
                    carm[1].abs_cavity_field()**2, 
                    label="{:0.4f}".format(disp_dets[carm[0]]), 
                    linewidth=1)
    disp_ax[2].plot(q_drives,
                    carm[1].g2(),
                    label="{:0.4f}".format(disp_dets[carm[0]]), 
                    linewidth=1)
    disp_ax[3].plot(q_drives,
                    carm[1].correlator(),
                    label="{:0.4f}".format(disp_dets[carm[0]]), 
                    linewidth=1)
    
# Set plot parameters

disp_ax[0].set_xlim([min(q_drives), max(q_drives)])
disp_ax[0].set_ylim([min(disp_alphas), max(disp_alphas)/2])

disp_ax[0].set_title('Semiclassical', 
                     loc='right', 
                     fontdict={'fontsize': 16, 
                               'verticalalignment': 'bottom'})
disp_ax[1].set_title('Quantum', 
                     loc='right', 
                     fontdict={'fontsize': 16, 
                               'verticalalignment': 'bottom'})

disp_ax[2].set_title('Quantum', 
                     loc='right', 
                     fontdict={'fontsize': 16, 
                               'verticalalignment': 'bottom'})


disp_ax[0].set_ylabel('$\left|A\\right|^2$')
disp_ax[1].set_ylabel('$\left|A\\right|^2$')
disp_ax[2].set_ylabel('$\\langle a \\sigma_- \\rangle - \\langle a \\rangle \\langle \\sigma_-\\rangle$')
disp_ax[3].set_ylabel('$g^{(2)}(0)$')

disp_ax[3].set_xlabel('$\\xi$')

disp_ax[0].legend(loc='best')
disp_ax[1].legend(loc='best')
disp_ax[2].legend(loc='best')

plt.tight_layout()

plt.savefig('{}.pdf'.format(disp_dets[0]))

plt.show()
