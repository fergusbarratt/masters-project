import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='ticks', context='talk', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 12, 8

import numpy as np
import qutip as qt

import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

q_disp_drives = np.linspace(0.015, 0.025, 5)

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

# Make all the axes
disp_fig = plt.figure()
gs = gridspec.GridSpec(8, 12)
disp_ax = [None, None, None]
disp_ax[0] = disp_fig.add_subplot(gs[0:7, :])

## Semiclassical 
# determine drives from A
disp_alphas = np.linspace(0, 100, 1000)
disp_det = 0.005
disp_drives = np.array([disp_drive(A, 
                  cavity_freq, 
                  cavity_freq+disp_det, 
                  -1, 
                  c_q_det, 
                  coupling_strength, 
                  kappa_disp) for A in disp_alphas])

disp_ax[0].plot(disp_drives, disp_alphas)
# disp_ax[0].set_xscale('log')
disp_ax[0].set_xlim([0, 0.5])
# disp_ax[0].set_yscale('log')
# Set plot parameters
disp_ax[0].set_title('Semiclassical', 
                     loc='right', 
                     fontdict={'fontsize': 12, 
                               'verticalalignment': 'bottom'})
disp_ax[0].set_xlabel('$\\xi$')
disp_ax[0].set_ylabel('$\left|A\\right|^2$')

sb.despine(disp_fig)
plt.savefig('scurve|det:{}.pdf'.format(disp_det))
plt.show()
