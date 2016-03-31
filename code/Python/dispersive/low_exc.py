import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='notebook', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 12, 8

import numpy as np
import qutip as qt

import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

q_disp_drives = np.linspace(0.1, 0.8, 5)
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

disp_alphas = np.linspace(0, 350, 200)
disp_detrange = np.linspace(-0.01, 0.035, 280)
disp_drives = np.array([[disp_drive(A, 
                  cavity_freq, 
                  cavity_freq+d, 
                  -1, 
                  c_q_det, 
                  coupling_strength, 
                  kappa_disp) for d in disp_detrange]
                              for A in disp_alphas])

disp_fig = plt.figure()
gs = gridspec.GridSpec(8, 12)
disp_ax = disp_fig.add_subplot(gs[:, :-1])
cbar_ax = disp_fig.add_subplot(gs[:, -1])

mble_disp = disp_ax.contour(disp_detrange, disp_alphas, disp_drives, 
                        levels=q_disp_drives, 
                        linewidths=1.0,
                        cmap='jet')

disp_ax.set_title('Semiclassical',
                  loc='right')
disp_ax.set_xlabel('$\\omega_c-\\omega_d$')
disp_ax.set_ylabel('$\\left| A \\right|$')
# Wider colorbar lines
disp_cbar = plt.colorbar(mble_disp, cax=cbar_ax, 
             label='Drive Strength')
disp_childs = disp_cbar.ax.get_children()
disp_childs[0].set_linewidths(20)

plt.savefig('low_exc.pdf')
plt.show()
