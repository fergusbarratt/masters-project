import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='talk', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 12, 8
import numpy as np

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
disp_alphas = np.linspace(3, 6000, 2000)
disp_dets = np.linspace(0.002, 0.0002, 6)
disp_drives = np.asarray([np.array([disp_drive(A, 
                  cavity_freq, 
                  cavity_freq+dd, 
                  -1, 
                  c_q_det, 
                  coupling_strength, 
                  kappa_disp) for A in disp_alphas]) for dd in disp_dets])
for da in enumerate(disp_drives):
    disp_ax[0].plot(da[1], disp_alphas, label="{:0.4f}".format(disp_dets[da[0]]), c=sb.color_palette('viridis')[da[0]], linewidth=1)
disp_ax[0].legend(loc='best')
# disp_ax[0].set_xscale('log')
disp_ax[0].set_xlim([0.22, 0.6])
disp_ax[0].set_ylim([3, 1000])
disp_ax[0].set_yscale('log')
# Set plot parameters
disp_ax[0].set_title('Semiclassical', 
                     loc='right', 
                     fontdict={'fontsize': 12, 
                               'verticalalignment': 'bottom'})
disp_ax[0].set_xlabel('$\\xi$')
disp_ax[0].set_ylabel('$\left|A\\right|^2$')

plt.savefig('{}.pdf'.format(disp_dets[0]))
plt.show()
