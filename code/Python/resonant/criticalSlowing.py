import matplotlib.pyplot as plt
from matplotlib import gridspec
import seaborn as sb
sb.set(style='whitegrid', context='notebook', rc={'image.cmap': 'viridis'})
import matplotlib as mpl
mpl.rcParams['figure.figsize'] = 17, 5

import numpy as np
import qutip as qt
from scipy.special import jv
import scipy.integrate as integr

import time

import importlib
import quantumoptics as qo
qo = importlib.reload(qo)

### Parameters
matrix_size = 5

# Dispersive
c_q_det = (8.1831-10.5665)
kappa_disp = 0.0024
gamma = 0
coupling_strength = 0.3347
cavity_freq = 10.5665
sigmaz=-1
q_disp_drives = np.linspace(0.02, 0.05, 2)
q_disp_detrange = np.linspace(-0.1, 0.1, 100)

xiC1= (abs(c_q_det)*kappa_disp)**(3/2)/(3**(3/4)*coupling_strength**2)
chi0=coupling_strength**2/abs(c_q_det)

def dxidA(A, drive_freq):
    return np.sqrt(((cavity_freq-sigmaz*coupling_strength**2/np.sqrt(coupling_strength**2*(sigmaz+A**2)*2+c_q_det**2))**2-drive_freq**2)**2+drive_freq**2*kappa_disp**2)/cavity_freq+(A**2*sigmaz*coupling_strength**4*(cavity_freq-sigmaz*coupling_strength**2/np.sqrt(coupling_strength**2 *(sigmaz+A**2)*2+c_q_det**2))/np.sqrt(((cavity_freq-sigmaz*coupling_strength**2/np.sqrt(coupling_strength**2*(sigmaz+A**2)*2+c_q_det**2))**2-cavity_freq**2)**2+(cavity_freq**2) *kappa_disp**2)/(coupling_strength**2*((sigmaz+A**2)*2.0+c_q_det**2)**(3.0/2.0))*((cavity_freq-sigmaz*coupling_strength**2/np.sqrt(coupling_strength**2*(sigmaz+A**2)*2.0+c_q_det**2))**2-cavity_freq**2)*4.0)/cavity_freq

diff_alphas = np.linspace(0.1, 100, 250)
diff_freqs = cavity_freq+(30*q_disp_detrange*chi0+0.03)
diff_drive_diffs = np.array([[dxidA(alpha, freq) for alpha in diff_alphas] for freq in diff_freqs])
plt.contour(diff_freqs, diff_alphas/xiC1, diff_drive_diffs.T, levels=np.linspace(-0.02, 0.02, 10), linewidths=0.8, figsize=(15, 5), cmap='jet')
plt.xlabel('Frequency')
plt.ylabel('Alpha')
plt.yscale('log')
plt.xscale('log')
plt.title('Zero dxidA contour')
plt.vlines(cavity_freq, 0, 100000)
plt.show()
