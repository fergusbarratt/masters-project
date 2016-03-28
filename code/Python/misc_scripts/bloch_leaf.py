import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid', rc={'figure.figsize':(10, 4)})

g = 0.3347
omega_c = 10.5665
omega_q = 8.1831
kappa = 0.0012
gamma = 0.0001

omega_ds = np.arange(10.55, 10.65, 0.0001)
n = np.arange(0.001, 100, 0.05)

n0 = (gamma**2)/(8*g**2)
I_cavs = n/n0

[omega_d, I_cav] = np.meshgrid(omega_ds, I_cavs)

delta_cd = omega_d-omega_c
delta_qd = omega_d-omega_q

phi = 2*delta_cd/kappa
Delta = 2*delta_qd/gamma
C = (2*g**2)/(kappa*gamma)

term = 2*C/(1+Delta**2 + I_cav)
I_d = I_cav * ((1+term)**2 + (phi - Delta*term)**2)
epsilon = np.sqrt(kappa**2*n0*I_d/4)

sp = 0.03
# conts = np.linspace(0.01-sp, 0.01+sp, 10)
conts = np.linspace(0, 1, 100)
plt.contour(omega_ds, n, epsilon, conts, cmap='inferno', linewidths=1)

plt.colorbar()
plt.xlabel('$\\omega_d-\\omega_c$')
plt.ylabel('$A^2$')
plt.title('Semiclassical: Optical Bloch Equations', loc='right')
plt.savefig('disp_optical_bloch_equations.pdf')
plt.show()


