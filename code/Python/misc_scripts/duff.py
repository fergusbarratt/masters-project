import qutip as qt
import numpy as np
import time as t
import quantumoptics as qo
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='talk',
        style='whitegrid',
        palette='pastel',
        rc={'image.cmap': 'inferno', 'figure.figsize':(9, 9)})


xvec, yvec = np.linspace(-12,12, 100), np.linspace(-12, 12, 100)

c_freq = 10
anh = 0.02
kappa = 1
drive = 15
dets_0 = np.linspace(-10, 10, 200)
dets_1 = -3.47
sys_0 = qo.QuantumDuffingOscillator(drive,
                                    c_freq+dets_0,
                                    c_freq,
                                    anh,
                                    120,
                                    [kappa])

sys_1 = qo.QuantumDuffingOscillator(drive,
                                    c_freq+dets_1,
                                    c_freq,
                                    anh,
                                    120,
                                    [kappa])

fig, axes = plt.subplots(sys_1.length+1, 1)

try:
    axes[0]
except:
    axes = [axes]

axes[-1].plot(dets_0, sys_0.abs_cavity_field())

plt.savefig('''Images/c_freq={c_freq}.anh={anh}.kappa={kappa}:{timestamp}.pdf'''.format(c_freq=c_freq,
       anh=anh,
       kappa=kappa,
       timestamp=t.time()))

for ax in enumerate(axes[:-1]):
    mble = ax[1].contourf(xvec, yvec, sys_0.qps(
      xvec,
      yvec)[ax[0]], 50)

    ax[1].set_title('''$\\omega_c = {}$ $\\omega_d = {}$ $ \\xi = {}$ $ \\chi = {}$'''.format(*sys_0.params[ax[0]]),
       loc='right',
       fontdict={'verticalalignment':'bottom', 'fontsize': 16})

plt.savefig('''Images/c_freq={c_freq}.anh={anh}.kappa={kappa}:{timestamp}.pdf'''.format(c_freq=c_freq,
       anh=anh,
       kappa=kappa,
       timestamp=t.time()))

plt.show()
