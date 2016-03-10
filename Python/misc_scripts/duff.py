import qutip as qt
import numpy as np
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
drive = 6
dets = 0
sys_0 = qo.QuantumDuffingOscillator(drive, 
                                    c_freq+dets, 
                                    c_freq, 
                                    anh, 
                                    100, 
                                    [kappa])

# plt.plot(dets, sys_0.abs_cavity_field())

# plt.savefig('''Images/c_freq={c_freq}.anh={anh}.kappa={kappa}.pdf'''.format(c_freq=c_freq, 
#        anh=anh, 
#        kappa=kappa))

# plt.show()

### Q functions
fig, axes = plt.subplots(sys_0.length, 1)

try:
    axes[0]
except:
    axes = [axes]

for ax in enumerate(axes):
    mble = ax[1].contour(xvec, yvec, sys_0.qps(xvec, yvec)[ax[0]], 50)
    plt.colorbar(mble, ax=ax[1])

# plt.savefig('''Images/Q: c_freq={c_freq}.anh={anh}.kappa={kappa}.pdf'''.format(c_freq=c_freq, 
        # anh=anh, 
        # kappa=kappa))
plt.show()
