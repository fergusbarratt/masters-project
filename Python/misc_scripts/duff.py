import qutip as qt
import numpy as np
import quantumoptics as qo
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='talk')

xvec, yvec = np.linspace(-5, 5, 100), np.linspace(-5, 5, 100)

sys_0 = qo.QuantumDuffingOscillator(6, 1, 40, [])

plt.contour(sys_0.qps(xvec, yvec)[0])
plt.show()


#jc_sys_1 = qo.JaynesCummingsSystem([4.5, 5.5], 10, 10, 10, [1, 1], 10, 40)

#xvec, yvec = np.linspace(-10, 10, 100), np.linspace(-10, 10, 100)

#fig, ax = plt.subplots(2)
#ax[0].contourf(xvec, yvec, jc_sys_1.qps(xvec, yvec, tr='cavity')[0], 100, cmap='viridis')
#ax[1].contourf(xvec, yvec, jc_sys_1.qps(xvec, yvec, tr='cavity')[1], 100, cmap='viridis')
#plt.show()
