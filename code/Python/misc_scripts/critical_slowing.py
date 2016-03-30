import numpy as np
import matplotlib.pyplot as plt
import quantumoptics as qo
import seaborn as sns
sns.set(context='talk', rc={'figure.figsize':(10, 5)})

t = np.linspace(0, 8, 200)
drive_1 = 10
drive_2 = 4.5 
det_1 = det_2 = 0

carm_1 = qo.TimeDependentJaynesCummingsModel(drive_range=drive_1,
                                             omega_qubit_range=10,
                                             omega_cavity_range=10,
                                             omega_drive_range=10,
                                             c_op_params=[1, 1],
                                             g=10,
                                             N=40,
                                             tlist=t)

carm_2 = qo.TimeDependentJaynesCummingsModel(drive_range=drive_2,
                                             omega_qubit_range=10,
                                             omega_cavity_range=10,
                                             omega_drive_range=10,
                                             c_op_params=[1, 1],
                                             g=10,
                                             N=40,
                                             tlist=t)
y_1 = carm_1.mesolve([carm_1.num])
y_2 = carm_2.mesolve([carm_2.num])
plt.plot(t, y_1.expect[0])
plt.plot(t, y_2.expect[0])
plt.legend(['drive = {}, det = {}'.format(drive_1, det_1), 'drive = {}, det = {}'.format(drive_2, det_2)])
plt.xlabel('time')
plt.ylabel('$A^2$')
plt.title('Critical Slowing', loc='right', fontdict={'fontsize': 15})
plt.savefig('critical_slowing.pdf')
plt.show()
