import quantumoptics as qo
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(style='whitegrid', rc={'figure.figsize': (5, 5)})
import numpy as np

matrix_size = 20
q_res_detrange = 0
q_res_drive = np.linspace(4.5, 5.4, 3)
params = qo.JaynesCummingsParameters(10, matrix_size).det_params( drive_strengths=q_res_drive, drive_cavity_detunings=q_res_detrange, qubit_cavity_detunings=0, c_op_params=[1, 1], omega_cavity=10)

carm_sys = qo.SteadyStateJaynesCummingsModel(*params)

xvec = np.linspace(-10, 10, 100)
yvec = np.linspace(-10, 10, 100)
qfuncs = carm_sys.qps(np.linspace(-10, 10, 100), np.linspace(-10, 10, 100))

fig, ax = plt.subplots(len(qfuncs), 1, sharex=True)

for qfunc in enumerate(qfuncs):
    ax[qfunc[0]].contour(xvec, yvec, qfunc[1])

plt.show()
