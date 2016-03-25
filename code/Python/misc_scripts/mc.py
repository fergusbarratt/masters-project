import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import quantumoptics as qo
import seaborn as sns
sns.set()
N = 10
c_q_det = (8.1831-10.5665)
kappa_disp = 0.0024
gamma = 0
coupling_strength = 0.3347
cavity_freq = 10.5665
sigmaz=-1

bis_sys = qo.JaynesCummingsParameters(coupling_strength, N).det_params(
        drive_strengths=0.020,
        drive_cavity_detunings=0.028,
        qubit_cavity_detunings=c_q_det,
        c_op_params=[kappa_disp],
        omega_cavity=cavity_freq)

bishop_system = qo.TimeDependentJaynesCummingsModel(*bis_sys, np.linspace(0, 100, 100), noisy=True)
ode_res = bishop_system.trajectory(draw=True)
plt.show()
