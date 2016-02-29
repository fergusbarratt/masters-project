import numpy as np
import qutip as qt
import matplotlib.pyplot as plt
import quantumoptics as qo
import seaborn as sns
sns.set(context='paper', style='white')

q_res_drives = 10
q_res_detrange = 0.01 
t_list = np.linspace(0, 10, 20)
c_params = qo.JaynesCummingsParameters(10, 10).det_params(
              drive_strengths=q_res_drives,
              drive_cavity_detunings=q_res_detrange,
              qubit_cavity_detunings=0,
              c_op_params=[1, 1],
              omega_cavity=10)

carmichael_system = qo.TimeDependentJaynesCummingsModel(*c_params, t_list)
carmichael_system.trajectory(draw=True)
