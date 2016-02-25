import quantumoptics as qo
import numpy as np

matrix_size = 20
q_res_detrange = 0
q_res_drive=np.linspace(4.5, 5.4, 30)

params = qo.JaynesCummingsParameters(10, matrix_size).det_params( drive_strengths=q_res_drive, drive_cavity_detunings=q_res_detrange, qubit_cavity_detunings=0, c_op_params=[1, 1], omega_cavity=10)

carm_sys = qo.SteadyStateJaynesCummingsModel(*params)

carm_sys.draw_qps(save=True, form='gif', xvec=np.linspace(-8, 8, 100), yvec=np.linspace(-10, 4, 100), plottype='cf')
