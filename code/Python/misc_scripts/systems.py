from quantumoptics import *
import matplotlib.pyplot as plt
import cmath
from scipy.optimize import fixed_point, newton

def Ghz(num):
    return num*10**9

c_q_det = -2*np.pi*Ghz(1)
kappa = 2*np.pi*Ghz(0.001)
gamma = 2*np.pi*Ghz(0.000000)
coupling_strength = 2*np.pi*Ghz(0.2)
matrix_size = 40

def norm_e(num):
    return num*(kappa/np.sqrt(2))

def norm_det(num, cavity_qubit_detuning):
    return num*abs(coupling_strength**2/cavity_qubit_detuning)

carmichael_parameters = JaynesCummingsParameters(10, matrix_size).det_params(
                            np.linspace(4, 6, 10), 0, 0, [1, 1])
carmichael_system = SteadyStateJaynesCummingsModel(
    *carmichael_parameters)

bishop_parameters = JaynesCummingsParameters(coupling_strength, matrix_size).det_params(
                              drive_strengths=[norm_e(e) for e in np.linspace(100, 100, 1)],
                              drive_cavity_detunings=[norm_det(det, c_q_det) for det in np.linspace(0, 1.5, 20)],
                              qubit_cavity_detunings=c_q_det,
                              c_op_params=[kappa, gamma],
                              omega_cavity=Ghz(9.07))

bishop_system = SteadyStateJaynesCummingsModel(
    *bishop_parameters)

def absAs(sys, xvec=np.linspace(-5, 5, 100), yvec=np.linspace(-5, 5, 100)):
    Qs = sys.qps(xvec, yvec)
    return [abs(np.sum(xvec.T.dot(Q))+1j*np.sum(Q.dot(yvec))) for Q in Qs]

plt.plot(np.linspace(0, 1.5, 20), absAs(bishop_system))
plt.show()
