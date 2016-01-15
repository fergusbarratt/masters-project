from quantumoptics import *

def norm_e(num):
    return num/np.sqrt(2)
def Ghz(num):
    return num*10**9
def norm_det(num, detuning):
    chi = coupling_strength**2/detuning
    return num/abs(chi)

matrix_size = 40
coupling_strength = Ghz(0.2)*2*np.pi

carmichael_parameters = JaynesCummingsParameters(coupling_strength, matrix_size).det_params(
                                  5, np.linspace(-3, 3, 20), 0, [1, 1])
carmichael_system = SteadyStateJaynesCummingsModel(
    *carmichael_parameters)

c_d_det = Ghz(-1)*2*np.pi
kappa = 2*np.pi*Ghz(0.001)

bishop_parameters = JaynesCummingsParameters(coupling_strength, matrix_size).det_params(
                              norm_e(1000), norm_det(0.5, c_d_det), c_d_det, [kappa, 1], Ghz(9.07))

bishop_system = SteadyStateJaynesCummingsModel(
    *bishop_parameters)

bishop_system.draw_qps()
