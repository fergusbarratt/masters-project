from quantumoptics import *
matrix_size = 40
coupling_strength = 10

carmichael_parameters = JaynesCummingsParameters(coupling_strength, matrix_size).det_params(
                                  5, np.linspace(-3, 3, 20), 0, [1, 1])
carmichael_system = SteadyStateJaynesCummingsModel(
    *carmichael_parameters)

bishop_parameters = JaynesCummingsParameters(coupling_strength, matrix_size).det_params(
                              5, np.linspace(-3, 3, 20), 10000, [1, 1])
bishop_system = SteadyStateJaynesCummingsModel(
    *bishop_parameters)

bishop_system.plot_exp()
