import qutip as qt
import numpy as np
import quantumoptics as qo

sys = qo.QuantumDuffingOscillator(5, 1, 40, [0])

print(list(sys.hamiltonian()))

input()
