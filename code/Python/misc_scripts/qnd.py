import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='talk', style='ticks')

def lorentz(width, strength, center, x):
    return strength/(width**2 + 1j*(x-center)**2)

width = 0.5
center_1 = 1
center_2 = -1 
strength = 0.2*(1/0.9)
x = np.linspace(-5, 5, 200)
plt.plot(x, [lorentz(width, strength, center_1, x) for x in x], label = "\sigma_z = -1")
plt.plot(x, [lorentz(width, strength, center_2, x) for x in x], label="\sigma_z = +1")
plt.legend(loc="best")
plt.ylim(0, 1)
plt.xlabel("Detuning, normalised to dispersive shift")
plt.ylabel("Transmission, arbitrary units")
plt.title('qubit dependent transmission spectrum', loc='right')
sns.despine()
plt.savefig('transmission.pdf')
plt.show()
