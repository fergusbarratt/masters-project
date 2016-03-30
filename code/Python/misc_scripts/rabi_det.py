import matplotlib.pyplot as plt 
import numpy as np
import seaborn as sns 
sns.set(palette='Paired', context='notebook')

delta = 1
f = 2
theta = np.sqrt(delta**2 + f**2)
t = np.linspace(0, 10, 200)
plt.plot(t, abs(-1j*f/theta*np.sin(theta*t))**2, label="$\\left|C_0(t)\\right|^2$")
plt.plot(t, abs((np.cos(theta*t) - 1j*delta/theta*np.sin(theta*t)))**2, label="$\\left|C_1(t)\\right|^2$")
plt.title('Rabi Flopping: with detuning', loc='right')
plt.xlabel('Time')
plt.ylabel('probability amplitude')
plt.legend(loc='best')
plt.savefig('rabi_detuned.pdf')
plt.show()
