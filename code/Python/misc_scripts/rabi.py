import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
sns.set(context='notebook', palette='Paired')

t = np.linspace(0, 10, 100)
plt.plot(t, np.sin(t)**2, label='$\\left| C_1(t) \\right|^2$')
plt.plot(t, np.cos(t)**2, label='$\\left| C_2(t) \\right|^2$')
plt.ylabel('Probability Amplitude')
plt.xlabel('Time')
plt.title('Rabi Flopping', loc='right')
plt.legend(loc=4)
plt.savefig('rabi.pdf')
plt.show()
