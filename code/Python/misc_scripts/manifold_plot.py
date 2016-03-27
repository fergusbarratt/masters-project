import numpy as np
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

x = np.linspace(-10, 10, 300)
y = np.linspace(-10, 10, 300)
x, y = np.meshgrid(x, y)
z = -4*x**3 - 2*y*x

ax.plot_surface(x, y, z)
plt.show()
