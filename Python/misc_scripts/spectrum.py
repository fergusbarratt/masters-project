import matplotlib
matplotlib.use('TkAgg')
import numpy as np
import matplotlib.pyplot as plt


def spectrum(W, W0, gamma, E, g, kappa):
    np.asarray(W)
    Y1 = (abs(E) / g)**4 * (((2 * (0.5 * (kappa + gamma / 2))**3 / np.pi) / ((0.5 * (kappa + (gamma / 2)))**2 + (W - W0 + g) ** 2)**2) + 2 * (0.5 * ((kappa+(gamma/2)))**3/np.pi)/((0.5 * (kappa + gamma / 2))**2 + (W-W0-g)**2)**2)

    Y2 = 0.5*((0.5*(kappa+gamma/2)/np.pi)/(0.25*(kappa+gamma/2)**2+(W-W0+g)**2) + (0.5*(kappa+gamma/2)/np.pi)/(0.25*(kappa+gamma/2)**2+(W - W0 - g)**2));

    return (W, Y1, W, Y2)

plt.plot(*spectrum(np.linspace(80, 120, 100), 100, 5, 5, 10, 1))
plt.show()
