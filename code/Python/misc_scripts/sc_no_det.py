import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(context='notebook')

pr = sns.xkcd_rgb["pale red"]
db = sns.xkcd_rgb["denim blue"]
dg = sns.xkcd_rgb["leaf green"]
dy = sns.xkcd_rgb["dark yellow"]
p = sns.xkcd_rgb["purple"]


def drive(alpha, kappa, g):
    return np.abs(1j*alpha*(kappa-1j*(g/(2*np.abs(alpha)))))
def drive_min(alpha, kappa, g):
    return np.abs(1j*alpha*(kappa+1j*(g/(2*np.abs(alpha)))))

g, kappa = 10, 1

# abss = 5
abss = np.linspace(0, 5, 200) 
phas = np.pi/2
# phas = np.linspace(0, np.pi)
sm_drives = np.linspace(0, 5, 200)

alphas = abss*np.exp(1j*phas)
drives = [drive(al, kappa, g) for al in alphas]
drives_min = [drive_min(al, kappa, g) for al in alphas]

fig, ax = plt.subplots(2, 1, sharex=True)

try:
    len(ax)
except:
    ax = [ax]

ax[0].plot(sm_drives, np.sqrt(1-(2*sm_drives/g)**2), c=db)
bl = ax[0].plot(drives, np.zeros_like(drives), c=db, label='inversion')

ax[0].plot(drives, np.abs(alphas), c=pr)
re = ax[0].plot(sm_drives, np.zeros_like(sm_drives), c=pr, label='field')

ax[0].legend(loc='best')

ax[0].set_title('resonant qubit inversion and cavity field amplitude', loc='right')
ax[0].set_ylabel('$\\alpha, \zeta$')

ax[0].set_ylim(-0.5, 5)
ax[0].set_xlim(0, 7)

phas_ = np.linspace(-2*np.pi, 2*np.pi)
abss_ = 5

alphas_ = abss*np.exp(1j*phas)
drives = np.asarray([drive(al, kappa, g) for al in alphas])

ax[1].plot(drives, np.angle(alphas_), c=dy)
ax[1].plot(drives, -np.angle(alphas_), c=dg)
ax[1].plot(sm_drives, np.zeros_like(sm_drives)+0.03, c=dy)
ax[1].plot(sm_drives, np.zeros_like(sm_drives)-0.03, c=dg)
ax[1].set_title('field phase', loc='right')
ax[1].set_ylabel('$\\theta = arg(\\alpha)$')
ax[1].set_xlabel('drive $\\xi$')
plt.savefig('sc_no_det.pdf')
plt.show()
