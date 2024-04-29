import matplotlib.pyplot as plt
import numpy as np
plt.style.use('fast')


# Definitions/Reading results
wVec       = np.genfromtxt(r'wVec.txt', delimiter=',')
Psi        = np.genfromtxt(r'PsiXY.txt', delimiter=',')
tMax       = 100.
dt         = 0.1
t          = [i*dt for i in range(int(tMax/dt))]

plt.plot(t, [Psi[i][0] for i in range(0, len(Psi))], label = "Real part")
plt.plot(t, [Psi[i][1] for i in range(0, len(Psi))], label = "Imaginary part")

plt.xlabel(r"time: $t$")
plt.ylabel(r"$\psi(0.1, 0.0)(t)$")
plt.title(r"$\psi(t)$ at $(x, y) = (0.1, 0.0)$")
plt.grid(True)
plt.xlim(0., 100.)
plt.legend()
plt.tight_layout()
plt.savefig(f'Wavefunction(dt={dt},T={tMax}).png', dpi = 500)