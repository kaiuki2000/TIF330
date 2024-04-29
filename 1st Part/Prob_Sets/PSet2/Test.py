import numpy as np
from scipy.fftpack import fft2, ifft2

# Parameters
L = 10  # Length of space domain (in meters)
N = 512  # Number of grid points in each dimension
dx = L/N  # Grid spacing (in meters)
dt = 0.01  # Time step (in seconds)
T = 10  # Total time to evolve system (in seconds)

# Initialize wavefunction
x = np.linspace(-L/2, L/2, N)
y = np.linspace(-L/2, L/2, N)
xx, yy = np.meshgrid(x, y)
psi = np.exp(-((xx)**2 + (yy)**2)/2)
kx = np.fft.fftfreq(N, dx) * 2*np.pi
ky = np.fft.fftfreq(N, dx) * 2*np.pi
kxx, kyy = np.meshgrid(kx, ky)
k2 = kxx**2 + kyy**2
V = np.zeros_like(psi)

# Evolution loop
t = 0
while t < T:
    # Kinetic energy term
    psi = np.exp(-1j*k2*dt/2)*psi

    # Potential energy term (V=0)
    psi = np.exp(1j*V*dt)*psi

    # Kinetic energy term
    psi = np.exp(-1j*k2*dt/2)*psi

    # Normalize wavefunction
    psi = psi / np.sqrt(np.sum(np.abs(psi)**2*dx**2))

    # Update time
    t += dt

# Plot final wavefunction
import matplotlib.pyplot as plt
plt.imshow(np.abs(psi)**2, extent=[-L/2, L/2, -L/2, L/2])
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()
plt.savefig("Testing.png")