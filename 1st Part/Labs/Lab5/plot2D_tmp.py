import matplotlib.pyplot as plt
import numpy as np
import struct
f = open('plot2D_tmp.dat', 'rb')
Ny = struct.unpack('I', f.read(4))[0]
Nx = struct.unpack('I', f.read(4))[0]
v = np.empty(shape=(Ny, Nx))
for y in range(Ny):
    for x in range(Nx):
        v[y][x] = struct.unpack('d', f.read(8))[0]
fig, ax = plt.subplots()
vAbsMax = np.maximum(abs(np.amin(v)), abs(np.amax(v)))
plot = ax.imshow(v, interpolation='none', extent=[0,1,0,1], origin='lower', cmap='plasma', vmin = 0, vmax = vAbsMax)
ax.set_aspect('equal')
fig.colorbar(plot, ax=ax, location='right')
plt.savefig('xpx9900.png')
