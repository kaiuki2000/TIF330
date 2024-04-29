import matplotlib.pyplot as plt
import numpy as np
import struct

experimentName = 'out1'

infile = open(experimentName + '.txt', 'r')
experimentTxt = infile.readlines()
infile.close()
Nx = int(experimentTxt[0])
Ny = int(experimentTxt[1])

print(experimentName, ": Nx = ", Nx, ", Ny = ", Ny)
f = open(experimentName + '.dat', "rb")

for frame in range(10):
    fig, ax = plt.subplots(1,2)

    Ey = np.zeros(shape=(Ny, Nx))
    
    for y in range(Ny):
        for x in range(Nx):
            Ey[Ny - 1 - y, x] = struct.unpack('d', f.read(8))[0]
    
    ax[0].set_title('$E_y$')
    ax[0].set_aspect(1.0)
    ax[0].set_xlabel('x [cm]')
    ax[0].set_ylabel('y [cm]')
    im = ax[0].imshow(Ey, cmap='RdBu', interpolation='none', extent=(0, 1, 0, 1), vmax=2, vmin=-2)
    
    density = np.zeros(shape=(Ny, Nx))

    for y in range(Ny):
        for x in range(Nx):
            density[Ny - 1 - y, x] = struct.unpack('d', f.read(8))[0]
    
    ax[1].set_title('$density$')
    ax[1].set_aspect(1.0)
    ax[1].set_xlabel('x [cm]')
    ax[1].set_ylabel('y [cm]')
    im = ax[1].imshow(density, cmap='YlOrBr', interpolation='none', extent=(0, 1, 0, 1), vmax=5, vmin=0)

    
    plt.savefig(experimentName + "_" + format(frame, '03d') + '.png', dpi=200)
    plt.close()
    print('Frame', frame, 'is completed.')
