import matplotlib.pyplot as plt 
import struct

experimentName = 'e1'

infile = open(experimentName + '.txt', 'r')
experimentTxt = infile.readlines()
infile.close()
numberOfIterations = int(experimentTxt[0])
gridSize = int(experimentTxt[1])

u = [0]*gridSize

x = [0]*gridSize
for i in range(gridSize):
    x[i] = 0 + (i + 0.5)*1.0/gridSize

f = open(experimentName + '.dat', "rb")

for frame in range(numberOfIterations):
    for i in range(gridSize):
        u[i] = struct.unpack('d', f.read(8))[0]
    fig, ax = plt.subplots()
    ax.plot(x, u)
    plt.savefig(experimentName + "_" + format(frame, '03d') + '.png', dpi=100)
    plt.close()
    print('Frame', frame, 'is completed.')