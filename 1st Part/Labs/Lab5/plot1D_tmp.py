import matplotlib.pyplot as plt
import struct
f = open('plot1D_tmp.dat', 'rb')
N = struct.unpack('I', f.read(4))[0]
y = [0]*N
for i in range(N):
    y[i] = struct.unpack('d', f.read(8))[0]
fig, ax = plt.subplots()
ax.plot(y)
plt.savefig('ExFig9900.png')
