import matplotlib.pyplot as plt
import numpy as np

with open('drum_sol.txt', 'r') as f:
    l = [[float(num) for num in line.split(',')] for line in f]

# Definitions
h = 1/100; tau = 1e-5

t = [i*tau for i in range(1000)]
x = [-1 + (j * 2*h) for j in range(200)]

T, X = np.meshgrid(t, x)
U = np.array(l)

ax = plt.axes(projection='3d')
ax.plot_surface(T, X, U.T, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title("Solution (E_y(x, t))")
ax.set_xlabel("Time (t)")
ax.set_ylabel("Position (x)")

plt.savefig(f"drum_Plot_h={h:.5f}_tau={tau:.5f}.png", dpi = 350)