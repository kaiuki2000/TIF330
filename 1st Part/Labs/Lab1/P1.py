import matplotlib.pyplot as plt
import numpy as np

with open('Solution.txt', 'r') as f:
    l = [[float(num) for num in line.split(',')] for line in f]

# Definitions
h = 0.1; tau = 0.00005

t = [i*tau for i in range(int(10/tau))] # This isn't very consistent, but I'll leave it as it is, for now.
x = [j*h for j in range(int(1/h))]

T, X = np.meshgrid(t, x)
U = np.array(l)

ax = plt.axes(projection='3d')
ax.plot_surface(T, X, U.T, rstride=1, cstride=1,
                cmap='viridis', edgecolor='none')
ax.set_title("Solution (u(x, t))")
ax.set_xlabel("Time (t)")
ax.set_ylabel("Position (x)")

plt.savefig(f"Solution_Plot_h={h}_tau={tau}.png", dpi = 350)