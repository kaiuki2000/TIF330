import matplotlib.pyplot as plt
import numpy as np

# Definitions/Reading results
with open('xVec.txt', 'r') as f:
    x = np.array([[float(num) for num in line.split(',')] for line in f][0])
with open('yVec.txt', 'r') as f:
    y = np.array([[float(num) for num in line.split(',')] for line in f][0])

# Plotting result
plt.plot(x, y, label = 'Solution')
plt.title(f"Catenary")
plt.ylabel(f"y")
plt.xlabel("x")
plt.grid(True)
plt.legend()
plt.savefig(f"Catenary.png", dpi = 350)