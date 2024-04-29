import matplotlib.pyplot as plt
import numpy as np

dx = 0.1
x = [i*dx for i in range(11)]
y = np.genfromtxt(r'Lab4.txt', delimiter=',')
for i in range(len(y)):
    plt.plot(x, y[i], label = 'Solution')
    plt.legend()
    plt.savefig(f'Graph{i}.png', dpi = 350)
    plt.clf()