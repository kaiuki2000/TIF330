import matplotlib.pyplot as plt
import numpy as np
import imageio.v2 as imageio
from pathlib import Path

SIZE = 2048

# Definitions
h         = 2/(2 * float(SIZE) - 2); tau = 4/(SIZE**2)
x_E       = [-1 + (j * 2*h) for j in range(SIZE)]
x_B       = [(-1 + h) + (j * 2*h) for j in range(SIZE - 1)]
Data      = {}

ids       = ['E_y_1', 'B_z_1',
             'E_y_2', 'B_z_2',
             'E_y_3', 'B_z_3', 'TimeVec']

cases     = ['./Periodic_BC', './PEC_BC', './Perfect_BC']

filenames = ['./E_Periodic_BC.txt', './B_Periodic_BC.txt',
             './E_PEC_BC.txt'     , './B_PEC_BC.txt',
             './E_Perfect_BC.txt' , './B_Perfect_BC.txt', 'timeVec.txt']

for folder in cases:
    Path(f"{folder}").mkdir(parents=True, exist_ok=True)

for file, id in zip(filenames, ids):
    with open(file, 'r') as f:
        Data[id] = [[float(num) for num in line.split(',')] for line in f]

for ind in range(1,4):
    for i in range(0, len(Data[f'E_y_{ind}'])):
        if((i % 40 == 0) or (i == (len(Data[f'E_y_{ind}']) - 1))):
            plt.plot(x_E, Data[f'E_y_{ind}'][i], label = 'E_y')
            plt.plot(x_B, Data[f'B_z_{ind}'][i], label = 'B_z')
            plt.title(f"E_y(x, t = {Data['TimeVec'][0][i]:.5f}) and B_z(x, t = {Data['TimeVec'][0][i]:.5f})")
            plt.ylabel(f"Field value")
            plt.xlabel("Position (x)")
            plt.ylim(-1, 1)
            plt.legend()
            plt.savefig(f"./{cases[(ind - 1)]}/Img{int(i/40)}.png", dpi = 100)
            plt.savefig(f"./{cases[(ind - 1)]}/Img{int(i/40)}.pdf", dpi = 100)
            plt.clf()

    # Small piece of code to join all the 127 images into a gif, for easier visualization. Takes some minutes to run.
    filenames = [f'./{cases[(ind - 1)]}/Img{n}.png' for n in range(0, int(i/40) + 1)]
    with imageio.get_writer(f'./{cases[(ind - 1)]}/Movie_E+B.gif', mode='I') as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)