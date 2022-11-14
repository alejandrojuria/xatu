# Script to plot bands from the text file generated from the 
# system.cpp code

import matplotlib.pyplot as plt
import numpy as np

# --------------------- Plot bands ---------------------
file = open("eigval.out", "r")
lines = file.readlines()
kpoints = []
energies = []
for line in lines:
    line = line.split()
    kEnergy = np.array([float(i) for i in line])
    energies.append(kEnergy)

energies = np.array(energies)
kpoints = np.array(kpoints)

print(f"Gap: {energies[99, 4] - energies[99, 3]} eV")

for n, band in enumerate(np.transpose(energies)):
    if n == 13 or n == 12:
        plt.plot(band, 'r-')
    else:
        plt.plot(band, 'g-')

nk = 100
xticks = [0, nk - 1, 2*nk - 2, 3*nk - 3]
xlabels = [r"$\Gamma$", r"K", "M", r"$\Gamma$"]
plt.xticks(labels=xlabels, ticks=xticks)
plt.show()