# Script to plot bands from the text file generated from the 
# system.cpp code

import matplotlib.pyplot as plt
import numpy as np

# --------------------- Plot bands ---------------------
file = open("bands_from_file", "r")
lines = file.readlines()
kpoints = []
energies = []
for line in lines:
    line = line.split('\t')
    kpoints.append(float(line[0]))
    kEnergy = np.array([float(i) for i in line[1:-1]])
    energies.append(kEnergy)

energies = np.array(energies)
kpoints = np.array(kpoints)

for n, band in enumerate(np.transpose(energies)):
    if n == 13 or n == 12:
        plt.plot(kpoints, band, 'r-')
    else:
        plt.plot(kpoints, band, 'g-')

nk = 50
xticks = [0, nk, 2*nk, 3*nk + 1]
xlabels = [r"$\Gamma$", r"K", "M", r"$\Gamma$"]
plt.xticks(labels=xlabels, ticks=xticks)
plt.show()