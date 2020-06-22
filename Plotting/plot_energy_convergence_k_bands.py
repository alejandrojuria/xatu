# Script to plot energy convergence of the exciton fundamental state
# written in the text file generated from the main.cpp code

import matplotlib.pyplot as plt
import numpy as np

file = open("energy_convergence_k_nbands", "r")
lines = file.readlines()
bands = []
kpoints = []
energies = []

for n, line in enumerate(lines):
    if (n == 0):
        continue
    line = line.split('\t')
    bands.append(float(line[0]))
    kpoints.append(float(line[1]))
    allEnergy = [float(i) for i in line[2:]]
    energies.append(allEnergy)


file.close()
energies = np.array(energies)

band0 = bands[0]
kpointsBand = []
energiesBand = []
colors = ['r-', 'g-', 'b-']
it = 0
for i in range(len(bands)):
    if (bands[i] != band0):
        band0 = bands[i]
        plt.plot(kpointsBand, energiesBand, colors[it])
        kpointsBand = []
        energiesBand = []
        it += 1
    kpointsBand.append(kpoints[i])
    energiesBand.append(energies[i,0])
plt.plot(kpointsBand, energiesBand, colors[it])

plt.title('Exciton energy convergence (w/ num. bands, num. k)')
plt.ylabel('$E$ (eV)', fontsize=15)
plt.xlabel('#k points', fontsize=15)
plt.axis('tight')
plt.grid(True, axis='both')
plt.legend(["2 bands", "4 bands", "6 bands"])

plt.show()