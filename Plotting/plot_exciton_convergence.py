# Script to plot energy convergence of the exciton fundamental state
# written in the text file generated from the main.cpp code

import matplotlib.pyplot as plt
import numpy as np

file = open("twobands_Ncell_10_to_500", "r")
lines = file.readlines()
kpoints = []
energies = []
for n, line in enumerate(lines):
    if (n == 0):
        continue
    line = line.split('\t')
    kpoints.append(float(line[0]))
    allEnergy = [float(i) for i in line[1:]]
    energies.append(allEnergy)

file.close()
energies = np.array(energies)
print(kpoints)

plt.plot(kpoints, energies[:, 0], 'g+')
plt.plot(kpoints, energies[:, 1], 'r+')
plt.plot(kpoints, energies[:, 2], 'b+')

plt.plot(kpoints, energies[:, 0], 'g-')
plt.plot(kpoints, energies[:, 1], 'r-')
plt.plot(kpoints, energies[:, 2], 'b-')

plt.title('Exciton energy convergence')
plt.ylabel('$E$ (eV)', fontsize=15)
plt.xlabel('#k points', fontsize=15)
plt.axis('tight')
plt.grid(True, axis='both')
plt.legend(["Total energy", "TB energy", "V energy"])

plt.show()