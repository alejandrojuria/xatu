# Script to plot energy convergence of the exciton fundamental state
# written in the text file generated from the main.cpp code

import matplotlib.pyplot as plt
import numpy as np

file = open("exciton_convergence", "r")
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

file = open("exciton_convergence_10_to_400_V0=0", "r")
lines = file.readlines()
kpoints1 = []
energies1 = []
for n, line in enumerate(lines):
    if (n == 0):
        continue
    line = line.split('\t')
    kpoints1.append(float(line[0]))
    allEnergy = [float(i) for i in line[1:]]
    energies1.append(allEnergy)

file.close()
energies1 = np.array(energies1)

file = open("exciton_convergence_10_to_1000_Vcontour", "r")
lines = file.readlines()
kpoints2 = []
energies2 = []
for n, line in enumerate(lines):
    if (n == 0):
        continue
    line = line.split('\t')
    kpoints2.append(float(line[0]))
    allEnergy = [float(i) for i in line[1:]]
    energies2.append(allEnergy)

file.close()
energies2 = np.array(energies2)

plt.plot(kpoints, energies[:, 0], 'g+')
plt.plot(kpoints1, energies1[:, 0], 'r+')
plt.plot(kpoints2, energies2[:, 0], 'b+')

plt.plot(kpoints, energies[:, 0], 'g-')
plt.plot(kpoints1, energies1[:, 0], 'r-')
plt.plot(kpoints2, energies2[:, 0], 'b-')

plt.title('First exciton energy (all methods)')
plt.ylabel('$E$ (eV)', fontsize=15)
plt.xlabel('#k points', fontsize=15)
plt.axis('tight')
plt.grid(True, axis='both')
plt.legend(["Full discrete transform", "V(k=0)=0", "V(k=0) contour"])

plt.show()