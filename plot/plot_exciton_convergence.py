# Script to plot energy convergence of the exciton fundamental state
# written in the text file generated from the main.cpp code

import matplotlib.pyplot as plt
import numpy as np

gap = 3.625*2

file = open("./hBN/e_conv_approx_1_20", "r")
file2 = open("./hBN/e_conv_exact_1_15", "r")
lines = file.readlines()

ncell = []
energies = []
for line in lines[1:]:
    if line[0] == "#":
        continue
    line = line.split('\t')
    ncell.append(float(line[0]))
    allEnergy = [float(i) for i in line[1:-2]]
    energies.append(allEnergy)

file.close()

firststate = []
for energy in energies:
    firststate.append(energy[0] - gap)

lines = file2.readlines()
ncell2 = []
energies2 = []
for line in lines[1:]:
    if line[0] == "#":
        continue
    line = line.split('\t')
    ncell2.append(float(line[0]))
    allEnergy = [float(i) for i in line[1:-2]]
    energies2.append(allEnergy)

exactstate = []
for energy in energies2:
    exactstate.append(energy[0] - gap)

exciton_expected_e = -np.ones(26)*(1.932)

plt.plot(ncell, firststate, 'g+')
plt.plot(ncell, firststate, 'g-')

plt.plot(ncell2, exactstate, 'b+')
plt.plot(ncell2, exactstate, 'b-')

ncell = range(1, 27)
plt.plot(ncell, exciton_expected_e, "r-")

plt.title('hBN exciton energy convergence')
plt.ylabel(r'$E_X$ (eV)', fontsize=15)
plt.xlabel('#N. cell', fontsize=15)
plt.ylim((-2, 0))
plt.axis('tight')
plt.grid(True, axis='both')
plt.legend(["_", "Approximate", "_", "Exact", "Reference energy"])

plt.show()