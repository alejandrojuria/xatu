# Script to plot the edge exciton spectrum

import matplotlib.pyplot as plt
import numpy as np

N = 15

file = open("exciton_spectrum_exact", "r")
file.readline() # Skip header

lines = file.readlines()
Qpoints = []
energies = []
refEnergy = 0
for n, line in enumerate(lines):
    line = line.split('\t')
    Qpoints.append(float(line[0]))
    QEnergy = [float(i) for i in line[1:-1]]
    if(n == 0):
        refEnergy = QEnergy[0]
    QEnergy = np.array(QEnergy)
    energies.append(QEnergy)

Qpoints = np.array(Qpoints)
energies = np.array(energies) - refEnergy

plt.plot(Qpoints, energies[:, 0], 'b-')
plt.plot(-Qpoints, energies[:, 0], 'b-')

for i in range(len(QEnergy)-1):
    plt.plot(Qpoints, energies[:, i + 1], 'g+')
    plt.plot(-Qpoints, energies[:, i + 1], 'g+')

adjustAxis = False
if adjustAxis:
    plt.ylim(top=0.3, bottom=0.0)
    plt.xlim(right=0.25, left=-0.25)

plt.title('Exciton spectrum (both edges, nk = 800)')
plt.xlabel('$Q (A^{-1})$', fontsize = 13)
plt.ylabel('$\epsilon (eV)$', fontsize = 13)

#for i in [2*(N+1)*5 - 0, 2*(N+1)*5 + 2]:
#    plt.plot(kpoints, energies[:, i], color='green')

plt.show()
