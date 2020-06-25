# Script to plot bands from the text file generated from the 
# zigzag.cpp code

import matplotlib.pyplot as plt
import numpy as np

N = 15
dimTB = int(2*(N+1)*8)
iband = 2*(N+1)*5 - 8
fband = 2*(N+1)*5 + 8

# --------------------- Plot bands ---------------------
file = open("bands_N15_k500_nosoc.txt", "r")
lines = file.readlines()
kpoints = []
energies = []
for line in lines:
    line = line.split('\t')
    kpoints.append(float(line[0]))
    kEnergy = [float(i) for i in line[1:-1]]
    energies.append(kEnergy)

energies = np.array(energies)
kpoints = np.array(kpoints)

plt.figure()
for i in range(iband, fband):
    plt.plot(kpoints , energies[:, int(i)], 'g-')
plt.title("Energy bands")
plt.xlabel(r"$k(A^{-1})$")
plt.ylabel(r"$\epsilon (eV)$")

if(True):
    file_spin = open("spin_bands")
    lines = file_spin.readlines()
    spinV = []
    spinC = []
    spinC_2 = []
    spinC_3 = []
    spinxV = []
    spinxC = []
    for line in lines:
        line = line.split("\t")
        spinV.append(float(line[1]))
        spinC.append(float(line[2]))
        spinC_2.append(float(line[3]))
        spinC_3.append(float(line[4]))
        spinxV.append(float(line[5]))
        spinxC.append(float(line[6]))

    combined_data = np.array([spinV,spinC])
    _min, _max = np.amin(combined_data), np.amax(combined_data)
    plt.figure()
    fig, ax = plt.subplots(1,2)
    ax[1].plot(kpoints, spinC_2, 'r-')
    ax[1].plot(kpoints, spinC, 'c-')
    ax[0].set_ylabel("$\epsilon (eV)$", fontsize = 13)
    ax[0].set_xlabel("$k(A^{-1})$", fontsize = 13)
    ax[1].set_ylabel("$<S_z> (\hbar = 1)$", fontsize = 13)
    ax[1].set_xlabel("$k(A^{-1})$", fontsize = 13)
    ax[1].legend(["Valence", "Conduction"])
    fig.suptitle("Expected spin value $<S_z>$")



plt.show()