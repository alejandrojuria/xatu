import matplotlib.pyplot as plt 
import numpy as np 


file = open('exciton_bulk_spin_2b_z1em8', 'r')

lines = file.readlines()

nk = []
spinX = []
for line in lines:
    line = line.split('\t')
    nk.append(int(line[0]))
    spinX.append(float(line[1]))

file = open('exciton_bulk_energy_2b_z1em8', 'r')
file.readline() # Skip first line

lines = file.readlines()
energy = []
for line in lines:
    line = line.split('\t')
    energy.append(float(line[1]))

fig, ax = plt.subplots(1, 2)
ax[0].plot(nk, spinX, 'g')
ax[1].plot(nk, energy, 'g')

file = open('exciton_bulk_spin_4b_z1em8', 'r')

lines = file.readlines()

nk = []
spinX = []
for line in lines:
    line = line.split('\t')
    nk.append(int(line[0]))
    spinX.append(float(line[1]))

file = open('exciton_bulk_energy_4b_z1em8', 'r')
file.readline() # Skip first line

lines = file.readlines()
energy = []
for line in lines:
    line = line.split('\t')
    energy.append(float(line[1]))

ax[0].plot(nk, spinX, 'b')
ax[1].plot(nk, energy, 'b')

file = open('exciton_bulk_spin_6b_z1em8', 'r')

lines = file.readlines()

nk = []
spinX = []
for line in lines:
    line = line.split('\t')
    nk.append(int(line[0]))
    spinX.append(float(line[1]))

file = open('exciton_bulk_energy_6b_z1em8', 'r')
file.readline() # Skip first line

lines = file.readlines()
energy = []
for line in lines:
    line = line.split('\t')
    energy.append(float(line[1]))

ax[0].plot(nk, spinX, 'r')
ax[0].title.set_text(r'$S^T_z$ convergence with nk')
ax[0].set_xlabel('Number of k points')
ax[0].set_ylabel(r'$S^T_z (\hbar = 1)$')

ax[1].plot(nk, energy, 'r')
ax[1].title.set_text(r'$\epsilon$ convergence with nk')
ax[1].set_xlabel('Number of k points')
ax[1].set_ylabel(r'$\epsilon (eV)$')
plt.legend(["2 bands", "4 bands", "6 bands"])


plt.show()