import matplotlib.pyplot as plt 
import numpy as np 


file = open('spectrum_bulk_spin_evolution_2', 'r')

lines = file.readlines()

nk = []
spinX = []
for line in lines:
    line = line.split('\t')
    nk.append(int(line[0]))
    spinX.append(float(line[1]))

file = open('spectrum_bulk_energy_evolution_2', 'r')
file.readline() # Skip first line

lines = file.readlines()
energy = []
for line in lines:
    line = line.split('\t')
    energy.append(float(line[1]))

plt.figure()
plt.plot(nk, spinX, 'g')
plt.title(r'$S^T_z$ convergence with nk')
plt.xlabel('Number k points')
plt.ylabel(r'$S^T_z (\hbar = 1)$')

plt.figure()
plt.plot(nk, energy, 'g')
plt.title(r'$S^T_z$ convergence with nk')
plt.xlabel('Number k points')
plt.ylabel(r'$S^T_z (\hbar = 1)$')


plt.show()