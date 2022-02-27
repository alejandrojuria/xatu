import matplotlib.pyplot as plt 
import numpy as np 

file = open("eigval_convergence_reduced_mesh_f2.out", "r")

for _ in range(4):
    file.readline()

lines = file.readlines()

screening = []
E, T, V = [], [], []
for line in lines:
    line = line.split("\t")
    screening.append(float(line[0]))
    E.append(float(line[1]))
    T.append(float(line[2]))
    V.append(float(line[3]))

plt.plot(screening, E, 'g-')
plt.plot(screening, T, 'r-')
plt.plot(screening, V, 'b-')
plt.legend(["Total energy", "Band energy", "Potential energy"])

plt.xlabel(r"$r_0 (\AA)$, ")
plt.ylabel(r'$E_X (eV)$')
plt.title("Exciton energy convergence with screening length")
plt.show()
