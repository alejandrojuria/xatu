# Script to plot spin eigenvalues obtained from bands program (libzigzag
# implementation)

import matplotlib.pyplot as plt
import numpy as np
import math

file = open("spin_eigenvalues", "r")
lines = file.readlines()
kpoints = []
eigenvalues = []
for line in lines:
    line = line.split('\t')
    kpoints.append(float(line[0]))
    auxEig = [float(i) for i in line[1:-1]]
    eigenvalues += auxEig

file.close()

itp = 0
itm = 0
for val in eigenvalues:
    if(val > 0):
        itp += 1
    elif(val < 0):
        itm -= 1

print(itp)
print(itm)


# Plotting
plt.plot(np.sort(eigenvalues), 'ro')
plt.title('Ground-state projected spin ($P_z$) eigenvalues')
plt.ylabel('$P_z$ eigenvalue ($\hbar$ units)', fontsize=12)
plt.xlabel('State index', fontsize=12)
plt.axis('tight')

plt.show()