# Script to plot energy convergence of the exciton fundamental state
# in conjunction with edge bands

import matplotlib.pyplot as plt
import numpy as np
import math

N = 15

file = open("class_bulk_wf", "r")
lines = file.readlines()
kpoints = []
coefs = []
for n, line in enumerate(lines):
    line = line.split('\t')
    kpoints.append(float(line[0]))
    coefs.append(float(line[1]))

print(len(kpoints))
print(len(coefs))


file.close()

fig, ax1 = plt.subplots()

width = 2*math.pi/((n+1)*4.5332)
ax1.bar(kpoints, coefs, width, align='center', color="blue")

#coefsEdge = np.zeros(499)
#coefsEdge[86] = 5.
#coefsEdge[499 - 86] = 5.

#width = 4*2*math.pi/((n+1)*4.5332)
#ax1.bar(kpoints, coefsEdge, width, align='center', color="red")
#ax1.plot(kpoints, coefs, 'b-')
#ax1.plot(kpoints, coefs, 'b+')

plt.title('Bulk exciton & edge e-h pair')
ax1.set_ylabel('$|A(k)|^2$', fontsize=15)
ax1.set_xlabel('$k (A^{-1})$', fontsize=15)

ax1.legend(["Bulk exciton", "Edge e-h pair"])

#filebands = open("bands_N15_k500.txt", "r")
#lines = filebands.readlines()
#kpoints = []
#energies = []
"""
for line in lines:
    line = line.split('\t')
    kpoints.append(float(line[0]))
    kEnergy = np.array([float(i) for i in line[1:-1]])
    energies.append(kEnergy)

energies = np.array(energies)
filebands.close()

ax2 = ax1.twinx()

kpoints = np.array(kpoints)
for i in [2*(N+1)*5 - 2, 2*(N+1)*5 + 1]:
    ax2.plot(kpoints, energies[:, int(i)], color='green')
ax2.set_ylabel(r"$\epsilon (eV)$")
#Q = 0.0366
#ax2.plot(kpoints - Q, energies[:, 2*(N+1)*5], color='blue')
#ax2.plot(kpoints - Q, energies[:, 2*(N+1)*5 + 1], color='blue')
"""

fig.tight_layout()

plt.show()