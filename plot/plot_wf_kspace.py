# Script to plot energy convergence of the exciton fundamental state
# written in the text file generated from the main.cpp code

import matplotlib.pyplot as plt
import numpy as np
import math

file = open("exciton_bulk_wf", "r")
lines = file.readlines()
kpoints = []
coefs = []
for n, line in enumerate(lines):
    line = line.split('\t')
    kpoints.append(float(line[0]))
    coefs.append(float(line[1]))

file.close()
print(coefs)

width = 2*math.pi/((n+1)*4.5332)
plt.bar(kpoints, coefs, width, align='center', color="blue")

plt.title('First exciton wavefunction in k-space (Ncell = 500)')
plt.ylabel('$|A(k)|^2$', fontsize=15)
plt.xlabel('$k (A^{-1})$', fontsize=15)
plt.axis('tight')

plt.show()