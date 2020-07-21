import numpy as np
import matplotlib.pyplot as plt

file = open("ribbon_width_maxE", "r")

lines = file.readlines()
N = []
maxE = []
for line in lines:
    line = line.split('\t')
    N.append(float(line[0]))
    maxE.append(float(line[1]))

plt.plot(N, maxE, 'r-+')
plt.title("Max energy available to non-int. e-h edge pairs")
plt.ylabel(r"$\epsilon(eV)$")
plt.xlabel("N (ribbon width)")
plt.show()