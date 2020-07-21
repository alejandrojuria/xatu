import numpy as np
import matplotlib.pyplot as plt

file = open("exciton_convergence_ribbon_width_2bands", "r")
file.readline()

lines = file.readlines()
E = []
dos = []
for line in lines:
    line = line.split('\t')
    E.append(float(line[0]))
    dos.append(float(line[1]))

plt.plot(E, dos, 'r-+')
plt.title("Bulk exciton energy w/ ribbon width (2 bands)")
plt.ylabel(r"$\epsilon (eV)$")
plt.xlabel(r"N")
plt.show()