import matplotlib.pyplot as plt 
import numpy as np 

file = open("exciton_convergence", "r")
lines = file.readlines()

nk = []
E = []
for line in lines[1:]:
    line = line.split("\t")
    nk.append(int(line[0]))
    E.append(float(line[1]))

plt.plot(nk, E, 'g-')
plt.xlabel("nk")
plt.ylabel(r'$E_X (eV)$')
plt.title("Exciton energy convergence with nk")
plt.show()
