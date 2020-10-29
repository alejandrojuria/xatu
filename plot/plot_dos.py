import matplotlib.pyplot as plt 
import numpy as np

file = open("dos_eh_pair", "r")
lines = file.readlines()
energies = []
dos = []
for line in lines:
    line = line.split("\t")
    energies.append(float(line[0]))
    dos.append(float(line[1]))

energies = np.array(energies)
dos = np.array(dos)

plt.plot(energies, dos, "b-")
plt.title("Edge DOS (nk = 1001)")
plt.xlabel("$E (eV)$")
plt.ylabel(r'$\rho(E)$')
plt.show()