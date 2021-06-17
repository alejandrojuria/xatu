import matplotlib.pyplot as plt 
import numpy as np 

file = open('band_derivative', 'r')
lines = file.readlines()
kpoints = []
derivative = []
derivative2 = []

for line in lines:
    line = line.split('\t')
    kpoints.append(float(line[0]))
    derivative.append(float(line[1]))
    derivative2.append(float(line[2]))


kpoints = np.array(kpoints)
derivative = np.array(derivative)
derivative2 = np.array(derivative2)

plt.plot(kpoints, derivative, 'g-')
plt.plot(kpoints, derivative, 'b+')
#plt.plot(kpoints, derivative2, 'g-')
#plt.plot(kpoints, derivative2, 'b+')
plt.legend(['v'])
plt.title("Derivative of edge valence band")
plt.xlabel("k")
plt.ylabel("$\epsilon'$")
plt.show()