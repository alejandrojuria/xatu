import matplotlib.pyplot as plt 
import numpy as np 

file = open('overlap', 'r')
lines = file.readlines()
kpoints = []
overlapv = []
overlapc = []
overlapp = []
overlapp2 = []

for line in lines:
    line = line.split('\t')
    kpoints.append(float(line[0]))
    overlapv.append(float(line[1]))
    overlapc.append(float(line[2]))
    overlapp.append(float(line[3]))
    overlapp2.append(float(line[4]))


kpoints = np.array(kpoints)
overlapv = np.array(overlapv)
overlapc = np.array(overlapc)
overlapp = np.array(overlapp)
overlapp2 = np.array(overlapp2)


plt.plot(kpoints, overlapv, 'g-')
plt.plot(kpoints, overlapc, 'b-')
plt.plot(kpoints, overlapp, 'r-')
plt.plot(kpoints, overlapp2, 'c-')

plt.legend(['v', 'c'])
plt.title("Overlap between |k=0> and |k'> (edge states)")
plt.xlabel("k'")
plt.ylabel("Overlap")
plt.show()