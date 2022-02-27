# Script to plot bands from the text file generated from the 
# system.cpp code

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def rotate_C3(position):
    theta = 2*np.pi/3
    C3_rotation = np.array([[np.cos(theta), -np.sin(theta), 0],
                            [np.sin(theta),  np.cos(theta), 0],
                            [            0,              0, 1]])
    rotated_position = np.dot(C3_rotation, position)

    return rotated_position

# --------------------- Plot bands ---------------------
file = open("./k_wf.out", "r")
file.readline() # Skip header
lines = file.readlines()
kpoints = []
coefs = []
for line in lines:
    if line[0] == "#":
        continue
    line = line.split('\t')
    kpoint = [float(num) for num in line[0:3]]
    kpoints.append(kpoint)
    coefs.append(float(line[-1]))

coefs = np.array(coefs)
kpoints = np.array(kpoints)
for coef in coefs:
    if coef > 0.4:
        print(coef)

#rotated_kpoints = np.zeros(kpoints.shape)
#for i, kpoint in enumerate(kpoints):
#    rotated_kpoints[i, :] = rotate_C3(kpoint)

nk = int(np.sqrt(len(kpoints[:, 0])))
#x = kpoints[:, 0].reshape(nk, nk)
#y = kpoints[:, 1].reshape(nk, nk)
#z = coefs.reshape(nk, nk)

plt.figure()
#print(kpoints)
#plt.hexbin(kpoints[:, 0], kpoints[:, 1], gridsize=(50, 50), C=coefs, cmap="Reds")
plt.scatter(kpoints[:, 0], kpoints[:, 1], c=coefs)
plt.colorbar()
#plt.scatter(rotated_kpoints[:, 0], rotated_kpoints[:, 1], c=coefs)
#plt.scatter(rotated_kpoints[:, 0], rotated_kpoints[:, 1], c=coefs, cmap="Reds")


#fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
#ax.plot_trisurf(kpoints[:, 0], kpoints[:, 1], coefs)

plt.show()