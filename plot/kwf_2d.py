import matplotlib.pyplot as plt
import numpy as np

file = open("kwf.out", "r")
lines = file.readlines()
kpoints = []
coefs = []
states = []
for line in lines:
    if line[0] == "#":
        states.append(coefs)
        coefs = []
        continue
    elif line[0] == "k":
        kpoints = []
        continue
    line = line.split()
    kpoint = [float(num) for num in line[0:3]]
    kpoints.append(kpoint)
    coefs.append(float(line[-1]))

kpoints = np.array(kpoints)
file.close()

file = open("crystal.out", "r")
lines = file.readlines()
print(lines[0].split())
K = np.array([float(val) for val in lines[0].split()])
rotatedK = np.array([float(val) for val in lines[1].split()])
rotatedK2 = np.array([float(val) for val in lines[2].split()])
M = np.array([float(val) for val in lines[3].split()])
file.close()

xHex = [K[0], -rotatedK2[0], rotatedK[0], -K[0], rotatedK2[0], -rotatedK[0], K[0]]
yHex = [K[1], -rotatedK2[1], rotatedK[1], -K[1], rotatedK2[1], -rotatedK[1], K[1]]

xIBZ = [0, -rotatedK[0], M[0], 0]
yIBZ = [0, -rotatedK[1], M[1], 0]

fig, ax = plt.subplots(1, 2, sharey=True)

state = np.array(states[0]) + np.array(states[1])
ax[0].scatter(kpoints[:, 0], kpoints[:, 1], c=state, marker='H')
ax[0].plot(xHex, yHex, 'r-')
ax[0].plot(xIBZ, yIBZ, 'r-')
ax[0].text(-rotatedK[0], -rotatedK[1], "K", fontsize=20, color='red')
ax[0].text(M[0], M[1], "M", fontsize=20, color='red')
ax[0].text(-0.3, -0.1, rf"$\Gamma$", fontsize=20, color='red')
ax[0].axis('equal')




plt.show()