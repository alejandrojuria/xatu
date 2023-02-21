import matplotlib.pyplot as plt
import numpy as np

file = open("kwf.out", "r")
lines = file.readlines()
kpoints = []
ks = []
coefs = []
states = []
for line in lines:
    if line[0] == "k":
        continue
    if line[0] == "#":
        states.append(np.array(coefs))
        ks.append(kpoints)
        kpoints = []
        coefs = []
        continue
    line = line.split()
    kpoint = [float(num) for num in line[0:3]]
    kpoints.append(kpoint)
    coefs.append(float(line[-1]))

kpoints = np.array(ks[0])
file.close()

fig, ax = plt.subplots(1, 5, sharey=True, figsize=(15, 3), dpi=100)

degeneracies = [2, 1, 2, 1, 2]
wfs = []
counter = 0
for deg in degeneracies:
    state = np.zeros(states[0].shape)
    for i in range(deg):
        state += states[counter + i]
    counter += deg
    wfs.append(state)

for i, axis in enumerate(ax):
    axis.scatter(kpoints[:, 0], kpoints[:, 1], c=wfs[i], s=2)



plt.show()