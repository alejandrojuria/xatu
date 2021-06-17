# Script to represent the output from ./wavefunction program, i.e. the real-space 
# exciton wavefunction

import matplotlib.pyplot as plt
import numpy as np
import math

# ---------------- Read data from file -----------------------
file = open("rs_wf_soc_0_to_5_fixed2", "r")
holeLine = file.readline()

holeInfo = holeLine.split("\t")
holePosition = [float(holeInfo[0]), float(holeInfo[1])]

lines = file.readlines()
xArray = []
yArray = []
eDensity = []
fig, ax = plt.subplots(2, 5, sharex=True, sharey=True)
it = 0
previous_line = '0'

for line in lines:
    if line[0] == "#":
        ax[it%2, it//2].scatter(xArray, yArray, c=eDensity)
        ax[it%2, it//2].scatter(holePosition[0], holePosition[1], c="r", label="Hole")
        maxX = max(xArray)
        xArray = []
        yArray = []
        eDensity = []
        it += 1
        continue
    if previous_line[0] == "#":
        holeInfo = line.split("\t")
        holePosition = [float(holeInfo[0]), float(holeInfo[1])]
    else:
        lineData = line.split("\t")
        xArray.append(float(lineData[0]))
        yArray.append(float(lineData[1]))
        eDensity.append(float(lineData[2]))

    previous_line = line

#plt.colorbar()
plt.title("Real-space exciton w.f. (4 bands, w/o SOC)")
ax[0,0].set_xlabel("$x(A)$", fontsize=13)
ax[0,0].set_ylabel("$y(A)$", fontsize=13)
# plt.legend()
plt.show()