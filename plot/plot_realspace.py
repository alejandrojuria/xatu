# Script to represent the output from ./wavefunction program, i.e. the real-space 
# exciton wavefunction

import matplotlib.pyplot as plt
import numpy as np
import math

# ---------------- Read data from file -----------------------
#file = open("rs_wf_first_301nk_4bands_noedge", "r")
file = open("bands_TB_wavefunction", "r")
#holeLine = file.readline()

#holeInfo = holeLine.split("\t")
#holePosition = [float(holeInfo[0]), float(holeInfo[1])]

lines = file.readlines()
xArray = []
yArray = []
eDensity = []

for line in lines:
    lineData = line.split("\t")
    xArray.append(float(lineData[0]))
    yArray.append(float(lineData[1]))
    eDensity.append(float(lineData[2]))

# -----------------------------------------------------------
# ---------------------- Format data ------------------------
#x0 = xArray[0]
#npts = 0
#for x in xArray:
#    if (x == x0):
#        npts += 1
#    else:
#        break

#yArray = yArray[0:npts]
#auxXarray = [x0]
#for x in xArray:
#    if (x != x0):
#        x0 = x
#        auxXarray.append(x0)
#xArray = auxXarray

# eDensity = np.transpose(np.array(eDensity).reshape(npts, npts))

# -----------------------------------------------------------
# ------------------------ Plotting -------------------------

#h = plt.pcolormesh(xArray, yArray, eDensity)
#plt.colorbar(h)
#plt.title("Real-space exciton wavefunction (2 bands, no edges)")


h = plt.scatter(xArray, yArray, c=eDensity)
maxX = max(xArray)
#plt.scatter(holePosition[0], holePosition[1], c="r", label="Hole")
plt.colorbar(h)
plt.title("Real-space exciton w.f. (4 bands, w/ edges)")
#plt.ylim((-maxX/2, maxX/2))
plt.xlabel("$x(A)$", fontsize=13)
plt.ylabel("$y(A)$", fontsize=13)
#plt.legend()
plt.show()