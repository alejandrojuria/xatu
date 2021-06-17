# Script to represent the output from ./wavefunction program, i.e. the real-space 
# exciton wavefunction

import matplotlib.pyplot as plt
import numpy as np
import math

# ---------------- Read data from file -----------------------
#file = open("./rswf_bi_w8", "r")
#file = open("./bloch_state_K", "r")
#file = open("./exciton_bulk_wf_N15_nosoc", "r")
file = open("./hbn_rswf_approx_N30_hMid", "r")
holeLine = file.readline()

holeInfo = holeLine.split("\t")
holePosition = [float(holeInfo[0]), float(holeInfo[1])]

lines = file.readlines()
xArray = []
yArray = []
eDensity = []

for line in lines:
    if line[0] == "#":
        continue
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

#plt.hexbin(xArray, yArray, gridsize=(20,20), C=eDensity)


h = plt.scatter(xArray, yArray, c=eDensity)
maxX = max(xArray)
plt.scatter(holePosition[0], holePosition[1], c="r", label="Hole")
plt.colorbar(h)
plt.title("Real-space exciton w.f.")
plt.ylim((-maxX/2 + holePosition[1], maxX/2 + holePosition[1]))
plt.xlabel("$x(A)$", fontsize=13)
plt.ylabel("$y(A)$", fontsize=13)
plt.legend()
plt.show()