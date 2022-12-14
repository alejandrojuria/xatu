import matplotlib.pyplot as plt
import numpy as np

# ---------------- Read data from file -----------------------
file = open("./rswf.out", "r")

lines = file.readlines()
xArray = []
yArray = []
eDensity = []
states = []

for line in lines:
    if line[0] == "#":
        states.append(eDensity[1:])            
        xArray, yArray, eDensity = [], [], []
        continue
    
    lineData = line.split()
    xArray.append(float(lineData[0]))
    yArray.append(float(lineData[1]))
    eDensity.append(float(lineData[2]))
states.append(eDensity[1:]) # Append last state

# Remove first line which corresponds to hole position
holePosition = [xArray[0], yArray[0]]
xArray = np.array(xArray[1:])
yArray = np.array(yArray[1:])

# -----------------------------------------------------------

degeneracies = [1, 1, 1, 1, 1, 1]
wfs = []
counter = 0
for degeneracy in degeneracies:
    wf = np.array(states[counter])
    for stateIdx in range(counter + 1, degeneracy + counter):
        wf += np.array(states[stateIdx])
    wfs.append(wf)
    counter += degeneracy

fig, ax = plt.subplots(2, 3, sharex=True, sharey=True, dpi=200, figsize=(6,4))

for n, state in enumerate(wfs):
    state /= np.max(state)
    j = int(n%3)
    i = int(n/3)

    ax[i, j].scatter(xArray - holePosition[0], yArray - holePosition[1], c=state, cmap="Greens", s=30 - 30*(state - 1)**10)
    ax[i, j].scatter(0, 0, c="r", label="Hole")

for i, axis_array in enumerate(ax):
    for axis in axis_array:
        axis.set_xlim([-20, 20])
        axis.set_ylim([-20, 20])
        axis.set_aspect("equal", "box")
        axis.tick_params(axis='both', which='both', labelsize=13)
        for site in ['top', 'bottom', 'left', 'right']:
            axis.spines[site].set_linewidth(2)
        if i == 1:
            axis.set_xlabel(r"$x(\AA)$", fontsize=13)

ax[0, 2].set_xlabel(r"$x(\AA)$", fontsize=13)
ax[0, 2].set_xticks([-10, 0, 10])
ax[0, 0].set_ylabel(r"$y(\AA)$", fontsize=13)
ax[1, 0].set_ylabel(r"$y(\AA)$", fontsize=13)

ax[0, 0].legend()
plt.subplots_adjust(wspace=0.1, hspace=0.1)
plt.savefig("hbn_rswf.png", bbox_inches='tight')
plt.show()