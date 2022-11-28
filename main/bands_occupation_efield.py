# Script to plot bands 


import matplotlib.pyplot as plt
import numpy as np

plt.rc('font', family='serif')


N = 19
dimTB = int(2*(N+1)*8)

file = open("eigval.out", "r")
lines = file.readlines()
energies = []
left, right = [], []
for line in lines:
    line = line.split()
    line = [float(value) for value in line]
    energy = line[0::3]
    left_occ = line[1::3]
    right_occ = line[2::3]
    energies.append(energy)
    left.append(left_occ)
    right.append(right_occ)

energies = np.array(energies, dtype=object)
left = np.array(left, dtype=object)
right = np.array(right, dtype=object)
occupation = left - right

iband = 2*(N+1)*2
fband = 2*(N+1)*8
nk = len(energies[:, 0])
k = range(nk)
fig, ax = plt.subplots(1, 1, dpi=200, figsize=(5.5, 4))
fermiEnergy = np.sort(energies.reshape(-1, ))[nk*(2*(N+1)*5) - 1]
energies -= fermiEnergy
for i in range(iband, fband):
    if i in [2*(N+1)*5 - 2, 2*(N+1)*5 - 1, 2*(N+1)*5, 2*(N+1)*5 + 1]:
        fig = ax.scatter(k, energies[:, int(i)], c=occupation[:, int(i)], cmap="seismic", vmin=-1, vmax=1)
    else:
        ax.scatter(k, energies[:, int(i)], c="g")
        continue
cbar = plt.colorbar(fig, ax=ax, extend="both")
cbar.ax.tick_params(labelsize=14)
ax.set_ylim([-0.4, 0.4])

ax.set_xticks([0, nk/2, nk-1])
ax.set_xticklabels([r"K", r"$\Gamma$", "K"], fontsize=14)
ax.set_ylabel(r"$\varepsilon$ (eV)", fontsize=14)

ax.tick_params(axis='both', which='both', labelsize=14)
for site in ['top', 'bottom', 'left', 'right']:
    ax.spines[site].set_linewidth(3)


plt.savefig("ribbon_bands_N20_efield.png", bbox_inches="tight")
plt.show()

