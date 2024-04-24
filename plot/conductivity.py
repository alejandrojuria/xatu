import matplotlib.pyplot as plt
import numpy as np

# -----------------------------------------------------------
# Configure plot
# -----------------------------------------------------------

exciton_file = "hBN_ex.dat"
sp_file      = "hBN_sp.dat"
fontsize     = 15
markersize   = 110  # Adjust size of points in plot
latex        = True # Set to True if latex is installed
if latex:
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif"
    })
plt.rcParams["axes.labelsize"] = fontsize
plt.rcParams["axes.titlesize"] = fontsize
plt.rcParams["xtick.labelsize"] = fontsize
plt.rcParams["ytick.labelsize"] = fontsize
plt.rcParams["legend.fontsize"] = fontsize
plt.rcParams["font.size"] = fontsize

# -----------------------------------------------------------
# First extract conductivies from files
# -----------------------------------------------------------

file = open(exciton_file, "r")
energy, sigma_xx, sigma_yy, sigma_sp_xx, sigma_sp_yy = [], [], [], [], []
for line in file.readlines():
    line = [float(value) for value in line.split()]
    if not line: # Skip reading exciton velocity matrix elements
        break
    energy.append(line[0])
    sigma_xx.append(line[1])
    sigma_yy.append(line[5])
    
file = open(sp_file, "r")
for line in file.readlines():
    line = [float(value) for value in line.split()]
    sigma_sp_xx.append(line[1])
    sigma_sp_yy.append(line[5])

# -----------------------------------------------------------
# Plot states
# -----------------------------------------------------------

plt.plot(energy, sigma_xx, "g-", label="xx")
plt.plot(energy, sigma_yy, "r-", label="yy")
plt.plot(energy, sigma_sp_xx, "b--", label="sp_xx")
plt.plot(energy, sigma_sp_yy, "c--", label="sp_yy")
plt.xlabel("E (eV)")
plt.ylabel(r"$\sigma$ ($e^2/\hbar$)")

plt.legend(fontsize=fontsize/1.5)

plt.show()