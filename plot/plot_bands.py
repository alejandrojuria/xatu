# Script to plot bands from the text file generated from the 
# zigzag.cpp code

import matplotlib.pyplot as plt
import numpy as np

N = 15
dimTB = int(2*(N+1)*8)

plotBands = True
plotEdgesPosition = False
plotSpin = True

# --------------------- Plot bands ---------------------
if(plotBands):
    file = open("bands_N15_k500_nosoc.txt", "r")
    lines = file.readlines()
    kpoints = []
    energies = []
    for line in lines:
        line = line.split('\t')
        kpoints.append(float(line[0]))
        kEnergy = [float(i) for i in line[1:-1]]
        energies.append(kEnergy)

    energies = np.array(energies)
    kpoints = np.array(kpoints)

    iband = 2*(N+1)*5 - 6
    fband = 2*(N+1)*5 - 4
    plt.figure()
    for i in range(iband, fband):
        plt.plot(kpoints , energies[:, int(i)], 'g-')
    plt.plot(kpoints, energies[:, fband + 1], 'b-')
    plt.title('Valence band braiding')
    plt.xlabel(r'$k (A^{-1})$')
    plt.ylabel(r'$\epsilon (eV)$')
    
    plt.figure()
    gap = energies[:, int(2*(N+1)*5+1)] - energies[:, int(2*(N+1)*5-1)]
    #plt.plot(kpoints, gap, 'g-')
    prevGap = gap[0]
    for i in range(len(kpoints)):
        actualGap = gap[i]
        if(actualGap > prevGap):
            break
        prevGap = actualGap
    nk = i
    plt.plot(kpoints[nk-1:], gap[nk-1:], 'g-')
    rk = list(reversed(kpoints[:nk-1])) + kpoints[nk - 1]
    plt.plot(rk, gap[:nk-1], 'r-')
    plt.title("Gap")
    plt.legend(["Right branch", "Left branch"])
    plt.ylabel("Gap (eV)")
    plt.xlabel("$k (A^{-1})$")


# --------------------- Plot edge localization ---------------------
if(plotEdgesPosition):
    file2 = open("bands_position_edge")
    lines = file2.readlines()
    posv = []
    posc = []
    for line in lines:
        line = line.split("\t")
        posv.append(float(line[0]))
        posc.append(float(line[1]))

    plt.figure()
    plt.scatter(kpoints, energies[:, int(2*(N+1)*5-2)], c=posv)
    plt.scatter(np.array(kpoints), energies[:, int(2*(N+1)*5 + 1)], c=posc)
    plt.title("Edge state localization (left)")
    plt.xlabel("$k(A^{-1})$", fontsize = 13)
    plt.ylabel("$\epsilon (eV)$", fontsize = 13)
    plt.colorbar()

# --------------------- Plot band expected spin ---------------------
if(plotSpin):
    file_spin = open("bands_spin")
    lines = file_spin.readlines()
    spinV = []
    spinV3 = []
    spinV2 = []
    spinV4 = []
    spinxV = []
    spinxC = []
    for line in lines:
        line = line.split("\t")
        spinV.append(float(line[1]))
        spinV2.append(float(line[2]))
        spinV3.append(float(line[3]))
        spinV4.append(float(line[4]))
        spinxV.append(float(line[5]))
        spinxC.append(float(line[6]))



    fig, ax = plt.subplots(figsize=(10,5))
    axins = ax.inset_axes([0.05, 0.5, 0.205, 0.35])

    ax.plot(kpoints, spinV, 'r-')
    ax.plot(kpoints, spinV2, 'b-')
    ax.plot(kpoints, spinV3, 'g-')
    ax.plot(kpoints, spinV4, 'c-')
    axins.plot(kpoints, spinV, 'r-')
    axins.plot(kpoints, spinV2, 'b-')
    axins.plot(kpoints, spinV3, 'g-')
    axins.plot(kpoints, spinV4, 'c-')

    x1, x2, y1, y2 = 0.52, 0.54, 0.38, 0.4
    axins.set_xlim(x1, x2)
    axins.set_ylim(y1, y2)
    axins.set_xticklabels('')
    axins.set_yticklabels('')   

    ax.indicate_inset_zoom(axins) # Plots rectangle to show where inset is coming from

    ax.legend(['First valence band', 'First valence band deg.', 'Sec. valence band', 'Sec. valence band'])
    ax.title.set_text(r'$S_z$ for valence bands')
    ax.set_ylabel(r'$S_z (\hbar = 1)$')
    ax.set_xlabel(r'$k(A^{-1})$')

    '''
    combined_data = np.array([spinV,spinV3])
    _min, _max = np.amin(combined_data), np.amax(combined_data)
    plt.figure()
    fig, ax = plt.subplots(1,2)
    sca1 = ax[0].scatter(kpoints, energies[:, int(2*(N+1)*5 - 4)], c=spinV)
    sca1.set_clim(_min, _max)
    sca1 = ax[0].scatter(kpoints, energies[:, int(2*(N+1)*5 - 5)], c=spinV3)
    sca1 = ax[0].scatter(kpoints, energies[:, int(2*(N+1)*5 - 6)], c=spinV4)
    sca1 = ax[0].scatter(kpoints, energies[:, int(2*(N+1)*5) - 3], c=spinV2)
    sca1.set_clim(_min, _max)
    fig.colorbar(sca1, ax=ax[0])
    ax[1].plot(kpoints, spinV, 'r-')
    ax[1].plot(kpoints, spinV2, 'g-')
    ax[1].plot(kpoints, spinV3, 'c-')
    ax[1].plot(kpoints, spinV4, 'b-')
    ax[0].set_ylabel("$\epsilon (eV)$", fontsize = 13)
    ax[0].set_xlabel("$k(A^{-1})$", fontsize = 13)
    ax[1].set_ylabel("$<S_z> (\hbar = 1)$", fontsize = 13)
    ax[1].set_xlabel("$k(A^{-1})$", fontsize = 13)
    ax[1].legend(["Valence", "Conduction"])
    fig.suptitle("Expected spin value $<S_z>$")
'''

'''     plt.figure()
    plt.plot(kpoints, spinV, 'r-')
    plt.plot(kpoints, spinC, 'c-')
    plt.title("Expected spin value $S_z$")
    plt.xlabel("$k(A^{-1})$", fontsize = 13)
    plt.ylabel("$<S_z> (\hbar = 1)$", fontsize = 13)
    plt.legend(["Valence", "Conduction"])

    plt.figure()
    gap = energies[:, int(2*(N+1)*5 + 1)] - energies[:, int(2*(N+1)*5 - 2)]
    gap2 = energies[:, int(2*(N+1)*5 )] - energies[:, int(2*(N+1)*5 - 1)]

    plt.plot(kpoints, gap)
    plt.plot(kpoints, gap2, "g-")
    plt.title("Gap between edge bands")
    plt.xlabel("$k(A^{-1})$")
    plt.ylabel("$\Delta\epsilon (eV)$")
    plt.legend(["Left", "Right"])

    #plt.plot(kpoints, posv)
    #plt.plot(kpoints, posc)

    fig2, ax2 = plt.subplots(1,3)
    ax2[0].plot(kpoints, spinxV, "g-")
    ax2[0].plot(kpoints, spinxC, "b-")
    ax2[0].set_xlabel("$k(A^{-1})$")
    ax2[0].set_ylabel("$<S_x>$")
    ax2[0].legend(["Valence", "Conduction"])

    ax2[1].plot(kpoints, spinyV, "g-")
    ax2[1].plot(kpoints, spinyC, "b-")
    ax2[1].set_xlabel("$k(A^{-1})$")
    ax2[1].set_ylabel("$<S_y>$")
    ax2[1].legend(["Valence", "Conduction"])

    ax2[2].plot(kpoints, spinV, "g-")
    ax2[2].plot(kpoints, spinC, "b-")
    ax2[2].set_xlabel("$k(A^{-1})$")
    ax2[2].set_ylabel("$<S_z>$")
    ax2[2].legend(["Valence", "Conduction"])'''


plt.show()
