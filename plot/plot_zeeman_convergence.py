# Script to plot wavefunction convergence with zeeman term of the exciton fundamental state

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import math

N = 15

# file = open("edge_wf_Qneq0_conv", "r")
file = open("spectrum_bulk_wf_2bands", "r")
lines = file.readlines()
kpoints = []
coefs = []
it = 0
colors = cm.viridis(np.linspace(0, 1, 25))
fig, ax = plt.subplots(figsize=(10,5))
#axins = ax.inset_axes([0.05, 0.5, 0.205, 0.35])


for n, line in enumerate(lines):
    line = line.split('\t')
    if(line[0] == '#\n'):
        ax.plot(kpoints, coefs, color = colors[it])
        #axins.plot(kpoints, coefs, color = colors[it])
        kpoints = []
        coefs = []
        it += 1
    else:
        kpoints.append(float(line[0]))
        coefs.append(float(line[1]))


# sub region of the original image
#x1, x2, y1, y2 = 0.67, 0.69, -0.1, 1.8
#axins.set_xlim(x1, x2)
#axins.set_ylim(y1, y2)
#axins.set_xticklabels('')
#axins.set_yticklabels('')

#ax.indicate_inset_zoom(axins) # Plots rectangle to show where inset is coming from

plt.legend(["400", "1000", "2000", "4000", "6000", "10000"], title="#k", ncol = 2)
plt.title("Convergence of fundamental exciton w.f. ($Q>0$, left edge, $L^2$ norm)")
plt.ylabel("$|A(k)|^2$", fontsize = 13)
plt.xlabel("$k (A^{-1})$", fontsize = 13)
plt.show()
file.close()