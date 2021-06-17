import matplotlib.pyplot as plt
import numpy as np

file = open("bi_entanglement_spectrum", "r")
lines = file.readlines()
spectrum = []
for line in lines:
    line = line.split()
    band = [float(num) for num in line]
    spectrum.append(band)

spectrum = np.array(spectrum)

for band in spectrum.T:
    plt.plot(band, "g+")

plt.show()