import matplotlib.pyplot as plt 
import numpy as np 

file_one = open("transition_ribbon_2bands_100cell_same", "r")
file_two = open("transition_ribbon_2bands_200cell_same", "r")


lines = file_one.readlines()
trate_one = []
kpoints1 = []
for n, line in enumerate(lines):
    line = line.split()
    trate_one.append(1/float(line[1])*10**12)
    kpoints1.append(float(line[0]))


lines = file_two.readlines()
trate_two = []
kpoints2 = []
for n, line in enumerate(lines):
    line = line.split()
    trate_two.append(1/float(line[1])*10**12)
    kpoints2.append(float(line[0]))


plt.figure()
plt.yscale('log')
plt.plot(kpoints1, trate_one, 'b-+')
plt.plot(kpoints2, trate_two, 'g-+')
plt.xlabel("N")
plt.ylabel("Decay time (ps)")
plt.legend(["No sp", "Sp"])
plt.title("Decay time w/ ribbon width (2 bands)")

plt.tight_layout()
plt.show()

