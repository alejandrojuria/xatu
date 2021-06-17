import matplotlib.pyplot as plt 
import numpy as np 

file_one = open("transition_ribbon_2_300_noonsite_dos", "r")
file_two = open("transition_ribbon_4_300_noonsite_dos", "r")
file_three = open("transition_ribbon_6_300_noonsite_dos", "r")
file_four = open("transition_ribbon_4_300_noonsite_dos_lr", "r")


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

lines = file_three.readlines()
trate_three = []
kpoints3 = []
for n, line in enumerate(lines):
    line = line.split()
    trate_three.append(1/float(line[1])*10**12)
    kpoints3.append(float(line[0]))

lines = file_four.readlines()
trate_four = []
kpoints4 = []
for n, line in enumerate(lines):
    line = line.split()
    trate_four.append(1/float(line[1])*10**12)
    kpoints4.append(float(line[0]))


plt.figure()
plt.yscale('log')
plt.plot(kpoints1, trate_one, 'b-+')
plt.plot(kpoints2, trate_two, 'g-+')
plt.plot(kpoints3, trate_three, 'r-+')
#plt.plot(kpoints2, trate_four, 'y-+')

plt.xlabel("N")
plt.ylabel("Decay time (ps)")
plt.legend(["2 bands", "4 bands", "6 bands", "600 nk"])
plt.title("Decay time w/ ribbon width (4 bands, 600 nk)")

plt.tight_layout()
plt.show()

