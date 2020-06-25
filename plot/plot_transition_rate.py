import matplotlib.pyplot as plt 
import numpy as np 

file_one = open("transition_dumb", "r")
file_two = open("transition_dumb", "r")

lines = file_one.readlines()
trate_one = []
trate_two = []
kpoints = []
for line in lines:
    line = line.split("\t")
    trate_one.append(float(line[1]))
    kpoints.append(float(line[0]))

lines = file_two.readlines()
for line in lines:
    line = line.split("\t")
    trate_two.append(float(line[1]))

fig, ax1 = plt.subplots()
ax1.plot(kpoints, trate_one, 'b+-')
ax1.set_xlabel("nk")
ax1.set_ylabel("Transition rate", color='tab:blue')

ax2 = ax1.twinx()
ax2.plot(kpoints, trate_two, 'g+-')
ax2.set_ylabel("Transition rate", color='tab:green')
fig.legend(["Same edge", "Different edges"], loc='upper center')

fig.tight_layout()
plt.show()

