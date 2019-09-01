import numpy as np
from copy import copy
import matplotlib
import matplotlib.pyplot as plt

# TODO: this file is highly targeted in dimensions, make more robust ASAP!

fig, axs = plt.subplots(2,1, figsize=(6,9))

n_scaling_lines = open("output/data/ie_solver_n_scaling.txt", "r").readlines()
e_scaling_lines = open("output/data/ie_solver_e_scaling.txt", "r").readlines()

n_x_vals = []
n_y_vals = []

for line in n_scaling_lines:
	linesplit = line.split(",")
	n_x_vals.append(int(linesplit[0]))
	n_y_vals.append(float(linesplit[1]))

e_x_vals = []
e_y_vals = []
for line in e_scaling_lines:
	linesplit = line.split(",")
	e_x_vals.append(float(linesplit[0]))
	e_y_vals.append(float(linesplit[1]))
	

# RUNTIME VS N PLOT

axs[0].set_xlabel("Number of discretization nodes")
axs[0].set_ylabel("Time (s)")
axs[0].set_title("Skel time vs N")
axs[0].loglog(n_x_vals, n_y_vals)

# RUNTIME VS EPSILON PLOT
# axs[1].set_xlabel("ID error tolerance")
# axs[1].set_ylabel("Time (s)")
# axs[1].semilogx(e_x_vals, e_y_vals)
# axs[1].set_title("Skel time vs ID tolerance")

plt.show()
