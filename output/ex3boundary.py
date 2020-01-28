import numpy as np
import sys
from copy import copy
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
# TODO: this file is highly targeted in dimensions, make more robust ASAP!


# Ex 3
is_channel_plot = False
ARROW_LENGTH = None
BORDER_WIDTH = 8
HEAD_WIDTH = None
QUIVER_RES_X = None
QUIVER_RES_Y = None
BOUNDARY_RES = 5
ZOOM = 1
TICK_LABEL_SIZE = 40
TICKS = [-1,0,1]
OUTPUT_FILE = "ex3a.eps"
# config.num_boundary_points = pow(2, 12);
# config.domain_size = 200;

print("args: {ZOOM} ")
fig, ax = plt.subplots(figsize=(14,14))
ax.axis("off")

MASKED_VALUE = 11111.1
if(len(sys.argv) > 1):
  ZOOM = float(sys.argv[1])

boundary_lines = open("output/data/ie_solver_boundary.txt", "r").readlines()


# Boundary plot

colors=["green", "orange", "blue"]
color_idx=0


for i in range(0,len(boundary_lines)-BOUNDARY_RES,BOUNDARY_RES):
	pixel = boundary_lines[i].split(",")
	pixel = [float(pixel[0]), float(pixel[1])]
	next_pixel = boundary_lines[i+BOUNDARY_RES].split(",")
	next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
	if((pixel[0]-next_pixel[0])**2+(pixel[1]-next_pixel[1])**2>0.1**2):
		color_idx+=1
		continue
	ax.plot([pixel[0], next_pixel[0]], [pixel[1],next_pixel[1]], \
		linewidth=BORDER_WIDTH, color=colors[color_idx])


xl, xr = ax.get_xlim()
yl, yr = ax.get_ylim()
l = min(xl,yl)-0.01
r = max(xr,yr)+0.01
print(l,r)
c = (r+l)/2.

ax.set_xlim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)
ax.set_ylim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)
ax.text(c-0.1,c+0.4,r"$\Omega$", fontsize=50, color="black")
ax.text(c+1,c+1,r"$\Gamma_0$", fontsize=50, color="green")
ax.text(c+0.6,c-0.6,r"$\Gamma_2$", fontsize=50, color="orange")
ax.text(c-0.8,c-0.6,r"$\Gamma_1$", fontsize=50, color="blue")

ax.scatter([c-0.84, c+0.77],[c-0.05,c-0.05], s=80, color="black")

plt.savefig("boundary.eps", bbox_inches="tight", format="eps")
plt.show()
