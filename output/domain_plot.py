import numpy as np
import sys
from copy import copy
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
# TODO: this file is highly targeted in dimensions, make more robust ASAP!

# Ex 1
# is_channel_plot = False
# ARROW_LENGTH = 0.4
# BORDER_WIDTH = 8
# HEAD_WIDTH = 3
# QUIVER_RES_X = 10
# QUIVER_RES_Y = 10
# BOUNDARY_RES = 5
# ZOOM = 1
# TICK_LABEL_SIZE = 40
# TICKS = [0.2, 0.9, 1.6]
# OUTPUT_FILE = "ex1.eps"
# config.num_boundary_points = pow(2, 12);
# config.domain_size = 200;

# Ex 2 
is_channel_plot = True
ARROW_LENGTH = 0.4
BORDER_WIDTH = 8
HEAD_WIDTH = 5
QUIVER_RES_X = 20
QUIVER_RES_Y = 10
BOUNDARY_RES = 1
ZOOM = 1.65
TICK_LABEL_SIZE = 20
TICKS = [0.25, 1, 1.75]
OUTPUT_FILE = "ex2.eps"
# config.num_boundary_points = pow(2, 12);
# config.domain_size = 200;

# Ex 3
# is_channel_plot = False
# ARROW_LENGTH = None
# BORDER_WIDTH = 8
# HEAD_WIDTH = None
# QUIVER_RES_X = None
# QUIVER_RES_Y = None
# BOUNDARY_RES = 5
# ZOOM = 1
# TICK_LABEL_SIZE = 40
# TICKS = [-1,0,1]
# OUTPUT_FILE = "ex3a.eps"
# config.num_boundary_points = pow(2, 12);
# config.domain_size = 200;


# Ex 4
# is_channel_plot = False
# ARROW_LENGTH = 0.4
# BORDER_WIDTH = 8
# HEAD_WIDTH = 3
# QUIVER_RES_X = 10
# QUIVER_RES_Y = 10
# BOUNDARY_RES = 5
# ZOOM = 1
# TICK_LABEL_SIZE = 40
# TICKS = [0.15, 0.75, 1.35]
# OUTPUT_FILE = "ex3b.eps"
# config.num_boundary_points = pow(2, 12);
# config.domain_size = 400;


print("args: {ZOOM} ")
fig, ax = plt.subplots(figsize=(14,14))
ax.axis("off")

MASKED_VALUE = 11111.1
if(len(sys.argv) > 1):
  ZOOM = float(sys.argv[1])


solution_lines = open("output/data/ie_solver_solution.txt", "r").readlines()
boundary_lines = open("output/data/ie_solver_boundary.txt", "r").readlines()

is_stokes = (len(solution_lines[0].split(",")) == 4)
if is_stokes:
	CMAP = copy(matplotlib.cm.viridis)
else:
	CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad("white",1.)
solution_dim = int(np.sqrt(len(solution_lines)))
solution_grid = np.array([[MASKED_VALUE for x in range(solution_dim)] for y in range(solution_dim)])

X, Y, U, V = [], [], [], []
min_sol_x, min_sol_y, max_sol_x, max_sol_y = 10,10,-10,-10
for i in range(solution_dim):
	for j in range(solution_dim):
		linesplit = [float(n) for n in solution_lines[i+solution_dim*j].split(',')]
		min_sol_x = min(min_sol_x, (linesplit[0]))
		max_sol_x = max(max_sol_x, (linesplit[0]))
		min_sol_y = min(min_sol_y, (linesplit[1]))
		max_sol_y = max(max_sol_y, (linesplit[1]))

		mag = np.sqrt((linesplit[2])**2 + (linesplit[3])**2) if is_stokes \
			else linesplit[2]
		solution_grid[i][j] = mag if mag!=0 else MASKED_VALUE
		if(is_stokes):
			if(i % QUIVER_RES_X != 0 or j % QUIVER_RES_Y != 0 \
				or np.sqrt((linesplit[2])**2 + (linesplit[3])**2)<0.1):
				continue
			X.append((linesplit[0]))
			Y.append((linesplit[1]))
			U.append((linesplit[2]))
			V.append((linesplit[3]))

solution_grid = np.ma.masked_where(solution_grid == MASKED_VALUE, solution_grid)
imsh = ax.imshow(solution_grid,
	extent=[min_sol_x, max_sol_x, min_sol_y, max_sol_y], origin="lower", \
	cmap=CMAP, interpolation="bilinear") # for ex3a: , vmin=-1.5, vmax=1.5)
if(is_stokes):
	quiver_scale = (10 / ARROW_LENGTH ) * ZOOM
	ax.quiver(X,Y,U,V, color="white",scale=quiver_scale, headwidth=HEAD_WIDTH)

# Boundary plot
for i in range(0,len(boundary_lines)-BOUNDARY_RES,BOUNDARY_RES):
	pixel = boundary_lines[i].split(",")
	pixel = [float(pixel[0]), float(pixel[1])]
	next_pixel = boundary_lines[i+BOUNDARY_RES].split(",")
	next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
	if((pixel[0]-next_pixel[0])**2+(pixel[1]-next_pixel[1])**2>0.1**2):
		continue
	ax.plot([pixel[0], next_pixel[0]], [pixel[1],next_pixel[1]], \
		linewidth=BORDER_WIDTH, color="black")

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.5)
cbar= plt.colorbar(imsh, cax=cax, ticks=TICKS)
cbar.ax.tick_params(labelsize=TICK_LABEL_SIZE)

# circle = plt.Circle(((xr-xl)/2.,(yr-yl)/2.), 0.5, color='g', fill=True, linewidth=2,linestyle="--")
# plt.gcf().gca().add_artist(circle)
# plt.savefig("ex3.png") #, format="eps")
xl, xr = ax.get_xlim()
yl, yr = ax.get_ylim()
l = min(xl,yl)-0.01
r = max(xr,yr)+0.01
ax.set_xlim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)

if(is_channel_plot):
	yl-=0.2
	yr+=0.2
	ax.set_ylim((yl - (yr+yl)/2.)/ZOOM + (yr+yl)/2., (yr - (yr+yl)/2.)/ZOOM + (yr+yl)/2.)
else:
	ax.set_ylim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)

plt.savefig(OUTPUT_FILE, bbox_inches="tight", format="eps")
plt.show()
