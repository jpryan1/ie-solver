import numpy as np
import sys
from copy import copy
import matplotlib
from mpl_toolkits.axes_grid1 import make_axes_locatable
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
# TODO: this file is highly targeted in dimensions, make more robust ASAP!

# Ex 2
is_channel_plot = True
ARROW_LENGTH = 0.4
BORDER_WIDTH = 5
HEAD_WIDTH = 5
QUIVER_RES_X = 20
QUIVER_RES_Y = 10
BOUNDARY_RES = 5
ZOOM = 1.65
OUTPUT_FILE = "ex2subplot.eps"
# config.num_boundary_points = pow(2, 14);
# config.domain_size = 200;

fig, axs = plt.subplots(1,3, figsize=(24,8))
MASKED_VALUE = 11111.1

CMAP = copy(matplotlib.cm.viridis)
CMAP.set_bad("white",1.)

for ax_idx in range(3):
	axs[ax_idx].axis("off")

	solution_lines = open("output/data/ex2_sol" + str(ax_idx+1)+".txt", "r").readlines()
	boundary_lines = open("output/data/ex2_bound" + str(ax_idx+1)+".txt", "r").readlines()
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
			mag = np.sqrt((linesplit[2])**2 + (linesplit[3])**2) 
			solution_grid[i][j] = mag if mag!=0 else MASKED_VALUE
			if(i % QUIVER_RES_X != 0 or j % QUIVER_RES_Y != 0 \
				or np.sqrt((linesplit[2])**2 + (linesplit[3])**2)<0.1):
				continue
			X.append((linesplit[0]))
			Y.append((linesplit[1]))
			U.append((linesplit[2]))
			V.append((linesplit[3]))

	solution_grid = np.ma.masked_where(solution_grid == MASKED_VALUE, solution_grid)
	imsh = axs[ax_idx].imshow(solution_grid,
		extent=[min_sol_x, max_sol_x, min_sol_y, max_sol_y], origin="lower", \
		cmap=CMAP, interpolation="bilinear")
	quiver_scale = (10 / ARROW_LENGTH ) * ZOOM
	axs[ax_idx].quiver(X,Y,U,V, color="white",scale=quiver_scale, headwidth=HEAD_WIDTH)

	# Boundary plot
	for i in range(0,len(boundary_lines)-BOUNDARY_RES,BOUNDARY_RES):
		pixel = boundary_lines[i].split(",")
		pixel = [float(pixel[0]), float(pixel[1])]
		next_pixel = boundary_lines[i+BOUNDARY_RES].split(",")
		next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
		if((pixel[0]-next_pixel[0])**2+(pixel[1]-next_pixel[1])**2>0.1**2):
			continue
		axs[ax_idx].plot([pixel[0], next_pixel[0]], [pixel[1],next_pixel[1]], \
			linewidth=BORDER_WIDTH, color="black")
	xl, xr = axs[ax_idx].get_xlim()
	yl, yr = axs[ax_idx].get_ylim()
	l = min(xl,yl)-0.01
	r = max(xr,yr)+0.01
	axs[ax_idx].set_xlim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)
	axs[ax_idx].set_ylim((l - (r+l)/2.)/ZOOM + (r+l)/2., (r - (r+l)/2.)/ZOOM + (r+l)/2.)



plt.savefig(OUTPUT_FILE, bbox_inches="tight", format="eps")
plt.show()
