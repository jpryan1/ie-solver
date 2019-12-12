import numpy as np
import sys
from copy import copy
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
# TODO: this file is highly targeted in dimensions, make more robust ASAP!

print("args: {ZOOM} ")
fig = plt.figure(figsize=(14,14))

CMAP = copy(matplotlib.cm.viridis)
CMAP.set_bad('white', 1.)
MASKED_VALUE = 11111.1
ZOOM = 1
if(len(sys.argv) > 1):
  ZOOM = int(sys.argv[1])

quiver_scale = 40./ZOOM

solution_lines = open("output/data/ie_solver_solution.txt", "r").readlines()
boundary_lines = open("output/data/ie_solver_boundary.txt", "r").readlines()

is_stokes = (len(solution_lines[0].split(",")) == 4)
	
solution_dim = int(np.sqrt(len(solution_lines)))
solution_grid = np.array([[MASKED_VALUE for x in range(solution_dim)] for y in range(solution_dim)])

X, Y, U, V = [], [], [], []
quiver_res = 5
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
			if(i % quiver_res != 0 or j % quiver_res != 0 \
				or np.sqrt((linesplit[2])**2 + (linesplit[3])**2)<0.1):
				continue
			X.append((linesplit[0]))
			Y.append((linesplit[1]))
			U.append((linesplit[2]))
			V.append((linesplit[3]))

solution_grid = np.ma.masked_where(solution_grid == MASKED_VALUE, solution_grid)
imsh = plt.imshow(solution_grid,
	extent=[min_sol_x, max_sol_x, min_sol_y, max_sol_y], origin="lower", interpolation="bilinear")
if(is_stokes):
	plt.quiver(X,Y,U,V, color="white",scale=quiver_scale, headwidth=3)

# Boundary plot
bdry_res=5
for i in range(0,len(boundary_lines)-bdry_res,bdry_res):
	pixel = boundary_lines[i].split(",")
	pixel = [float(pixel[0]), float(pixel[1])]
	next_pixel = boundary_lines[i+bdry_res].split(",")
	next_pixel = [float(next_pixel[0]), float(next_pixel[1])]
	if((pixel[0]-next_pixel[0])**2+(pixel[1]-next_pixel[1])**2>0.1**2):
		continue
	plt.plot([pixel[0], next_pixel[0]], [pixel[1],next_pixel[1]], \
		linewidth=8, color="black")

plt.axis("off")
# circle = plt.Circle(((xr-xl)/2.,(yr-yl)/2.), 0.5, color='g', fill=True, linewidth=2,linestyle="--")
# plt.gcf().gca().add_artist(circle)
# plt.savefig("ex3.png") #, format="eps")
plt.colorbar(imsh)
xl, xr = plt.xlim()
yl, yr = plt.ylim()
plt.xlim((xl - (xr+xl)/2.)/ZOOM + (xr+xl)/2., (xr - (xr+xl)/2.)/ZOOM + (xr+xl)/2.)  
plt.ylim((yl - (yr+yl)/2.)/ZOOM + (yr+yl)/2., (yr - (yr+yl)/2.)/ZOOM + (yr+yl)/2.)  
plt.savefig("comparewidth.png")
plt.show()
