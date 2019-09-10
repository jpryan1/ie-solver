#  Reads from output/data - solution, tree, and boundary
#  saves as domain_plot.png
import numpy as np
import sys
from copy import copy
import matplotlib
import matplotlib.pyplot as plt
# TODO: this file is highly targeted in dimensions, make more robust ASAP!

WINDOW_SIZE = 140*5
IMAGE_SIZE = 100*5
CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad('lightgray', 1.)
MASKED_VALUE = 11111.1
fig, axs = plt.subplots(2,1, figsize=(8,12))


ZOOM = 1
if(len(sys.argv) > 1):
  ZOOM = int(sys.argv[1])
SHIFT = 0
if(len(sys.argv) > 2):
  SHIFT = int(sys.argv[2])*ZOOM

CENTER = WINDOW_SIZE/2.0
###########################################################
#
#							READING THE FILES
#
###########################################################
solution_lines = open("output/data/ie_solver_solution.txt","r").readlines()
tree_lines = open("output/data/ie_solver_tree.txt","r").readlines()
boundary_lines = open("output/data/ie_solver_boundary.txt","r").readlines()

boundary_points = []
for line in boundary_lines:
	linesplit = line.split(',')
	boundary_points.append([float(linesplit[0]), float(linesplit[1])])
solution_points = []
is_stokes = False
if len(solution_lines[0].split(",")) == 4:
	is_stokes = True
for line in solution_lines:
	linesplit = line.split(',')
	if(is_stokes):
		solution_points.append([float(linesplit[0]), float(linesplit[1]),
			float(linesplit[2]), float(linesplit[3])])
	else:
		solution_points.append([float(linesplit[0]), float(linesplit[1]),
			float(linesplit[2])])
# centerx centery sidelenth
quadtree_data = []
for line in tree_lines:
	linesplit = line.split(',')
	quadtree_data.append([float(linesplit[0]), float(linesplit[1]),
		float(linesplit[2]), linesplit[3], linesplit[4]])

###########################################################
#
#							SCALING THE PLOT
#
###########################################################
# Calculate the min and max point on boundary so we can scale properly
min_x = boundary_points[0][0]
max_x = boundary_points[0][0]
min_y = boundary_points[0][1]
max_y = boundary_points[0][1]
for pair in boundary_points:
	if pair[0] < min_x:
		min_x = pair[0]
	if pair[0] > max_x:
		max_x = pair[0]
	if pair[1] < min_y:
		min_y = pair[1]
	if pair[1] > max_y:
		max_y = pair[1]

dif_x = max_x-min_x
dif_y = max_y-min_y

# Translation needed so that the image is centered
# Note from John: I do not understand why this works, I should probably
# get on that
gamma=0
delta=0
if(dif_x>dif_y):
	delta = int((IMAGE_SIZE/2.0)*(1-(dif_y)/float(dif_x)))
else:
	gamma = int((IMAGE_SIZE/2.0)*(1-(dif_x)/float(dif_y)))
	
scale_factor = IMAGE_SIZE/max(dif_x,dif_y)
#	The points will undergoes a dilation and translation so that the
#	bounding box is [20,120]x[20,120].
def scaled_point(point):
	x = int(np.round( (point[0] - min_x)*(scale_factor)))
	y = int(np.round( (point[1] - min_y)*(scale_factor)))
	x += gamma + int((WINDOW_SIZE-IMAGE_SIZE)/2.0)
	y += delta + int((WINDOW_SIZE-IMAGE_SIZE)/2.0)
	return [x, y]

#
############################################################
#
#						DRAWING THE PLOT
#
############################################################
#
def draw_boundary(img, points, val):
	for point in points:
		pixel = scaled_point(point)
		
		for r in range(-1, 2):
			for c in range(-1,2):
			  x_zoom = (pixel[0] - CENTER)*ZOOM + CENTER + SHIFT
			  y_zoom = (pixel[1] - CENTER)*ZOOM + CENTER
			  x_coord = max(0,min(WINDOW_SIZE-1, x_zoom+r))
			  y_coord = max(0,min(WINDOW_SIZE-1, y_zoom+c))
			  img[int(x_coord)][int(y_coord)] = val


def draw_solution(img, points):
	for point in points:
		pixel = scaled_point(point[:2])
		if(np.isnan(point[2]) or point[2] == 0):
			img[pixel[0]][pixel[1]] = MASKED_VALUE
		else:
			for i in range(-8,9):
				for j in range(-8,9):
				  
				  x_zoom = (pixel[0] - CENTER)*ZOOM + CENTER + SHIFT
				  y_zoom = (pixel[1] - CENTER)*ZOOM + CENTER
				  if x_zoom+i < 0 or x_zoom+i > WINDOW_SIZE-1 or y_zoom+j < 0 or y_zoom+j > WINDOW_SIZE-1:
				    continue
				  x_coord = max(0,min(WINDOW_SIZE-1, x_zoom+i))
				  y_coord = max(0,min(WINDOW_SIZE-1, y_zoom+j))
				  img[int(x_coord)][int(y_coord)] = point[2]
				  


def get_quiver_data(points):
	# returns an array containing the four vecs necessary for a quiver plot
	# This can be replaced with vecs with colon index probably TODO
	X = []
	Y = []
	U = []
	V = []
	colors = []
	for point in points:
		pixel = scaled_point(point[:2])
		x_zoom = (pixel[0]  - CENTER)*ZOOM + CENTER+SHIFT
		y_zoom = (pixel[1]  - CENTER)*ZOOM + CENTER
		if x_zoom < 0 or x_zoom > WINDOW_SIZE-1 or y_zoom < 0 or y_zoom > WINDOW_SIZE-1:
		  continue
		x_coord = max(0,min(WINDOW_SIZE-1, x_zoom))
		y_coord = max(0,min(WINDOW_SIZE-1, y_zoom))
		
		X.append(x_coord)
		Y.append(y_coord)
		U.append(point[2])
		V.append(point[3])
		colors.append(point[2]**2 + point[3]**2)
	return [X, Y, U, V, colors]

def draw_box(img, side_length, bottom_left):
	top_right = [bottom_left[0] + side_length, bottom_left[1] + side_length]
	bottom_left = scaled_point(bottom_left)
	top_right = scaled_point(top_right)
	for x in range(bottom_left[0], top_right[0]+1):
		for y in [bottom_left[1], top_right[1]]:
			for r in range(-1, 2):
				for c in range(-1,2):
				  x_coord = max(0, min(x+r, len(img)-1))
				  y_coord = max(0, min(y+c, len(img[0])-1))
				  img[x_coord][y_coord] = 1.0

	for x in [bottom_left[0], top_right[0]]:
		for y in range(bottom_left[1], top_right[1]+1):
			for r in range(-1, 2):
				for c in range(-1,2):
				  x_coord = max(0, min(x+r, len(img)-1))
				  y_coord = max(0, min(y+c, len(img[0])-1))
				  img[x_coord][y_coord] = 1.0

def draw_quadtree(img, quadtree_data):
	# get corners
	for datum in quadtree_data:
		side_length = datum[2]
		bottom_left = datum[:2]
		draw_box(img, side_length, bottom_left)

###########################################################
#
#
#					MAIN CODE
#
#
###########################################################
# TREE DECOMPOSITION PLOT
tree_img = np.array([[MASKED_VALUE for x in range(WINDOW_SIZE)] for y in range(WINDOW_SIZE)])
draw_quadtree(tree_img, quadtree_data)
draw_boundary(tree_img, boundary_points, val=1.0)
tree_img = np.ma.masked_where(tree_img == MASKED_VALUE, tree_img)
axs[1].set_title("Quadtree")
axs[1].imshow(tree_img.T, cmap=CMAP, origin = "lower")

# SOLUTION PLOT
solution_img = np.array([[MASKED_VALUE for x in range(WINDOW_SIZE)] for y in range(WINDOW_SIZE)])
stokes_data = []
if(is_stokes):
	draw_boundary(solution_img, boundary_points, val=1.0)
	stokes_data = get_quiver_data(solution_points)
else:
	# draw_boundary(solution_img, boundary_points, val=1.0)
	draw_solution(solution_img, solution_points)
solution_img = np.ma.masked_where(solution_img == MASKED_VALUE, solution_img)
axs[0].set_title("Solution")
axs[0].imshow(solution_img.T, cmap=CMAP, origin = "lower")
if(is_stokes):
	axs[0].quiver(stokes_data[0], stokes_data[1], stokes_data[2], stokes_data[3],
		stokes_data[4], cmap = "Purples", #cmap='autumn',
		norm=matplotlib.colors.Normalize(vmin=0,vmax=1.), scale=20./ZOOM)

############################################################
#
#						INTERACTIVITY
#
############################################################
# def onclick(event):
# 	print("******************************\nCLICKED!")
# 	x = event.xdata
# 	y = event.ydata
# 	for datum in quadtree_data:
# 		qx, qy = scaled_point(datum[:2])
# 		if x-qx < scale_factor*datum[2] and x-qx>0 and y-qy>0 and y-qy < scale_factor*datum[2]:
# 			print(str(qx)+ ","+str(qy)+" SIZE: "+str(datum[2])+" ID: " + str(datum[3]) + " DOFS: "+str(datum[4]))
#   # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
#   #       ('double' if event.dblclick else 'single', event.button,
#   #        event.x, event.y, event.xdata, event.ydata))

# cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.savefig("domain_plot.png")
plt.show()
