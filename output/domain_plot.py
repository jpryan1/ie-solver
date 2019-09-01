import numpy as np
from copy import copy
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# TODO: this file is highly targeted in dimensions, make more robust ASAP!

IMG_SIZE = 140*5
CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad('lightgray', 1.)
MASKED_VALUE = 11111.1
is_stokes = False
fig, axs = plt.subplots(2,1, figsize=(6,9))

###########################################################
#
#
#
#
#							READING THE FILES
#
#
#
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
#
#							SCALING THE PLOT
#
#
###########################################################
# Calculate the min and max point on boundary so we can scale properly
min1 = boundary_points[0][0]
max1 = boundary_points[0][0]
min2 = boundary_points[0][1]
max2 = boundary_points[0][1]
for pair in boundary_points:
	if pair[0] < min1:
		min1 = pair[0]
	if pair[0] > max1:
		max1 = pair[0]
	if pair[1] < min2:
		min2 = pair[1]
	if pair[1] > max2:
		max2 = pair[1]

inner_size = IMG_SIZE-40.0*5
dif1 = max1-min1
dif2 = max2-min2
gamma=0
delta=0
if(dif1>dif2):
	delta = int((IMG_SIZE/2.0)*(1-(dif2)/float(dif1)))
else:
	gamma = int((IMG_SIZE/2.0)*(1-(dif1)/float(dif2)))
scale_factor = inner_size/max(dif1,dif2)
#	The points will undergoes a dilation and translation so that the
#	bounding box is [20,120]x[20,120]. 
def scaled_point(point, bound=True):
	inner_size = IMG_SIZE-40.0*5
	x = int(np.round( (point[0] - min1)*(scale_factor) + 20*5))
	y = int(np.round( (point[1] - min2)*(scale_factor) + 20*5))
	x += gamma
	y += delta
	#delete this soon please
	if(bound):
		x = max(0,min(139*5,x))
		y = max(0,min(139*5,y))
	return [x, y]

#
############################################################
#
#						DRAWING THE PLOT
#
#
############################################################
#
def draw_boundary(img, points, val):
	for point in points:
		pixel = scaled_point(point)
		
		for r in range(-1, 2):
			for c in range(-1,2):
				img[pixel[0]+r][pixel[1]+c] = val


def draw_solution(img, points):
	for point in points:
		pixel = scaled_point(point[:2])
		if(np.isnan(point[2]) or point[2] == 0):
			img[pixel[0]][pixel[1]] = MASKED_VALUE
		else:
			img[pixel[0]][pixel[1]] = point[2]
	 

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
		X.append(pixel[0])
		Y.append(pixel[1])
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
		
					img[x+r][y+c] = 1.0

	for x in [bottom_left[0], top_right[0]]:
		for y in range(bottom_left[1], top_right[1]+1):
			for r in range(-1, 2):
				for c in range(-1,2):
		
					img[x+r][y+c] = 1.0

def draw_quadtree(img, quadtree_data):
	# get corners
	for datum in quadtree_data:
		side_length = datum[2]
		bottom_left = datum[:2]
		draw_box(img, side_length, bottom_left)

###########################################################
#
#
#
#
#					MAIN CODE
#
#
#
#
###########################################################
# TREE DECOMPOSITION PLOT
tree_img = np.array([[MASKED_VALUE for x in range(IMG_SIZE)] for y in range(IMG_SIZE)])
draw_quadtree(tree_img, quadtree_data)
draw_boundary(tree_img, boundary_points, val=1.0)
tree_img = np.ma.masked_where(tree_img == MASKED_VALUE, tree_img)
axs[1].set_title("Quadtree")
axs[1].imshow(tree_img.T, cmap=CMAP, origin = "lower")

# SOLUTION PLOT
solution_img = np.array([[MASKED_VALUE for x in range(IMG_SIZE)] for y in range(IMG_SIZE)])
stokes_data = []
if(is_stokes):
	draw_boundary(solution_img, boundary_points, val=1.0)
	stokes_data = get_quiver_data(solution_points)
else:
	draw_boundary(solution_img, boundary_points, val=MASKED_VALUE)
	draw_solution(solution_img, solution_points)
solution_img = np.ma.masked_where(solution_img == MASKED_VALUE, solution_img)
axs[0].set_title("Solution")
axs[0].imshow(solution_img.T, cmap=CMAP, origin = "lower")
if(is_stokes):
	axs[0].quiver(stokes_data[0], stokes_data[1], stokes_data[2], stokes_data[3], 
		stokes_data[4], cmap = "Purples")

############################################################
#
#
#						INTERACTIVITY
#
#
############################################################
def onclick(event):
	print("******************************\nCLICKED!")
	x = event.xdata
	y = event.ydata
	for datum in quadtree_data:
		qx, qy = scaled_point(datum[:2], bound=False)
		if x-qx < scale_factor*datum[2] and x-qx>0 and y-qy>0 and y-qy < scale_factor*datum[2]:
			print(str(qx)+ ","+str(qy)+" SIZE: "+str(datum[2])+" ID: " + str(datum[3]) + " DOFS: "+str(datum[4]))
  # print('%s click: button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
  #       ('double' if event.dblclick else 'single', event.button,
  #        event.x, event.y, event.xdata, event.ydata))

cid = fig.canvas.mpl_connect('button_press_event', onclick)
plt.savefig("domain_plot.png")
#plt.show()
