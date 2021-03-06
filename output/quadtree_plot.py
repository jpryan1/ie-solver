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
CMAP = copy(matplotlib.cm.magma)
CMAP.set_bad('lightgray', 1.)
MASKED_VALUE = 11111.1
fig = plt.figure(figsize=(14,14))


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
tree_lines = open("output/data/ie_solver_tree.txt","r").readlines()
boundary_lines = open("output/data/ie_solver_boundary.txt","r").readlines()

boundary_points = []
for line in boundary_lines:
	linesplit = line.split(',')
	boundary_points.append([float(linesplit[0]), float(linesplit[1])])
solution_points = []

quadtree_data = []
current_level_data = []
for line in tree_lines[1:]:
	if "level " in line:
		quadtree_data.append(current_level_data)
		current_level_data = []
	else:
		linesplit = line.split(',')
		current_level_data.append([float(linesplit[0]), float(linesplit[1]),
			float(linesplit[2]), float(linesplit[3]), float(linesplit[4])])
num_levels = len(quadtree_data)
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


def draw_box(img, datum):
	side_length = datum[2]
	bottom_left = datum[:2]

	pixel_data = datum[3]
	if(pixel_data == 0.):
		pixel_data = MASKED_VALUE
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

	for x in range(bottom_left[0]+1, top_right[0]):
		for y in range(bottom_left[1]+1, top_right[1]):
			img[x][y] = pixel_data


def draw_quadtree(img, quadtree_data):
	# get corners
	for datum in quadtree_data:
		draw_box(img, datum)

###########################################################
#
#
#					MAIN CODE
#
#
###########################################################
# TREE DECOMPOSITION PLOT
current_shown_level = 0
tree_images = []
for i in range(num_levels):
	tree_img = np.array([[MASKED_VALUE for x in range(WINDOW_SIZE)] for y in range(WINDOW_SIZE)])
	draw_quadtree(tree_img, quadtree_data[i])
	draw_boundary(tree_img, boundary_points, val=1.0)
	tree_img = np.ma.masked_where(tree_img == MASKED_VALUE, tree_img)
	tree_images.append(tree_img)

plt.title("Quadtree")
# im = plt.imshow(tree_images[current_shown_level].T, cmap=CMAP, vmin=0., vmax=1., origin = "lower")
im = plt.imshow(tree_images[current_shown_level].T, cmap=CMAP, vmin=0., vmax=1., origin = "lower")
plt.colorbar(im)

############################################################
#
#						INTERACTIVITY
#
############################################################

def on_key(event):
	global current_shown_level, im

	if(event.key == "down" and current_shown_level < num_levels-1):
		current_shown_level += 1
		im.set_data(tree_images[current_shown_level].T)		
		plt.draw()
	elif(event.key == "up" and current_shown_level > 0):
		current_shown_level -= 1
		im.set_data(tree_images[current_shown_level].T)
		plt.draw()


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

cid = fig.canvas.mpl_connect('key_press_event', on_key)
plt.savefig("quadtree_plot.png")
plt.show()
