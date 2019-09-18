#  Reads from output/bake/sol
#  Shows animation and writes to cwd as movie.mp4
import numpy as np
from copy import copy
from os import listdir
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure(figsize=(12,12))
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

WINDOW_SIZE = 140*5
IMAGE_SIZE = 100*5
CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad('lightgray', 1.)
MASKED_VALUE = 11111.1

###########################################################
#
#							READING THE FILES
#
###########################################################
#
#
#
#
#    TODO fix colorscheme across frames
#
num_files = len(listdir("output/bake/sol/"))

is_stokes = False
if(len(open("output/bake/sol/0.txt","r").readlines()[0].split(","))==4):
  is_stokes = True

# First, store all of the data in an array
files_boundary_points= []
for i in range(num_files):
  boundary_lines = open("output/bake/boundary/"+str(i)+".txt","r").readlines()
  boundary_points = []
  for line in boundary_lines:
    linesplit = line.split(',')
    boundary_points.append([float(linesplit[0]), float(linesplit[1])])
  files_boundary_points.append(boundary_points)
	
files_solution_points = []
for i in range(num_files):
  solution_lines = open("output/bake/sol/"+str(i)+".txt","r").readlines()
  solution_points = []
  for line in solution_lines:
    linesplit = line.split(',')
    if(is_stokes):
      solution_points.append([float(linesplit[0]), float(linesplit[1]),
        float(linesplit[2]), float(linesplit[3])])
    else:
      solution_points.append([float(linesplit[0]), float(linesplit[1]),
        float(linesplit[2])])
    #endif
  #endfor
  files_solution_points.append(solution_points)


###########################################################
#
#							SCALING THE PLOT
#
###########################################################
min_x = 0
min_y = 0
max_x = 1
max_y = 1

dif_x = max_x-min_x
dif_y = max_y-min_y
# Translation needed so that the image is centered
gamma=0
delta=0
if(dif_x>dif_y):
	delta = int((WINDOW_SIZE/2.0)*(1-(dif_y)/float(dif_x)))
else:
	gamma = int((WINDOW_SIZE/2.0)*(1-(dif_x)/float(dif_y)))
	
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
def draw_solution(img, points):
  for point in points:
    pixel = scaled_point(point[:2])
    if(np.isnan(point[2]) or point[2] == 0):
      img[pixel[0]][pixel[1]] = MASKED_VALUE
    else:
      img[pixel[0]][pixel[1]] = point[2]


def draw_boundary(img, points, val):
	for point in points:
		pixel = scaled_point(point)
		for r in range(-1, 2):
			for c in range(-1,2):
				img[pixel[0]+r][pixel[1]+c] = val

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



###########################################################
#
#
#					MAIN CODE
#
#
###########################################################


images = []
quivers = []
boundaries = []
for i in range(len(files_solution_points)):
  solution_points = files_solution_points[i]
  boundary_points = files_boundary_points[i]
  image = np.array([[MASKED_VALUE for x in range(WINDOW_SIZE)] for y in range(WINDOW_SIZE)])
  draw_boundary(image, boundary_points, val=1.0)
  if(is_stokes):
    stokes_data = get_quiver_data(solution_points)
    quivers.append(stokes_data)
  else:
    draw_solution(image, solution_points)
  image = np.ma.masked_where(image == MASKED_VALUE, image)

  images.append(image)

stokes_plot = 0
image_plot = plt.imshow(images[0].T, cmap=CMAP, animated=True, origin = "lower")
patches = [image_plot]
if is_stokes:
  stokes_data = quivers[0]
  stokes_plot = plt.quiver(stokes_data[0], stokes_data[1], stokes_data[2], stokes_data[3],
    stokes_data[4], cmap = "Purples")
  patches.append(stokes_plot)

idx = 0
# define updating function (just looks at next data in array, plots)
def animate(i):
    global idx, num_files, is_stokes
    idx += 1
    if idx == num_files:
      idx = 0
    if(is_stokes):
      stokes_data = quivers[idx]
      stokes_plot.set_UVC(stokes_data[2], stokes_data[3], stokes_data[4])
    image_plot.set_array(images[idx].T)
    return patches
ani = animation.FuncAnimation(fig, animate, interval=50, blit=True)
# ani.save('movie.mp4', writer=writer)
plt.show()

