#  Reads from output/bake/sol
#  Shows animation and writes to cwd as movie.mp4
import sys
import numpy as np
from copy import copy
from os import listdir
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure(figsize=(12,12))
RESOLUTION = 3
WINDOW_SIZE = 140*5
IMAGE_SIZE = 100*5
CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad('lightgray', 1.)
MASKED_VALUE = 11111.1

print("args: {ZOOM} {X_SHIFT}")
ZOOM = 1
if(len(sys.argv) > 1):
  ZOOM = int(sys.argv[1])
SHIFT = 0
if(len(sys.argv) > 2):
  SHIFT = int(sys.argv[2])*ZOOM

quiver_normalizer = matplotlib.colors.Normalize(vmin=0,vmax=1.)
quiver_scale = 40./ZOOM

CENTER = WINDOW_SIZE/2.0

Writer = animation.writers['ffmpeg']
writer = Writer(fps=30, metadata=dict(artist='Me'), bitrate=1800)

###########################################################
#
#							READING THE FILES
#
###########################################################
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
min_x = files_boundary_points[0][0][0]
max_x = files_boundary_points[0][0][0]
min_y = files_boundary_points[0][0][1]
max_y = files_boundary_points[0][0][1]
for pair in files_boundary_points[0]:
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
gamma=0
delta=0
if(dif_x>dif_y):
  delta = int((IMAGE_SIZE/2.0)*(1-(dif_y)/float(dif_x)))
else:
  gamma = int((IMAGE_SIZE/2.0)*(1-(dif_x)/float(dif_y)))
  
scale_factor = IMAGE_SIZE/max(dif_x,dif_y)
# The points will undergoes a dilation and translation so that the
# bounding box is [20,120]x[20,120].
def scaled_point(point):
  x = int(np.round( (point[0] - min_x)*(scale_factor)))
  y = int(np.round( (point[1] - min_y)*(scale_factor)))
  x += gamma + int((WINDOW_SIZE-IMAGE_SIZE)/2.0)
  y += delta + int((WINDOW_SIZE-IMAGE_SIZE)/2.0)
  return [x, y]

#
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
      for i in range(-RESOLUTION+1,RESOLUTION):
        for j in range(-RESOLUTION+1,RESOLUTION):
          
          x_zoom = (pixel[0] - CENTER)*ZOOM + CENTER + SHIFT
          y_zoom = (pixel[1] - CENTER)*ZOOM + CENTER
          if (x_zoom+i < 0 or x_zoom+i > WINDOW_SIZE-1 or 
              y_zoom+j < 0 or y_zoom+j > WINDOW_SIZE-1):
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
    if (x_zoom < 0 or x_zoom > WINDOW_SIZE-1 or 
        y_zoom < 0 or y_zoom > WINDOW_SIZE-1):
      continue
    x_coord = max(0,min(WINDOW_SIZE-1, x_zoom))
    y_coord = max(0,min(WINDOW_SIZE-1, y_zoom))
    
    X.append(x_coord)
    Y.append(y_coord)
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
    stokes_data[4], cmap = "Purples", norm=quiver_normalizer, scale=quiver_scale)
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

ani = animation.FuncAnimation(fig, animate, interval=150, blit=True)
# ani.save('testnewpxy.mp4', writer=writer)
plt.show()

