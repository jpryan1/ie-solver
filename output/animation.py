import numpy as np
from copy import copy
from os import listdir
import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig = plt.figure()
Writer = animation.writers['ffmpeg']
writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

IMG_SIZE = 140
CMAP = copy(matplotlib.cm.hot)
CMAP.set_bad('lightgray', 1.)
MASKED_VALUE = 11111.1
min1 = 0
min2 = 0
max1 = 1
max2 = 1

def scaled_point(point):
  inner_size = IMG_SIZE-40.0
  dif1 = max1-min1
  dif2 = max2-min2
  x = int(np.round( (point[0] - min1)*(inner_size/dif1) + 20))
  y = int(np.round( (point[1] - min2)*(inner_size/dif2) + 20))
  #delete this soon please
  x = max(0,min(139,x))
  y = max(0,min(139,y))
  return [x, y]


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


def draw_solution(img, points):
  for point in points:
    pixel = scaled_point(point[:2])
    if(np.isnan(point[2]) or point[2] == 0):
      img[pixel[0]][pixel[1]] = MASKED_VALUE
    else:
      img[pixel[0]][pixel[1]] = point[2]

solution_images = []

num_files = len(listdir("output/bake/sol/"))


is_stokes = False
if(len(open("output/bake/sol/0.txt","r").readlines()[0].split(","))==4):
  is_stokes = True


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

  if(is_stokes):
    stokes_data = get_quiver_data(solution_points)
    solution_images.append(stokes_data)
  else:
    img = np.array([[MASKED_VALUE for x in range(IMG_SIZE)] for y in range(IMG_SIZE)])
    draw_solution(img, solution_points)
    # draw_boundary(img1, boundary_points, val=MASKED_VALUE)
    img = np.ma.masked_where(img == MASKED_VALUE, img)
    solution_images.append(img)

im = 0
if is_stokes:
  stokes_data = solution_images[0]
  im = plt.quiver(stokes_data[0], stokes_data[1], stokes_data[2], stokes_data[3], 
    stokes_data[4], cmap = "Purples")
else: 
  im = plt.imshow(solution_images[0].T, cmap=CMAP, animated=True)

idx = 0

def updatefig(*args):
    global idx, num_files, is_stokes
    idx += 1

    if idx == num_files:
      idx = 0
    
    if(is_stokes):
      stokes_data = solution_images[idx]
      im.set_UVC(stokes_data[2], stokes_data[3], stokes_data[4]) 
    else:
      im.set_array(solution_images[idx].T)
    return im,

ani = animation.FuncAnimation(fig, updatefig, interval=100, blit=True)
ani.save('movie.mp4', writer=writer)
plt.show()

