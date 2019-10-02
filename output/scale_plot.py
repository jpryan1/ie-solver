import matplotlib.pyplot as plt
import numpy as np
import sys

if(len(sys.argv)<2):
    print("Pass path as arg please")
    exit()

path = sys.argv[1]

num_dofs = []
num_skel_dofs = []
time_skel = []
time_solve = []

lines = open(path+"1.txt").readlines()
for line in lines:
  if "num_dofs" in line:
    num_dofs.append(int(line.split(" ")[1])/5.)
  if "num_skel_dofs" in line:
    num_skel_dofs.append(int(line.split(" ")[1])/5.)
  if "skeletonize" in line:
    time_skel.append(float(line.split(" ")[2])/5.)  
  if "solve" in line:
    time_solve.append(float(line.split(" ")[2])/5.)
  

for trial in range(2, 6):
  stat_counter = -1
  lines = open(path+str(trial)+".txt")
  for line in lines:
    if "num_dofs" in line:
      stat_counter += 1
      num_dofs[stat_counter] += int(line.split(" ")[1])/5.
    if "num_skel_dofs" in line:
      num_skel_dofs[stat_counter] += int(line.split(" ")[1])/5.
    if "skeletonize" in line:
      time_skel[stat_counter] += float(line.split(" ")[2])/5.
    if "solve" in line:
      time_solve[stat_counter] += float(line.split(" ")[2])/5.


# num_dofs = [num_dofs[i] for i in range(1, len(num_dofs),2)]
num_skel_dofs = [num_skel_dofs[i] for i in range(1, len(num_skel_dofs),2)]
time_skel = [time_skel[i] for i in range(1, len(time_skel),2)]
time_solve = [time_solve[i] for i in range(1, len(time_solve),2)]

fig, axs = plt.subplots(1,3, figsize=(20,10))
fig.suptitle("Updated, not parallelized")
axs[0].plot(num_dofs, num_skel_dofs)
axs[0].set_xlabel("Num dofs")
axs[0].set_ylabel("Num skel dofs")
axs[1].plot(num_dofs, time_skel)
axs[2].plot(num_dofs, time_solve)
axs[1].set_xlabel("Num dofs")
axs[1].set_ylabel("Seconds")
axs[2].set_xlabel("Num dofs")
axs[2].set_ylabel("Seconds")
# plt.show()
plt.savefig(path+"timeplot.png")