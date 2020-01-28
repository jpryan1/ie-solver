import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys


font = {'size'   : 18}

matplotlib.rc('font', **font)
fig, axs = plt.subplots(1,2, figsize=(18,8))

num_dofs = []
time_skel_noupdate = []
time_solve_noupdate = []

time_skel_update = []
time_solve_update = []

lines = open("scale_plots/ex1_data.txt").readlines()

for i in range(16,21):
  j=i-16
  num_dofs.append(2**i)
  trial_lines = lines[(2+13*j):(2+13*(j+1))]
  trial_nums = [float(line.split(" ")[2]) for line in trial_lines[1:]]
  time_skel_noupdate.append((trial_nums[0] + trial_nums[4] + trial_nums[8])/3.)
  time_solve_noupdate.append((trial_nums[1] + trial_nums[5] + trial_nums[9])/3.)
  
  time_skel_update.append((trial_nums[2] + trial_nums[6] + trial_nums[10])/3.)
  time_solve_update.append((trial_nums[3] + trial_nums[7] + trial_nums[11])/3.)


LW=3
MS=10
axs[0].plot(num_dofs,time_skel_noupdate,":k^", label="Factor time", linewidth=LW, markersize=MS)
axs[0].plot(num_dofs,[0.00002*x for x in num_dofs],":b^", label="t ~ n", linewidth=LW, markersize=MS)
axs[0].plot(num_dofs,[(0.00005*x)**2 for x in num_dofs],":r^", label="t ~ n^2", linewidth=LW, markersize=MS)
axs[0].set_xlabel("Number of points")
axs[0].set_ylabel("Time (s)")
axs[0].set_title("Factorization time")
axs[0].legend()

axs[1].plot(num_dofs,time_solve_noupdate,":k^", label="Solve time", linewidth=LW, markersize=MS)
axs[1].plot(num_dofs,[0.0000007*x for x in num_dofs],":b^", label="t ~ n", linewidth=LW, markersize=MS)
axs[1].plot(num_dofs,[(0.00001*x)**2 for x in num_dofs],":r^", label="t ~ n^2", linewidth=LW, markersize=MS)
axs[1].set_xlabel("Number of points")
axs[1].set_ylabel("Time (s)")
axs[1].set_title("Linear solve time")
axs[1].legend()



for ax in axs:
  ax.set_xscale('log', basex=10)
  ax.set_yscale('log', basey=10)
  
  
plt.tight_layout()
plt.savefig("scale_plots/SM_scale_plot.eps",  bbox_inches="tight", format="eps")
# plt.savefig("SM_scale_plot.png")
plt.show()