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


axs[0].plot(num_dofs,time_skel_noupdate, ":r^", label="Factor from scratch", linewidth=LW, markersize=MS)
axs[0].plot(num_dofs,time_skel_update, ":b^", label="Update factorization", linewidth=LW, markersize=MS)
# axs[1,0].plot(num_dofs,[0.000003*x for x in num_dofs],":b^", label="y ~ x", linewidth=LW, markersize=MS)
# axs[1,0].plot(num_dofs,[(0.00002*x)**2 for x in num_dofs],":r^", label="y ~ x^2", linewidth=LW, markersize=MS)
axs[0].legend()
axs[0].set_xlabel("Number of points")
axs[0].set_ylabel("Time (s)")
axs[0].set_title("Factoring speedup after update")
axs[0].set_xticks([i*100000 for i in range(1,10,2)])
axs[0].set_xticklabels(["100k", "300k", "500k", "700k", "900k", ])


axs[1].plot(num_dofs,time_solve_noupdate,":r^", label="Solve w/ original factorization", linewidth=LW, markersize=MS)

axs[1].plot(num_dofs,time_solve_update, ":b^", label="Solve w/ updated factorization", linewidth=LW, markersize=MS)
# axs[1,1].plot(num_dofs,[0.000003*x for x in num_dofs],":b^", label="y ~ x", linewidth=LW, markersize=MS)
# axs[1,1].plot(num_dofs,[(0.00002*x)**2 for x in num_dofs],":r^", label="y ~ x^2", linewidth=LW, markersize=MS)
axs[1].legend()
axs[1].set_xlabel("Number of points")
axs[1].set_ylabel("Time (s)")
axs[1].set_title("Linear solve with initial/updated factorization")
axs[1].set_xticks([i*100000 for i in range(1,10,2)])
axs[1].set_xticklabels(["100k", "300k", "500k", "700k", "900k", ])


for ax in axs:
  ax.set_xscale('log', basex=10)
  ax.set_yscale('log', basey=10)
  
  
plt.tight_layout()
# plt.savefig("ex1_speedup_plot.png")
plt.savefig("scale_plots/ex1_speedup_plot.eps", format="eps")
plt.show()