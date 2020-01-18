import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys


font = {'size'   : 18}

matplotlib.rc('font', **font)
fig, axs = plt.subplots(1,2, figsize=(18,8))

skel_times_update = []
skel_times_noupdate = []
lines = open("scale_plots/ex2_data.txt").readlines()
slowlines = open("scale_plots/ex2_slow_data.txt").readlines()

for i in range(0,len(lines), 2):
  skel_times_update.append(float(lines[i].split(" ")[2]))


lines = open("scale_plots/ex2_data.txt").readlines()

for i in range(0,len(slowlines), 2):
  skel_times_noupdate.append(float(slowlines[i].split(" ")[2]))



LW=3
MS=10
axs[0].plot(range(len(skel_times_update)),skel_times_update,"r^", linewidth=LW, markersize=MS)

axs[0].set_xlabel("Update number")
axs[0].set_ylabel("Factor time (s)")
axs[0].set_title("Updating - 100 geometry changes")


elapsed_fast = [sum(skel_times_update[:(i+1)]) for i in range(len(skel_times_update))]
elapsed_slow = [sum(skel_times_noupdate[:(i+1)]) for i in range(len(skel_times_noupdate))]


axs[1].plot(range(len(elapsed_fast)),elapsed_fast,"r", label="Updating factorization", linewidth=LW)
axs[1].plot(range(len(elapsed_slow)),elapsed_slow,"b", label="Factorizing from scratch", linewidth=LW)

axs[1].set_xlabel("Update Number")
axs[1].set_ylabel("Elapsed Time (s)")
axs[1].set_title("Effect of updating on elapsed time")
axs[1].legend()


# for ax in axs:
#   ax.set_xscale('log', basex=2)
#   ax.set_yscale('log', basey=2)
  
plt.tight_layout()
plt.savefig("scale_plots/ex2_plot.eps", format="eps")
# plt.savefig("ex2_plot.png")