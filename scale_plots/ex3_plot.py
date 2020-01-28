import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys


font = {'size'   : 18}

matplotlib.rc('font', **font)
fig, axs = plt.subplots(1,2, figsize=(20,8))


times_fast = []
times_med = []
times_slow = []

timing_lines = open("scale_plots/ex3Atiming.txt","r").readlines()

fast_lines = timing_lines[1:(26*3+1)]
med_lines = timing_lines[(26*3+2):(26*6+2)]
slow_lines = timing_lines[(26*6+3):(26*9+3)]

N=25

for i in range(N):
  val = float(fast_lines[1+i])+float(fast_lines[27+i])+float(fast_lines[53+i])
  times_fast.append(val/3.)
  val = float(med_lines[1+i])+float(med_lines[27+i])+float(med_lines[53+i])
  times_med.append(val/3.)
  val = float(slow_lines[1+i])+float(slow_lines[27+i])+float(slow_lines[53+i])
  times_slow.append(val/3.)


LW=3
MS=10

axs[0].plot(range(N), \
  [sum(times_fast[:i]) for i in range(1,N+1)], \
  ":g^",label="Scheme A",\
  linewidth=LW, markersize=MS)
axs[0].plot(range(N), \
  [sum(times_med[:i]) for i in range(1,N+1)], \
  ":b^", label="Scheme B",\
  linewidth=LW, markersize=MS)
axs[0].plot(range(N), \
  [sum(times_slow[:i]) for i in range(1,N+1)], \
  ":r^",label="Scheme C",\
  linewidth=LW, markersize=MS)

# axs[0].set_title("Efficient multithreading for solution gradient estimation")
axs[0].set_title("Example 3a")
axs[0].legend()
axs[0].set_xlabel("Optimization step")
axs[0].set_ylabel("Elapsed time (s)")
  
  #####################################################
  
  
times_fast = []
times_med = []
times_slow = []

timing_lines = open("scale_plots/ex3Btiming.txt","r").readlines()

fast_lines = timing_lines[1:(26*3+1)]
med_lines = timing_lines[(26*3+2):(26*6+2)]
slow_lines = timing_lines[(26*6+3):(26*9+3)]

N=25

for i in range(N):
  val = float(fast_lines[1+i])+float(fast_lines[27+i])+float(fast_lines[53+i])
  times_fast.append(val/3.)
  val = float(med_lines[1+i])+float(med_lines[27+i])+float(med_lines[53+i])
  times_med.append(val/3.)
  val = float(slow_lines[1+i])+float(slow_lines[27+i])+float(slow_lines[53+i])
  times_slow.append(val/3.)


LW=3
MS=10

axs[1].plot(range(N), \
  [sum(times_fast[:i]) for i in range(1,N+1)], \
  ":g^",label="Scheme A",\
  linewidth=LW, markersize=MS)
axs[1].plot(range(N), \
  [sum(times_med[:i]) for i in range(1,N+1)], \
  ":b^", label="Scheme B",\
  linewidth=LW, markersize=MS)
axs[1].plot(range(N), \
  [sum(times_slow[:i]) for i in range(1,N+1)], \
  ":r^",label="Scheme C",\
  linewidth=LW, markersize=MS)

# axs[1].set_title("Efficient multithreading for solution gradient estimation")
axs[1].set_title("Example 3b")

axs[1].legend()
axs[1].set_xlabel("Optimization step")
axs[1].set_ylabel("Elapsed time (s)")
  
  
  
  
  
plt.tight_layout()
# plt.savefig("ex3_plot.png")
plt.savefig("scale_plots/ex3_plot.eps", format="eps")
