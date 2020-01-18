import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys


font = {'size'   : 18}

matplotlib.rc('font', **font)
fig, ax = plt.subplots(1,1, figsize=(10,8))

num_dofs = []
times_dumb_fast = []
times_dumb_med = []
times_dumb_slow = []

times_fast = []
times_med = []
times_slow = []


# dumb_fast1 = open("scale_plots/ex3_dumb_fast_1.txt").readlines()
# num_times = len(dumb_fast1)
# dumb_fast2 = open("scale_plots/ex3_dumb_fast_2.txt").readlines()
# dumb_fast3 = open("scale_plots/ex3_dumb_fast_3.txt").readlines()
# for i in range(num_times):
#   val1 = float(dumb_fast1[i].split(" ")[2])
#   val2 = float(dumb_fast2[i].split(" ")[2])
#   val3 = float(dumb_fast3[i].split(" ")[2])
#   times_dumb_fast.append((val1+val2+val3)/3.)

# dumb_med1 = open("scale_plots/ex3_dumb_med_1.txt").readlines()
# dumb_med2 = open("scale_plots/ex3_dumb_med_2.txt").readlines()
# dumb_med3 = open("scale_plots/ex3_dumb_med_3.txt").readlines()
# for i in range(len(dumb_med1)):
#   val1 = float(dumb_med1[i].split(" ")[2])
#   val2 = float(dumb_med2[i].split(" ")[2])
#   val3 = float(dumb_med3[i].split(" ")[2])
#   times_dumb_med.append((val1+val2+val3)/3.)

# dumb_slow1 = open("scale_plots/ex3_dumb_slow_1.txt").readlines()
# dumb_slow2 = open("scale_plots/ex3_dumb_slow_2.txt").readlines()
# dumb_slow3 = open("scale_plots/ex3_dumb_slow_3.txt").readlines()
# for i in range(len(dumb_slow1)):
#   val1 = float(dumb_slow1[i].split(" ")[2])
#   val2 = float(dumb_slow2[i].split(" ")[2])
#   val3 = float(dumb_slow3[i].split(" ")[2])
#   times_dumb_slow.append((val1+val2+val3)/3.)


fast1 = open("scale_plots/ex3_fast_1.txt").readlines()
num_times=len(fast1)
fast2 = open("scale_plots/ex3_fast_2.txt").readlines()
fast3 = open("scale_plots/ex3_fast_3.txt").readlines()
for i in range(num_times):
  val1 = float(fast1[i].split(" ")[2])
  val2 = float(fast2[i].split(" ")[2])
  val3 = float(fast3[i].split(" ")[2])
  times_fast.append((val1+val2+val3)/3.)

med1 = open("scale_plots/ex3_med_1.txt").readlines()
med2 = open("scale_plots/ex3_med_2.txt").readlines()
med3 = open("scale_plots/ex3_med_3.txt").readlines()
for i in range(num_times):
  val1 = float(med1[i].split(" ")[2])
  val2 = float(med2[i].split(" ")[2])
  val3 = float(med3[i].split(" ")[2])
  times_med.append((val1+val2+val3)/3.)
 
slow1 = open("scale_plots/ex3_slow_1.txt").readlines()
slow2 = open("scale_plots/ex3_slow_2.txt").readlines()
slow3 = open("scale_plots/ex3_slow_3.txt").readlines()
for i in range(num_times):
  val1 = float(slow1[i].split(" ")[2])
  val2 = float(slow2[i].split(" ")[2])
  val3 = float(slow3[i].split(" ")[2])
  times_slow.append((val1+val2+val3)/3.)





LW=3
MS=10

ax.plot(range(num_times), \
  [sum(times_fast[:i]) for i in range(1,num_times+1)], \
  ":g^",label="Scheme A",\
  linewidth=LW, markersize=MS)
ax.plot(range(num_times), \
  [sum(times_med[:i]) for i in range(1,num_times+1)], \
  ":b^", label="Scheme B",\
  linewidth=LW, markersize=MS)
ax.plot(range(num_times), \
  [sum(times_slow[:i]) for i in range(1,num_times+1)], \
  ":r^",label="Scheme C",\
  linewidth=LW, markersize=MS)

ax.set_title("Efficient multithreading for solution gradient estimation")
ax.legend()
ax.set_xlabel("Optimization step")
ax.set_ylabel("Elapsed time (s)")
  
plt.tight_layout()
# plt.savefig("ex3_plot.png")
plt.savefig("scale_plots/ex3_plot.eps", format="eps")
