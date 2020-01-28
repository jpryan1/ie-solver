import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys


font = {'size'   : 18}

matplotlib.rc('font', **font)
fig, axs = plt.subplots(1,2, figsize=(20,8))


ex3alines = open("scale_plots/ex3Aconvergence.txt", "r").readlines()
valsa = [float(line.split(" ")[5]) for line in ex3alines]


ex3blines = open("scale_plots/ex3Bconvergence.txt", "r").readlines()
valsb = [float(line.split(" ")[5]) for line in ex3blines]

besta=max(valsa)
bestb=max(valsb)

diffsa = [besta-val for val in valsa]
diffsb = [bestb-val for val in valsb]


LW=3
MS=10

axs[0].plot(range(len(valsa)), diffsa,  linewidth=LW, markersize=MS)
axs[1].plot(range(len(valsb)), diffsb,  linewidth=LW, markersize=MS)

axs[0].set_xlim(0,25)
axs[1].set_xlim(0,25)


axs[0].set_title("Example 3a Convergence")
axs[0].set_xlabel("Optimization step")
axs[0].set_ylabel(r"$x^*-x_n$")
axs[0].set_yscale('log', basex=10)

axs[1].set_title("Example 3b Convergence")
axs[1].set_xlabel("Optimization step")
axs[1].set_ylabel(r"$x^*-x_n$")
axs[1].set_yscale('log', basex=10)

plt.savefig("scale_plots/ex3convergence.eps", bbox_inches="tight", format="eps")
#plt.savefig("ex3tmp.png")
#plt.show()
