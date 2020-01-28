import matplotlib
matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt


font = {'size'   : 25}

matplotlib.rc('font', **font)

data = open("output/tst_ex4.txt", "r").readlines()
fig, ax = plt.subplots(1,1, figsize=(15,15))
angs = []
flows = []
for line in data:
  spl = line.split(" ")

  pre = [float(s) for s in spl[1:]]
  
  tmp = pre[0]
  while(tmp>2*np.pi):
    tmp -= 2*np.pi
  angs.append(tmp)
  flows.append(pre[1])
print(flows)
ax.set_xticks([i*(np.pi/2.) for i in range(5)])
ax.set_xticklabels(["0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$",r"$2\pi$"])

ax.plot(angs, flows, "ro", linewidth=10)
ax.set_xlabel("Source/sink angle")
ax.set_ylabel("Horizontal flow at origin ")

ax.set_title("Experiment 4")
plt.savefig("ex4landscape.png")
plt.savefig("ex4landscape.eps", format="eps")
# plt.show()
