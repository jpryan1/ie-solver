
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
font = {'size'   : 18}

matplotlib.rc('font', **font)
fig, axs = plt.subplots(1,2, figsize=(19,7))

data = open("output/ex2grads.txt", "r").readlines()

pointsx = []
pointsy = []
vals = []
for line in data:
  spl = line.split(" ")

  pre = [float(s) for s in spl]
  post = []
  for num in pre[:2]:
    tmp = num
    while(tmp>2*np.pi):
      tmp -= 2*np.pi
    post.append(tmp)
  pointsx.append(post[0])
  pointsy.append(post[1])
  vals.append(pre[2])
  
for i in range(len(pointsx)):
  if(pointsx[i]<pointsy[i]):
    pointsx[i] = pointsx[i] + 2*np.pi

sc = axs[0].scatter(pointsx,pointsy,c=vals, cmap=plt.cm.cool, s=300)
cbar=fig.colorbar(sc, ax=axs[0])
cbar.set_ticks([-0.03, 0, 0.03])
# axs[0].colorbar()
axs[0].set_xticks([i*(np.pi/2.) for i in range(1,8)])
axs[0].set_xticklabels([r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$",r"$5\pi/2$",r"$3\pi$",r"$7\pi/2$"])
axs[0].set_yticks([i*(np.pi/2.) for i in range(0,5)])
axs[0].set_yticklabels(["0",r"$\pi/2$",r"$\pi$",r"$3\pi/2$",r"$2\pi$"])

axs[0].set_xlabel(r"$\theta_1$")
axs[0].set_ylabel(r"$\theta_2$")
axs[0].set_title("Experiment 2")


data = open("output/ex3flows.txt", "r").readlines()

angs = []
flows = []
for line in data:
  spl = line.split(" ")

  pre = [float(s) for s in spl]
  
  tmp = pre[0]
  while(tmp>2*np.pi):
    tmp -= 2*np.pi
  angs.append(tmp)
  flows.append(pre[1])

axs[1].set_xticks([i*(np.pi/4.) for i in range(5)])
axs[1].set_xticklabels(["0", r"$\pi/4$", r"$\pi/2$", r"$3\pi/4$",r"$\pi$"])

axs[1].plot(angs, flows)
axs[1].set_xlabel("Fin Angle")
axs[1].set_ylabel("Corridor Flow")

axs[1].set_title("Experiment 3")
plt.savefig("landscape.png")
plt.show()