import matplotlib.pyplot as plt
import numpy as np

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
    post.append(int(np.round(tmp*100)))
  pointsx.append(post[0])
  pointsy.append(post[1])
  vals.append(pre[2])
  
for i in range(len(pointsx)):
  if(pointsx[i]<pointsy[i]):
    pointsx[i] = pointsx[i] + 2*np.pi*100


plt.scatter(pointsx,pointsy,c=vals, cmap=plt.cm.cool, s=300)
plt.colorbar()
plt.xlabel("Ang1 * 100")
plt.ylabel("Ang2 * 100")
plt.show()