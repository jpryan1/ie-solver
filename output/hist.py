import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm



txt = open("cond.txt","r").readlines()

widths = []
size = []
conds = []
for line in txt:
    split = line.split(',')
    widths.append(float(split[0]))
    size.append(float(split[1]))
    conds.append(float(split[2]))
    

widths_min = min(widths)
widths_max = max(widths)
size_min = min(size)
size_max = max(size)
conds_min = min(conds)
conds_max = max(conds)
xs = []
ys = []
colors = []
areas = [2 for x in range(len(widths))]

for i in range(len(widths)):
    x = (widths[i] - widths_min)/float(widths_max-widths_min)
    y = (size[i] - size_min)/float(size_max-size_min)
    xs.append(widths[i])
    ys.append(size[i])
    color = cm.hot((conds[i]-conds_min)/float(conds_max-conds_min))
    colors.append(color)

plt.scatter(xs, ys, s=areas, c=colors)
plt.xlabel("Width")
plt.ylabel("Dofs")
plt.savefig("scatter.png")