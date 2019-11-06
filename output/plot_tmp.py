import matplotlib.pyplot as plt


ns = [67500, 125000, 250000, 500000]

plt.title("Big skel block LU time")
plt.plot(ns, [0.19, 0.23, 0.35, 0.61], "-o")
plt.show()
