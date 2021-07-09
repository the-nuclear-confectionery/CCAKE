import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
import sys

filename = sys.argv[1]
outfilename = sys.argv[2]
colsToPlot = tuple(map(int,sys.argv[3:]))

data = np.loadtxt(filename, usecols=colsToPlot)
[x, y, f] = data.T

extent = np.min(x), np.max(x), np.min(y), np.max(y)
#plt.subplots(nrows=1, ncols=1) # same as two lines below
fig = plt.figure()
ax = plt.gca()

length = int(np.sqrt(f.size))
plt.imshow(f.reshape(length, length), cmap=plt.cm.viridis, interpolation='bicubic', extent=extent)
#ax.set_facecolor('black')

plt.show()