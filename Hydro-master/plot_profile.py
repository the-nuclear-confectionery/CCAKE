import matplotlib.pyplot as plt
import matplotlib.tri as tri
import numpy as np
from scipy.interpolate import griddata
import sys

#filename = sys.argv[1]
#colsToPlot = tuple(map(int,sys.argv[2:]))
filename = "C:/Users/Christopher Plumberg/Desktop/Research/UIUC/BSQ/Hydro-master/bsqsveprofile8_ev0.dat"
colsToPlot = (2, 3, 5)

data = np.loadtxt(filename, usecols=colsToPlot, skiprows=1)
[x, y, f] = data.T

#extent = x0, x1, y0, y1 = np.min(x), np.max(x), np.min(y), np.max(y)
extent = x0, x1, y0, y1 = -3.0, 3.0, -3.0, 3.0
#plt.subplots(nrows=1, ncols=1) # same as two lines below
fig = plt.figure()
ax = plt.gca()


#############
# version 1 #
#############

# Linearly interpolate the data (x, y) on a grid defined by (xi, yi).
xi = np.linspace(x0, x1, 1000)
yi = np.linspace(y0, y1, 1000)
triang = tri.Triangulation(x, y)
#interpolator = tri.LinearTriInterpolator(triang, f)
interpolator = tri.CubicTriInterpolator(triang, f, kind='geom')
#interpolator = tri.CubicTriInterpolator(triang, f, kind='min_E')
Xi, Yi = np.meshgrid(xi, yi)
zi = interpolator(Xi, Yi)

#dims = zi.shape
#zi = zi.flatten()
#zi[ np.where( zi is np.ma.masked ) ] = 0.0
#zi = zi.reshape(dims)
plt.imshow(zi, cmap=plt.cm.viridis, interpolation='bicubic', extent=extent)


'''
#################################################
# version 2 (basically same output as version 1 #
#################################################

xi = np.linspace(x0, x1, 1000)
yi = np.linspace(y0, y1, 1000)
Xi, Yi = np.meshgrid(xi, yi)
grid = griddata(data[:,[0,1]], data[:,2], (Xi, Yi), method='cubic')

plt.imshow(grid.T, cmap=plt.cm.viridis, interpolation='bicubic', extent=extent)
'''

ax.set_facecolor('black')

plt.show()