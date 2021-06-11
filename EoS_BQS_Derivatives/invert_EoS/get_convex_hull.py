'''import numpy as np
from scipy.spatial import ConvexHull


# add HDF functionality eventually
data       = np.loadtxt('EoS_small_grid.dat', usecols=(0,9,6,7,8))

hc         = 197.327
data[:,1] *= 0.001*data[:,0]**4/hc**3   # possibly change e to have same units as the rest?
data[:,2] *= data[:,0]**3/hc**3
data[:,3] *= data[:,0]**3/hc**3
data[:,4] *= data[:,0]**3/hc**3
data       = data[:,1:]
'''


from scipy.spatial import ConvexHull, Delaunay
import numpy as np
from numpy.linalg import det
from scipy.stats import dirichlet


import numpy as np
from scipy.optimize import linprog

def in_hull(points, x):
    n_points = len(points)
    n_dim = len(x)
    c = np.zeros(n_points)
    A = np.r_[points.T,np.ones((1,n_points))]
    b = np.r_[x, np.ones(1)]
    lp = linprog(c, A_eq=A, b_eq=b)
    return lp.success


def dist_in_hull(points, n):
    dims = points.shape[-1]
    hull = points[ConvexHull(points).vertices]
    deln = points[Delaunay(hull).simplices]

    vols = np.abs(det(deln[:, :dims, :] - deln[:, dims:, :])) / np.math.factorial(dims)    
    sample = np.random.choice(len(vols), size = n, p = vols / vols.sum())

    return np.einsum('ijk, ij -> ik', deln[sample], dirichlet.rvs([1]*(dims + 1), size = n))

data = np.random.randn(10,2)

import matplotlib.pyplot as plt

plt.scatter(data[:,0], data[:,1], color = 'blue')

hull = ConvexHull(data)




newdata = dist_in_hull(data, 100)

#plt.scatter(newdata[:,0], newdata[:,1], color = 'red')

minima = np.min(data, axis=0)
maxima = np.max(data, axis=0)

newdata2 = (np.mgrid[minima[0]:maxima[0]:10j,minima[1]:maxima[1]:10j]).reshape([2,10*10]).T

newdata2 = newdata2[list(map(lambda x: in_hull(data, x), newdata2))]

plt.scatter(newdata2[:,0], newdata2[:,1], color = 'purple')

hv = data[hull.vertices]
print(hv)

hv = np.append(hv, np.array([hv[0]]), axis=0)
#plt.plot(hv[:,0], hv[:,1], color = 'green', ls = '-')


plt.triplot(hv[:,0], hv[:,1], Delaunay(hv).simplices)




plt.show()

# construct convex hull of complete grid
#hull = ConvexHull(data)

# add some points (say, 50000) in the interior and Delaunay triangulate the whole mess
newdata = data[np.unique(np.concatenate((hull.vertices,np.random.choice(len(data),50000))))]
tri = Delaunay(newdata)




