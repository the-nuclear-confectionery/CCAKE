import numpy as np
from scipy.interpolate import interp1d, LinearNDInterpolator

grid=np.loadtxt('grid_ebsq_mag.dat')

zero_density = grid[np.where((grid[:,2]==0)&(grid[:,3]==0)&(grid[:,4]==0))]
slice_e = grid[np.where((grid[:,5]>=70)&(grid[:,5]<=80))]
slice_T = grid[np.where((grid[:,1]>=45)&(grid[:,1]<=125))]

interp = LinearNDInterpolator(slice_T[:,-4:], slice_T[:,1], rescale=True)
interp1D = interp1d(zero_density[:,-4], zero_density[:,1])

print('4D')
np.array(list(map( lambda x : interp((x,0,0,0)), np.arange(70,80,0.25))))
print('1D')
interp1D(np.arange(70,80,0.25))