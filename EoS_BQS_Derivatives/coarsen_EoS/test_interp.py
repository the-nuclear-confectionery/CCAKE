import numpy as np
from scipy.interpolate import interp1d, LinearNDInterpolator

grid=np.loadtxt('grid_ebsq_mag.dat')

zero_density = grid[np.where((grid[:,2]==0)&(grid[:,3]==0)&(grid[:,4]==0))];
slice_e = grid[np.where((grid[:,5]>=79990)&(grid[:,5]<=80010))]
slice_e.shape
#(280, 9)
print(np.amin(slice_e[:,1]), np.amax(slice_e[:,1]))
#445.0 470.0
slice_T = grid[np.where((grid[:,1]>=445)&(grid[:,1]<=475))]
slice_T.shape
#(354571, 9)

interp1D = interp1d(zero_density[:,-4], zero_density[:,1])
interp = LinearNDInterpolator(slice_T[:,-4:], slice_T[:,1], rescale=True)

interp1D(np.arange(79990,80010))
#array([470.13131512, 470.13275921, 470.1342033 , 470.13564739,
#       470.13709149, 470.13853558, 470.13997967, 470.14142376,
#       470.14286785, 470.14431194, 470.14575604, 470.14720013,
#       470.14864422, 470.15008831, 470.1515324 , 470.15297649,
#       470.15442059, 470.15586468, 470.15730877, 470.15875286])

np.array(list(map( lambda x : interp((x,0,0,0)), np.arange(79990,80010))))
#array([470.13131512, 470.13275921, 470.1342033 , 470.13564739,
#       470.13709149, 470.13853558, 470.13997967, 470.14142376,
#       470.14286785, 470.14431194, 470.14575604, 470.14720013,
#       470.14864422, 470.15008831, 470.1515324 , 470.15297649,
#       470.15442059, 470.15586468, 470.15730877, 470.15875286])