import numpy as np
from scipy.interpolate import interp1d, LinearNDInterpolator

#===============================================================================================
def get_interpolation1D(e0, rat, grid):
    print('In get_interpolation1D', flush=True)
    slice_e = grid[np.where(np.abs(grid[:,5] - e0) < rat*(np.amax(grid[:,5])-np.amin(grid[:,5])))]
    print('slice_e.shape =', slice_e.shape, flush=True)
    T0, T1 = np.amin(slice_e[:,1]), np.amax(slice_e[:,1])
    print('slice_T.shape =', slice_T.shape, T0, T1, flush=True)
    slice_T = grid[np.where((T0 <= grid[:,1]) & (grid[:,1] <= T1))]
    f = interp1d(grid[:,-4], grid[:,1])
    return f(e0)

#===============================================================================================
def get_interpolation4D(e0, b0, s0, q0, rat, grid):
    print('In get_interpolation4D', flush=True)
    grid = grid[np.where(np.abs(grid[:,5] - e0) < rat*(np.amax(grid[:,5])-np.amin(grid[:,5])))]
    [erange, brange, srange, qrange] \
             = rat*[np.amax(grid[:,5])-np.amin(grid[:,5]),
                    np.amax(grid[:,6])-np.amin(grid[:,6]),
                    np.amax(grid[:,7])-np.amin(grid[:,7]),
                    np.amax(grid[:,8])-np.amin(grid[:,8])]
    ce = np.abs(grid[:,5] - e0) <= erange
    cb = np.abs(grid[:,6] - b0) <= brange
    cs = np.abs(grid[:,7] - s0) <= srange
    cq = np.abs(grid[:,8] - q0) <= qrange
    slice_ebsq = grid[np.where(ce & cb & cs & cq)]
    print('slice_ebsq.shape =', slice_ebsq.shape, flush=True)
    T0, T1, mub0, mub1, mus0, mus1, muq0, muq1 \
        = np.amin(slice_ebsq[:,1]), np.amax(slice_ebsq[:,1]), \
          np.amin(slice_ebsq[:,2]), np.amax(slice_ebsq[:,2]), \
          np.amin(slice_ebsq[:,3]), np.amax(slice_ebsq[:,3]), \
          np.amin(slice_ebsq[:,4]), np.amax(slice_ebsq[:,4])
    cT   = (T0 <= grid[:,1]) & (grid[:,1] <= T1)
    cmub = (mub0 <= grid[:,2]) & (grid[:,2] <= mub1)
    cmus = (mus0 <= grid[:,3]) & (grid[:,3] <= mus1)
    cmuq = (muq0 <= grid[:,4]) & (grid[:,4] <= muq1)
    slice_Tmubmusmuq = grid[np.where(cT & cmub & cmus & cmuq)]
    print('slice_Tmubmusmuq.shape =', slice_Tmubmusmuq.shape, \
          T0, T1, mub0, mub1, mus0, mus1, muq0, muq1, flush=True)
    f = LinearNDInterpolator(slice_Tmubmusmuq[:,-4:], slice_Tmubmusmuq[:,1:5], rescale=True)
    return f(e0,b0,s0,q0)

#===============================================================================================
# load grid
print('Loading grid file...', flush=True)
grid=np.loadtxt('grid_ebsq_mag.dat')
grid[:,[3,4]] = grid[:,[4,3]]  # muS <--> muQ

zero_density = grid[np.where((grid[:,2]==0)&(grid[:,3]==0)&(grid[:,4]==0))]

print('4D')
np.array(list(map( lambda x : get_interpolation4D(x, 0, 0, 0, 0.01, grid), np.arange(79500,80500,50))))
print('1D')
np.array(list(map( lambda x : get_interpolation1D(x, 0, 0, 0, 0.01, grid), np.arange(79500,80500,50))))

for x in np.arange(79500,80500,50):
    print( x, get_interpolation1D(x, 0.01, grid), get_interpolation4D(x, 0, 0, 0, 0.01, grid), flush=True )
