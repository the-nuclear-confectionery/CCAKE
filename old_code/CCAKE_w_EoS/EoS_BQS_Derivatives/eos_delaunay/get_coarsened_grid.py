import numpy as np
from scipy.spatial import ConvexHull
from scipy.interpolate import griddata, LinearNDInterpolator

hc=197.327

def normalize(v):
    return np.amin(v), np.amax(v), (v-np.amin(v))/(np.amax(v)-np.amin(v))

data=np.loadtxt('../Thermodynamics_dense/EoS_Taylor_AllMu.dat', usecols=(0,1,2,3,6,7,8,9))
[T,muB,muQ,muS,b,s,q,e]=data.T
b*=T**3/hc**3; s*=T**3/hc**3; q*=T**3/hc**3; e*=T**4/hc**3;
loge=np.log(e); bt=b/(e/hc)**0.75; st=s/(e/hc)**0.75; qt=q/(e/hc)**0.75;
loge0,loge1,loge=normalize(loge); bt0,bt1,bt=normalize(bt); st0,st1,st=normalize(st); qt0,qt1,qt=normalize(qt)


hull = ConvexHull(np.c_[ loge, bt, st, qt ])
grid=np.c_[T,muB,muS,muQ,e,b,s,q,loge,bt,st,qt]
hullverts = grid[hull.vertices]

densgrid=np.unique(np.vstack((hullverts[:,[0,1,2,3,-4,-3,-2,-1]],grid[::1000,[0,1,2,3,-4,-3,-2,-1]])),axis=0)

interp = LinearNDInterpolator(densgrid[:,-4:], densgrid[:,0])

grid.shape
grid[15087080/2]
grid[15087080//2]
interp(8.62970736e-01,  6.21804422e-01,  2.89749791e-01,  3.38728023e-01)
interp(grid[0:10,-4:])
Test=interp(grid[:,-4:])
Test/T-1.0
np.amax(np.abs(Test/T-1.0))
np.isnan(Test).shape
np.where(np.isnan(Test)).shape
len(np.where(np.isnan(Test)))
np.where(np.isnan(Test))
np.where(~np.isnan(Test))
ok=np.where(~np.isnan(Test))
np.amax(np.abs(Test[ok]/T[ok]-1.0))