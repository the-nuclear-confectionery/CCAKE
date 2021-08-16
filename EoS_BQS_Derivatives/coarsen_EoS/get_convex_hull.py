import numpy as np
from scipy.spatial import ConvexHull

hc=197.327

def normalize(v):
    return np.amin(v), np.amax(v), (v-np.amin(v))/(np.amax(v)-np.amin(v))

data=np.loadtxt('../Thermodynamics_dense/EoS_Taylor_AllMu.dat', usecols=(0,1,2,3,6,7,8,9))
[T,muB,muQ,muS,b,s,q,e]=data.T
b*=T**3/hc**3; s*=T**3/hc**3; q*=T**3/hc**3; e*=T**4/hc**3;
loge=np.log(e); bt=b/(e/hc)**0.75; st=s/(e/hc)**0.75; qt=q/(e/hc)**0.75;
loge0,loge1,loge=normalize(loge); bt0,bt1,bt=normalize(bt); st0,st1,st=normalize(st); qt0,qt1,qt=normalize(qt)

hull = ConvexHull(np.c_[ loge, bt, st, qt ])
grid=np.c_[T,muB,muS,muQ,loge,bt,st,qt]

np.savetxt('hull.dat', grid[hull.vertices], fmt="%12.8f")