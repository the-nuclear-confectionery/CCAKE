import numpy as np
from scipy.spatial import ConvexHull

hc=197.327

def normalize(v):
    return np.amin(v), np.amax(v), (v-np.amin(v))/(np.amax(v)-np.amin(v))

# set up
data=np.loadtxt('../Thermodynamics_staggered/EoS_Taylor_AllMu.dat', usecols=(0,1,2,3,6,7,8,9))
[T,muB,muQ,muS,b,s,q,e]=data.T
b*=T**3/hc**3; s*=T**3/hc**3; q*=T**3/hc**3; e*=T**4/hc**3;
loge=np.log(e); bt=b/(e/hc)**0.75; st=s/(e/hc)**0.75; qt=q/(e/hc)**0.75;

# e, b, s, q; actual values
#hull = ConvexHull(np.c_[ e, b, s, q ])
grid=np.c_[range(len(T)),T,muB,muS,muQ,e,b,s,q]
np.savetxt('grid_staggered_ebsq_mag.dat', grid, fmt=("%d" + " %12.8f"*8))
#np.savetxt('hull_staggered_ebsq_mag.dat', grid[hull.vertices], fmt=("%d" + " %12.8f"*8))

print('hull = ConvexHull(np.c_[ loge, bt, st, qt ]); mag')
hull = ConvexHull(np.c_[ loge, bt, st, qt ])
grid=np.c_[range(len(T)),T,muB,muS,muQ,loge,bt,st,qt]
np.savetxt('grid_staggered_logebtstqt_mag.dat', grid, fmt=("%d" + " %12.8f"*8))
np.savetxt('hull_staggered_logebtstqt_mag.dat', grid[hull.vertices], fmt=("%d" + " %12.8f"*8))

loge0,loge1,loge=normalize(loge); bt0,bt1,bt=normalize(bt); st0,st1,st=normalize(st); qt0,qt1,qt=normalize(qt)
e0,e1,e=normalize(e); b0,b1,b=normalize(b); s0,s1,s=normalize(s); q0,q1,q=normalize(q)

print('hull = ConvexHull(np.c_[ e, b, s, q ]); norm')
hull = ConvexHull(np.c_[ e, b, s, q ])
grid=np.c_[range(len(T)),T,muB,muS,muQ,e,b,s,q]
np.savetxt('grid_staggered_ebsq_norm.dat', grid, fmt=("%d" + " %12.8f"*8))
np.savetxt('hull_staggered_ebsq_norm.dat', grid[hull.vertices], fmt=("%d" + " %12.8f"*8))

print('hull = ConvexHull(np.c_[ loge, bt, st, qt ]); norm')
hull = ConvexHull(np.c_[ loge, bt, st, qt ])
grid=np.c_[range(len(T)),T,muB,muS,muQ,loge,bt,st,qt]
np.savetxt('grid_staggered_logebtstqt_norm.dat', grid, fmt=("%d" + " %12.8f"*8))
np.savetxt('hull_staggered_logebtstqt_norm.dat', grid[hull.vertices], fmt=("%d" + " %12.8f"*8))

