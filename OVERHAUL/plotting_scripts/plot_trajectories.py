import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.stack([np.loadtxt(file, usecols=(1,5,6,7,8), skiprows=1) \
                 for file in sys.argv[1:]])

data = np.swapaxes(data, 0, 1)

plt.figure(figsize=(6,6), dpi=500)

pid = 7108

plt.plot( data[pid,:,0], data[pid,:,1] )
plt.plot( data[pid,:,0], data[pid,:,2] )
plt.plot( data[pid,:,0], data[pid,:,3] )
plt.plot( data[pid,:,0], data[pid,:,4] )

plt.savefig('particle_vs_t.png', dpi=500, bbox_inches='tight', pad_inches = 0)


'''
data = np.stack([np.loadtxt(file, usecols=(5,6,7,8), skiprows=1) \
                 for file in sys.argv[1:]])


data0 = data[0,:,0]
data = np.swapaxes(data, 0, 1)[np.where(data0 > 400.0)]

print(data.shape)

plt.figure(figsize=(4,4), dpi=500)

for particle in data:
    plt.plot( particle[:,1], particle[:,0] )

plt.savefig('T_vs_muB.png', dpi=500, bbox_inches='tight', pad_inches = 0)
print('Saved to T_vs_muB.png')

plt.figure(figsize=(4,4), dpi=500)

for particle in data:
    plt.plot( particle[:,2], particle[:,0] )

plt.savefig('T_vs_muS.png', dpi=500, bbox_inches='tight', pad_inches = 0)
print('Saved to T_vs_muS.png')

plt.figure(figsize=(4,4), dpi=500)

for particle in data:
    plt.plot( particle[:,3], particle[:,0] )

plt.savefig('T_vs_muQ.png', dpi=500, bbox_inches='tight', pad_inches = 0)
print('Saved to T_vs_muQ.png')
'''
