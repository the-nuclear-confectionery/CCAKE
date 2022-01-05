import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.stack([np.loadtxt(file, usecols=(5,6,7,8), skiprows=1) \
                 for file in sys.argv[1:]])

data0 = data[0,:,0]

data = np.swapaxes(data, 0, 1)[np.where(data0 > 484.0)]

print(data.shape)

plt.figure(figsize=(8,8), dpi=500)

for particle in data:
    plt.plot( particle[:,1], particle[:,0], 'r-' )

plt.savefig('T_vs_muB.png', dpi=500, bbox_inches='tight', pad_inches = 0)
