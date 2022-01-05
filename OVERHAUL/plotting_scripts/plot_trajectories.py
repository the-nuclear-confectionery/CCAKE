import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

data = np.stack([np.loadtxt(file, usecols=(5,6,7,8), skiprows=1) \
                 for file in sys.argv[1:]])

data = np.swapaxes(data, 0, 1)

print(data.shape)

data = data[np.where(data[0,0] > 150.0)]

print(data.shape)

#plt.figure(figsize=(8,8), dpi=500)

#plt.plot(  )