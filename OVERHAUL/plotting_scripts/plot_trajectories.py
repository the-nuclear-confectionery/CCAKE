import numpy as np
import sys

data = np.stack([np.loadtxt(file, usecols=(5,6,7,8), skiprows=1) \
                 for file in sys.argv[1:]])

print(data.shape)