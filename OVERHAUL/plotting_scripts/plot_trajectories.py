import numpy as np
import sys

data = np.stack([np.loadtxt(file, usecols=()) for file in sys.argv[1:]])

print(data.shape)