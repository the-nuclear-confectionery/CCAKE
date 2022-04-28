import numpy as np
import sys

# initialize arrays
print('Initializing arrays...')
prev_i = 0
previous_EoSs = np.genfromtxt(sys.argv[1], usecols=(-1), skip_header=1, dtype=str)
                           
# load all data
print('Loading all data...')
all_data = np.array([np.loadtxt(filename, usecols=tuple(range(1,6)), skiprows=1)
                     for filename in sys.argv[1:]])
print(all_data.shape)

for curr_i, filename in enumerate(sys.argv[2:], 1):
    print('  - differences between', prev_i, 'and', curr_i, ':')
    current_EoSs = np.genfromtxt(filename, usecols=(-1), skip_header=1, dtype=str)
    
    particles_to_print = np.where( current_EoSs != previous_EoSs )
    data_to_print = all_data[:, particles_to_print]
    
    for particle_to_print in particles_to_print:
        print(particle_to_print)
        print(all_data[:, particle_to_print].shape)
        print('\n\n\n')
        
    prev_i = curr_i
    previous_EoSs = current_EoSs