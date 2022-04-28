import numpy as np
import sys

# initialize arrays
print('Initializing arrays...')
prev_i = 0
previous_pids = np.genfromtxt(sys.argv[1], usecols=(0), skip_header=1, dtype=int)
previous_EoSs = np.genfromtxt(sys.argv[1], usecols=(-1), skip_header=1, dtype=str)
                           
# load all data
print('Loading all data...')
all_data = np.array([np.loadtxt(filename, usecols=tuple(range(2,17)), skiprows=1)
                     for filename in sys.argv[1:]])
print(all_data.shape)

for curr_i, filename in enumerate(sys.argv[2:]):
    print('  - differences between', prev_i, 'and', curr_i, ':')
    current_pids = np.genfromtxt(filename, usecols=(0), skip_header=1, dtype=int)
    current_EoSs = np.genfromtxt(filename, usecols=(-1), skip_header=1, dtype=str)
    
    particles_to_print = np.where( current_EoSs != previous_EoSs )
    data_to_print = all_data[:, particles_to_print]
    
    for particle_to_print in particles_to_print:
        print(particle_to_print)
        print(all_data[:, particle_to_print])
        print('\n\n\n')