import numpy as np
from operator import itemgetter
import sys

quantities = ['id','tau','x','y','p','T','muB','muS','muQ','e','rhoB','rhoS','rhoQ','s']
cols = dict(zip(quantities,range(len(quantities))))
chosenQuantities = ['tau','e','rhoB','rhoS','rhoQ','s']
chosenCols = itemgetter(*chosenQuantities)(cols)

# initialize arrays
print('Initializing arrays...')
prev_i = 0
previous_EoSs = np.genfromtxt(sys.argv[1], usecols=(-1), skip_header=1, dtype=str)
                           
# load all data
print('Loading all data...')
all_data = np.array([np.loadtxt(filename, skiprows=1, usecols=chosenCols)
                     for filename in sys.argv[1:]])

np.set_printoptions(precision=6, suppress=True)

for curr_i, filename in enumerate(sys.argv[2:], 1):
    print('  - differences between', prev_i, 'and', curr_i, ':')
    current_EoSs = np.genfromtxt(filename, usecols=(-1), skip_header=1, dtype=str)
    
    particles_to_print = np.array(np.where( current_EoSs != previous_EoSs )).flatten()
    print(particles_to_print)
    
    for particle_to_print in particles_to_print:
        print('-------------')
        print(particle_to_print)
        print(*chosenQuantities)
        print(str(all_data[:, particle_to_print]).replace(' [', '').replace('[', '').replace(']', ''))
        print('\n')
        
    prev_i = curr_i
    previous_EoSs = current_EoSs