import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

infilename = sys.argv[1]
outdirectory = sys.argv[2]

# HDF file structure
f = h5.File(infilename, 'r')
event = f['Event']
event_keys = list(event.keys())
n_timesteps = len(event_keys)

# Figure dimensions
width, height = 6, 6

# Axes ranges
#xmin, xmax, ymin, ymax = 2.5, 7.5, -7.5, -2.5

chosen_dpi = 200

#########################################################################################
def frame_to_array(i):
    frame = event[event_keys[i]]
    return [ frame.attrs['Time'], frame.attrs['e_2_X'], frame.attrs['e_2_P'] ]


#########################################################################################
def plot_eccentricities_vs_time():
    # Set up
    print('Building data...')
    data = np.array([frame_to_array(i) for i in range(n_timesteps)])
    
    #####################################
    # quantity vs. tau figure
    #####################################
    print('Plotting...')
    fig, ax = plt.subplots()
    fig.set_size_inches(width, height, forward=True)
    fig.set_dpi(chosen_dpi)
    
    ax.plot( data[:,0], data[:,1], 'b-', label=r'$\varepsilon_{2,X}$' )
    ax.plot( data[:,0], data[:,2], 'r--', label=r'$\varepsilon_{2,P}$' )
    
    plt.xlabel(r'$\tau$ (fm/$c$)')
        
    outfilename = outdirectory + '/eccentricities_vs_tau.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)



#########################################################################################
if __name__== "__main__":
    plot_eccentricities_vs_time()

