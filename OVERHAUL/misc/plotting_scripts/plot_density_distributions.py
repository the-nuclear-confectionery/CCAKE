import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm                                                                          
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import numpy as np
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable


# HDF file structure
f = h5.File(sys.argv[1], 'r')
event = f['Event']
n_timesteps = min([len(event.keys()),1000])
event_keys = list(event.keys())

# Figure dimensions
width, height = 6, 4

# Axes ranges
xmin, xmax, ymin, ymax = 2.5, 7.5, -7.5, -2.5
#xmin, xmax, ymin, ymax = -12, 12, -12, 12

chosen_colormap = cm.get_cmap('cool', 256)

chosen_dpi = 200

#########################################################################################
def frame_to_array(i, quantity):
    frame = event[event_keys[i]]
    q = np.array(frame[quantity])
    return np.c_[ np.full_like(q, frame.attrs['Time']), q ]


#########################################################################################
def plot_density_distribution_vs_time():
    # Set up
    print('Building data...')
    data = np.stack([frame_to_array(i, 'e') for i in range(n_timesteps)])
    
    data = data.reshape([len(data)//2,2])
    
    H, yedges, xedges = np.histogram2d(data[:,0], data[:,1])
    
    print(data.shape)
    
    #####################################
    # T vs. muB figure
    #####################################
    print('Plotting...')
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    ax1.pcolormesh(xedges, yedges, H, cmap='rainbow')
    
    outfilename = 'T_vs_muB.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)



#########################################################################################
if __name__== "__main__":
    plot_density_distribution_vs_time()



