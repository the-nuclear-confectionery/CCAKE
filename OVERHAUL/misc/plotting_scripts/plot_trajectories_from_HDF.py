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

selection = None

def get_selection(Tmin, Tmax):
    frame    = event[event_keys[0]]
    T0       = np.array(frame['T'])
    eos_tags = np.array(frame['labels'])
    #print(np.amin(T0),np.amax(T0))
    return np.where( (T0 >= Tmin) & (T0 <= Tmax) & (eos_tags == 0) )


def frame_to_array(i, Tmin, Tmax):
    frame = event[event_keys[i]]
    tau = frame.attrs['Time']
    T = np.array(frame['T'])
    return np.c_[ np.array(frame['muB'])[selection], T[selection] ]


#########################################################################################
def plot_all_particles():
    #####################################
    # Set up
    #####################################
    global selection
    Tmin, Tmax = 600, 650
    selection = get_selection(Tmin, Tmax)
    data = np.stack([frame_to_array(i, Tmin, Tmax) for i in range(n_timesteps)])

    data = np.swapaxes(data, 0, 1)
    
    # initial temperatures
    #T0 = event[event_keys[0]]['T']
    #T0 = T0[np.where( (T0 >= Tmin) & (T0 <= Tmax) )]
    
    minimum, maximum = Tmin, Tmax
    
    print(data.shape)
    
    #####################################
    # T vs. muB figure
    #####################################
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    for particle in data:
        r,g,b,a = chosen_colormap((particle[0,1]-minimum)/(maximum-minimum))
        plt.plot( particle[:,0], particle[:,1], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,0]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_B$ (MeV)')
    plt.ylabel(r'$T$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'T_vs_muB.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)



#########################################################################################
if __name__== "__main__":
    plot_all_particles()




