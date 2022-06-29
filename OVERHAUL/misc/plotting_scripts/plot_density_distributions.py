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
event_keys = list(event.keys())
n_timesteps = len(event_keys)

# Figure dimensions
width, height = 6, 6

# Axes ranges
xmin, xmax, ymin, ymax = 2.5, 7.5, -7.5, -2.5
#xmin, xmax, ymin, ymax = -12, 12, -12, 12

chosen_colormap = cm.get_cmap('cool', 256)

chosen_dpi = 200

labels = {'e': r'$e$ (MeV/fm$^3$)', 'B': r'$\rho_B$ (fm$^{-3}$)',
          'S': r'$\rho_S$ (fm$^{-3}$)', 'Q': r'$\rho_Q$ (fm$^{-3}$)'}

#########################################################################################
def space_log_ticks(xmin, xmax):
    minTick = int(np.floor(np.log10(xmin)))
    maxTick = int(np.ceil(np.log10(xmax)))
    ticks = np.arange(10**minTick,10**(minTick+1),10**minTick)
    ticks = ticks[ticks >= xmin]
    for i in range(minTick+1,maxTick):
        ticks = np.concatenate((ticks, np.arange(10**i,10**(i+1),10**i)))
    ticks = ticks[ticks <= xmax]
    return ticks

#########################################################################################
def frame_to_array(i, quantity):
    frame = event[event_keys[i]]
    q = np.array(frame[quantity])
    return np.c_[ np.full_like(q, frame.attrs['Time']), q, np.array(frame['e']) ]


#########################################################################################
def plot_density_distribution_vs_time(quantity):
    # Set up
    print('Building data...')
    data = np.stack([frame_to_array(i, quantity) for i in range(n_timesteps)])
    
    #print(data.shape)
    
    data = data.reshape([data.size//3,3])
    
    #print(np.amin(data[:,0]), np.amax(data[:,0]))
    
    ti, tf = np.amin(data[:,0]), np.amax(data[:,0])
    dt = (tf - ti) / (n_timesteps - 1.0)
    #print(ti, tf, tf - ti, n_timesteps, dt)
    timebins = np.arange(ti-dt, tf+dt, dt)
    #timebins = np.arange(0.55,13.15,0.05)
    timebins = 0.5*(timebins[1:]+timebins[:-1])
    
    # freeze-out cutoff [MeV/fm^3]
    eFO = 266.0
    data = data[ data[:,2] >= eFO ]
    
    # only look at positive densities, for simplicity
    data = data[ data[:,1] >= 1e-3 ]
    
    H, xedges, yedges = np.histogram2d(data[:,0], np.log(data[:,1]), bins=[timebins,250])
        
    #####################################
    # quantity vs. tau figure
    #####################################
    print('Plotting...')
    #fig, ax = plt.figure(figsize=(width, height), dpi=chosen_dpi)
    fig, ax = plt.subplots()
    fig.set_size_inches(width, height, forward=True)
    fig.set_dpi(chosen_dpi)
    
    #print(xedges.shape, yedges.shape, H.shape)

    #plt.pcolormesh(yedges, xedges, H.T, cmap='inferno')
    im = plt.imshow(H.T, cmap='inferno', interpolation='bicubic', origin='lower',\
               extent=[np.amin(xedges),np.amax(xedges),\
                       np.amin(yedges),np.amax(yedges)])
    
    # set y-axis ticks in an aesthetic way
    ymin, ymax = np.amin(yedges), np.amax(yedges)
    chosen_y_tick_values = space_log_ticks(np.exp(ymin), np.exp(ymax))
    chosen_y_tick_labels \
        = np.array(list(map(lambda y: \
                   r'$10^{y}$'.replace('y',str(int(np.log10(y)))) \
                   if np.isclose(np.log10(y), np.round(np.log10(y))) else '',\
                   chosen_y_tick_values)))
    
    isMajor = lambda y: np.isclose(np.log10(y), np.round(np.log10(y))) 
    #ax.set_yticks(np.log(chosen_y_tick_values))  # take log since y axis is log
    chosen_major_yticks = chosen_y_tick_values[isMajor(chosen_y_tick_values)]
    chosen_minor_yticks = chosen_y_tick_values[~isMajor(chosen_y_tick_values)]
    ax.set_yticks(np.log(chosen_major_yticks))
    ax.set_yticks(np.log(chosen_minor_yticks), minor=True)
    
    ax.set_yticklabels(chosen_y_tick_labels[isMajor(chosen_y_tick_values)])
    #ax.set_yticklabels(list(map(str,np.exp(yedges)[::50])))
    
    plt.xlabel(r'$\tau$ (fm/$c$)')
    plt.ylabel(labels[quantity])
    
    div = make_axes_locatable(ax)
    cax = div.append_axes('right', '5%', '5%')
    cbar = fig.colorbar(im, cax=cax, extend='both')
    cbar.set_label("Number of cells", fontsize=16)
    cbar.ax.tick_params(labelsize=14)
    
    outfilename = quantity + '_vs_tau.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)



#########################################################################################
if __name__== "__main__":
    plot_density_distribution_vs_time('e')
    plot_density_distribution_vs_time('B')
    plot_density_distribution_vs_time('S')
    plot_density_distribution_vs_time('Q')



