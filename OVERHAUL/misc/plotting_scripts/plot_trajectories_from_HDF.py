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

#########################################################################################
def get_selection(Tmin, Tmax):
    frame    = event[event_keys[0]]
    T0       = np.array(frame['T'])
    eos_tags = np.array(frame['labels'])
    #print(np.amin(T0),np.amax(T0))
    return np.where( (T0 >= Tmin) & (T0 <= Tmax) & (eos_tags == 0) )


#########################################################################################
def frame_to_array(i, Tmin, Tmax):
    frame = event[event_keys[i]]
    tau = frame.attrs['Time']
    T = np.array(frame['T'])
    return np.c_[ T[selection], np.array(frame['muB'])[selection],\
                  np.array(frame['muS'])[selection], np.array(frame['muQ'])[selection] ]


#########################################################################################
def plot_all_particles():
    #####################################
    # Set up
    #####################################
    global selection
    Tmin, Tmax = 400, 650
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
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,1], particle[:,0], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,1]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_B$ (MeV)')
    plt.ylabel(r'$T$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'T_vs_muB_Tmin'+str(Tmin)+'_Tmax'+str(Tmax)+'.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)

    #####################################
    # T vs. muS figure
    #####################################
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    for particle in data:
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,2], particle[:,0], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,2]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_S$ (MeV)')
    plt.ylabel(r'$T$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'T_vs_muS_Tmin'+str(Tmin)+'_Tmax'+str(Tmax)+'.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)

    #####################################
    # T vs. muQ figure
    #####################################
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    for particle in data:
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,3], particle[:,0], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,3]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_Q$ (MeV)')
    plt.ylabel(r'$T$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'T_vs_muQ_Tmin'+str(Tmin)+'_Tmax'+str(Tmax)+'.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)


    #####################################
    # muB vs. muS figure
    #####################################
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    for particle in data:
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,1], particle[:,2], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,2]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_S$ (MeV)')
    plt.ylabel(r'$\mu_B$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'muB_vs_muS_Tmin'+str(Tmin)+'_Tmax'+str(Tmax)+'.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)


    #####################################
    # muB vs. muQ figure
    #####################################
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    for particle in data:
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,1], particle[:,3], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,3]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_Q$ (MeV)')
    plt.ylabel(r'$\mu_B$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'muB_vs_muQ_Tmin'+str(Tmin)+'_Tmax'+str(Tmax)+'.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)

    #####################################
    # muS vs. muQ figure
    #####################################
    plt.figure(figsize=(width, height), dpi=chosen_dpi)

    for particle in data:
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,2], particle[:,3], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,3]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_Q$ (MeV)')
    plt.ylabel(r'$\mu_S$ (MeV)')
    plt.title(r'Pb+Pb $5.02$ TeV')
    norm = Normalize(vmin=minimum, vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'$T_0$ (MeV)', rotation=90)

    outfilename = 'muS_vs_muQ_Tmin'+str(Tmin)+'_Tmax'+str(Tmax)+'.png'
    plt.savefig(outfilename, dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to ' + outfilename)



#########################################################################################
if __name__== "__main__":
    plot_all_particles()




