import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm                                                                          
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import numpy as np
import sys

chosen_colormap = cm.get_cmap('cool', 256)

#########################################################################################
def plot_one_particle(pid):
    data = np.stack([np.loadtxt(file, usecols=(1,5,6,7,8), skiprows=1) \
                     for file in sys.argv[1:]])

    data = np.swapaxes(data, 0, 1)

    plt.figure(figsize=(6,6), dpi=100)

    plt.plot( data[pid,:,0], data[pid,:,1] )
    plt.plot( data[pid,:,0], data[pid,:,2] )
    plt.plot( data[pid,:,0], data[pid,:,3] )
    plt.plot( data[pid,:,0], data[pid,:,4] )

    plt.savefig('particle_vs_t.png', dpi=100, bbox_inches='tight', pad_inches = 0)


#########################################################################################
def plot_all_particles():
    data = np.stack([np.loadtxt(file, usecols=(5,6,7,8), skiprows=1) \
                     for file in sys.argv[1:]])
    
    chosen_dpi = 200

    print(data.shape)

    data0 = data[0,:,0]
    data = np.swapaxes(data, 0, 1)[np.where(data0 > 400.0)]
    data0 = data0[np.where(data0 > 400.0)]
    
    maximum, minimum = np.amax(data0), np.amin(data0)
    
    print("minimum =", minimum)
    print("maximum =", maximum)

    print(data.shape)

    plt.figure(figsize=(4,4), dpi=chosen_dpi)

    for iParticle, particle in enumerate(data):
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,1], particle[:,0], color=(r,g,b), alpha=0.3 )

    plt.savefig('T_vs_muB.png', dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to T_vs_muB.png')

    plt.figure(figsize=(4,4), dpi=chosen_dpi)

    for iParticle, particle in enumerate(data):
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,2], particle[:,0], color=(r,g,b), alpha=0.3 )

    plt.savefig('T_vs_muS.png', dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to T_vs_muS.png')

    plt.figure(figsize=(4,4), dpi=chosen_dpi)

    for iParticle, particle in enumerate(data):
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,3], particle[:,0], color=(r,g,b), alpha=0.3 )

    plt.savefig('T_vs_muQ.png', dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to T_vs_muQ.png')


#########################################################################################
if __name__== "__main__":
    plot_all_particles()
    #plot_one_particle(7108)

