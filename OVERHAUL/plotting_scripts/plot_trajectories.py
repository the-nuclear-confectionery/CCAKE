import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import cm                                                                          
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import numpy as np
import sys
from mpl_toolkits.axes_grid1 import make_axes_locatable

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
def plot_all_particles(Tmin, Tmax):
    data = np.stack([np.loadtxt(file, usecols=(5,6,7,8,30), skiprows=1) \
                     for file in sys.argv[1:]])
    
    chosen_dpi = 200

    print(data.shape)

    data0 = data[0,:,0]       # Temperature values at initial timestep
    dataFreeze = data[-1,:,4] # Freeze flag column at final timestep
    particleSelectionCriteria = (data0 <= float(Tmax))&(data0 >= float(Tmin))&(dataFreeze != 5)
    data = np.swapaxes(data, 0, 1)[np.where(particleSelectionCriteria)]
    #data0 = data0[np.where(particleSelectionCriteria)]
    
    maximum, minimum = np.amax(data0), np.amin(data0)
    
    print("minimum =", minimum)
    print("maximum =", maximum)

    data0 = data0[np.where(particleSelectionCriteria)]

    print(data.shape)

    plt.figure(figsize=(4,4), dpi=chosen_dpi)

    for iParticle, particle in enumerate(data):
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,1], particle[:,0], color=(r,g,b), alpha=0.3 )
        
    maxrange = np.amax(np.abs(data[:,:,1]))
    plt.xlim([-1.1*maxrange, 1.1*maxrange])
    plt.xlabel(r'$\mu_B$ (MeV)')
    plt.ylabel(r'$T$ (MeV)')
    norm = Normalize(vmin=minimum,vmax=maximum)
    sm = cm.ScalarMappable(cmap=chosen_colormap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm)    
    cbar.set_label(r'Initial $T$ (MeV)', rotation=90, fontsize=14)

    plt.savefig('T_'+str(Tmin)+'_to_'+str(Tmax)+'_vs_muB.png', \
                dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to T_'+str(Tmin)+'_to_'+str(Tmax)+'_vs_muB.png')

    plt.figure(figsize=(4,4), dpi=chosen_dpi)

    for iParticle, particle in enumerate(data):
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,2], particle[:,0], color=(r,g,b), alpha=0.3 )

    plt.savefig('T_'+str(Tmin)+'_to_'+str(Tmax)+'_vs_muS.png', \
                dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to T_'+str(Tmin)+'_to_'+str(Tmax)+'_vs_muS.png')

    plt.figure(figsize=(4,4), dpi=chosen_dpi)

    for iParticle, particle in enumerate(data):
        r,g,b,a = chosen_colormap((particle[0,0]-minimum)/(maximum-minimum))
        plt.plot( particle[:,3], particle[:,0], color=(r,g,b), alpha=0.3 )

    plt.savefig('T_'+str(Tmin)+'_to_'+str(Tmax)+'_vs_muQ.png', \
                dpi=chosen_dpi, bbox_inches='tight', pad_inches = 0)
    print('Saved to T_'+str(Tmin)+'_to_'+str(Tmax)+'_vs_muQ.png')


#########################################################################################
if __name__== "__main__":
    #plot_all_particles(600,700)
    #plot_all_particles(500,600)
    #plot_all_particles(400,500)
    #plot_all_particles(300,400)
    #plot_all_particles(200,300)
    #plot_all_particles(100,200)
    plot_all_particles(600,1000)
    #plot_one_particle(7108)

