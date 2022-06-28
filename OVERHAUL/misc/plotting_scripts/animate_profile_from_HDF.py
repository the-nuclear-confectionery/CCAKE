import h5py as h5
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap, LogNorm, SymLogNorm
import os, sys, time

#########################################################################################
# load arguments
infilename = sys.argv[1]
h = float(sys.argv[2])
quantity = sys.argv[3]
outfilename = sys.argv[4]

f = h5.File(infilename, 'r')
event = f['Event']
event_keys = list(event.keys())
n_timesteps = min([len(event_keys),10])

fig = plt.figure(figsize=(12,12), dpi=125)
ax = fig.add_subplot(111)
xmin, xmax, ymin, ymax = -15, 15, -15, 15
n = 51

knorm = 10.0/(7.0*np.pi*h*h)
fixed_maximum = True

#########################################################################################
cmap_energy \
    = LinearSegmentedColormap.from_list( 'energy',\
        ( (0.000, (1.000, 1.000, 1.000)),\
          (0.500, (0.6392156862745098, 0.45098039215686275, 0.3176470588235294)),\
          (1.000, (0.23137254901960785, 0.1607843137254902, 0.11372549019607843)) ) ) 

#########################################################################################
cmap_baryon \
    = LinearSegmentedColormap.from_list( 'baryon',\
        ( (0.000, (0.0, 0.5764705882352941, 0.5764705882352941)),\
          (0.500, (1.000, 1.000, 1.000)),\
          (1.000, (0.9568627450980393, 0.47843137254901963, 0.000)) ) )

#########################################################################################
cmap_strange \
    = LinearSegmentedColormap.from_list( 'strange',\
        ( (0.000, (0.5019607843137255, 0.0, 1.0)),\
          (0.500, (1.000, 1.000, 1.000)),\
          (1.000, (0.0, 0.5019607843137255, 0.000)) ) )

#########################################################################################
cmap_electric \
    = LinearSegmentedColormap.from_list( 'electric',\
        ( (0.000, (0.0, 0.0, 1.0)),\
          (0.500, (1.000, 1.000, 1.000)),\
          (1.000, (1.0, 0.0, 0.000)) ) )

#########################################################################################
colormap                                              \
    = {"energy_density":              cmap_energy,    \
       "baryon_density":              cmap_baryon,    \
       "strange_density":             cmap_strange,   \
       "electric_density":            cmap_electric,  \
       "temperature":                 plt.cm.inferno, \
       "baryon_chemical_potential":   cmap_baryon,    \
       "strange_chemical_potential":  cmap_strange,   \
       "electric_chemical_potential": cmap_electric}[quantity]

#########################################################################################
quantityLabel                                 \
    = {"energy_density":              'e',    \
       "baryon_density":              'rhoB', \
       "strange_density":             'rhoS', \
       "electric_density":            'rhoQ', \
       "temperature":                 'T',    \
       "baryon_chemical_potential":   'muB',  \
       "strange_chemical_potential":  'muS',  \
       "electric_chemical_potential": 'muQ'}[quantity]
    
#########################################################################################
plotLabel                                                      \
    = {"energy_density":              r'$e$ (MeV/fm$^3$)',     \
       "baryon_density":              r'$\rho_B$ (fm$^{-3}$)', \
       "strange_density":             r'$\rho_S$ (fm$^{-3}$)', \
       "electric_density":            r'$\rho_Q$ (fm$^{-3}$)', \
       "temperature":                 r'$T$ (MeV)',            \
       "baryon_chemical_potential":   r'$\mu_B$ (MeV)',        \
       "strange_chemical_potential":  r'$\mu_S$ (MeV)',        \
       "electric_chemical_potential": r'$\mu_Q$ (MeV)'}[quantity]

#########################################################################################
use_log_scale                                 \
    = {"energy_density":              True,   \
       "baryon_density":              True,   \
       "strange_density":             True,   \
       "electric_density":            True,   \
       "temperature":                 False,  \
       "baryon_chemical_potential":   False,  \
       "strange_chemical_potential":  False,  \
       "electric_chemical_potential": False}[quantity]


#########################################################################################
data = None
maximum = None
minimum = None
im = ax.imshow(np.zeros((n,n)))


#########################################################################################
def kernel(q):
    global knorm
    return np.piecewise(q, [q>=1, q<1], \
                        [lambda q: 0.25*knorm*(2.0-q)**3,\
                         lambda q: knorm*(1.0 - 1.5*q**2 + 0.75*q**3)])

#########################################################################################
def evaluate_field(r):
    global h
    neighbors = data[ (r[0]-data[:,0])**2+(r[1]-data[:,1])**2 <= 4.0*h**2 ]
    weights = kernel( np.sqrt( (r[0]-neighbors[:,0])**2+(r[1]-neighbors[:,1])**2 )/h )
    if np.sum(weights) < 1e-10:
        return 0
    else:
        return np.sum( neighbors[:,2]*weights ) / (np.sum(weights)+1e-10)

#########################################################################################
def init():
    im.set_data(np.zeros((n,n)))
    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    return [im]


#########################################################################################
def animate(i):
    global data, maximum, minimum, im, n
    print('Plotting frame', i, flush=True)
    tic = time.perf_counter()

    fig.clear()
    frame = event[event_keys[i]]
    tau = frame.attrs['Time']
    x = np.array(frame['x'])
    y = np.array(frame['y'])
    T = np.array(frame['T'])
    data = np.c_[ x, y, T ]
    
    toc = time.perf_counter()
    print(f"Set up frame in {toc - tic:0.4f} seconds", flush=True)

    X, Y = np.meshgrid( np.linspace(xmin, xmax, n), np.linspace(ymin, ymax, n) )
    tic = time.perf_counter()

    f = np.array([ evaluate_field(point) for point in np.c_[ X.flatten(), Y.flatten() ] ])

    toc = time.perf_counter()
    print(f"Generated field grid in {toc - tic:0.4f} seconds", flush=True)

    
    if i==0 or not fixed_maximum:
        maximum = np.amax(np.abs(f))
        #minimum = np.amin(f[np.abs(f)>0.0])
        minimum = {"energy_density": np.amin(f[np.abs(f)>0.0]), \
                   "baryon_density": -maximum, \
                   "strange_density": -maximum, \
                   "electric_density": -maximum, \
                   "temperature": np.amin(f[np.abs(f)>0.0]), \
                   "baryon_chemical_potential": -maximum, \
                   "strange_chemical_potential": -maximum, \
                   "electric_chemical_potential": -maximum}[quantity]
    
    
    extent = xmin, xmax, ymin, ymax
    tic = time.perf_counter()
        
    #im = ax.imshow(f.reshape(n, n)+1e-10, cmap=colormap,\
    #                norm=LogNorm(vmin=minimum+1e-15, vmax=maximum),\
    #                interpolation='bicubic', extent=extent)
    if use_log_scale:
        if quantity == "energy_density" or quantity == "temperature":
            im = ax.imshow(f.reshape(n, n)+1e-10, cmap=colormap,\
                           norm=LogNorm(vmin=minimum+1e-15, vmax=maximum),\
                           interpolation='bicubic', extent=extent)
        else:
            # set range around zero to be linear instead of logarithmic
            linthresh = 0.01*np.amax(np.abs(f))
            im = ax.imshow(f.reshape(n, n), cmap=colormap,\
                           norm=SymLogNorm(linthresh, vmin=minimum, vmax=maximum),\
                           interpolation='bicubic', extent=extent)
    else:
        im = ax.imshow(f.reshape(n, n), cmap=colormap,\
                       vmin=minimum, vmax=maximum,\
                       interpolation='bicubic', extent=extent)
    
    toc = time.perf_counter()
    print(f"Generated frame in {toc - tic:0.4f} seconds")

    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.text(0.075, 0.925, r'$\tau = %(t)5.2f$ fm$/c$'%{'t': tau}, \
            {'color': 'white', 'fontsize': 24}, transform=ax.transAxes,
            horizontalalignment='left', verticalalignment='top')
    ax.set_xlabel(r'$x$ (fm)', fontsize=16)
    ax.set_ylabel(r'$y$ (fm)', fontsize=16)
    cbar = fig.colorbar(im, ax=ax)
    cbar.set_label(plotLabel, fontsize=16)

    return [im]



#########################################################################################
def main():
    global outfilename
    plt.margins(0, 0)
    
    ani = animation.FuncAnimation(fig, animate, np.arange(n_timesteps), \
                                  init_func=init, blit=True)
    
    FFwriter = animation.FFMpegWriter(fps=25, extra_args=['-vcodec', 'libx264'])
    ani.save(outfilename, writer=FFwriter)

    print('Finished everything.')

    return 0

  
#########################################################################################
if __name__== "__main__":
    main()

