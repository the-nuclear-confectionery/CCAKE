import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
import os, sys

#cmap_energy \
#    = LinearSegmentedColormap.from_list( 'energy',\
#        ( (0.000, (1.000, 1.000, 1.000)),\
#          (0.250, (1.000, 1.000, 1.000)),\
#          (0.500, (1.000, 1.000, 1.000)),\
#          (0.750, (0.6392156862745098, 0.45098039215686275, 0.3176470588235294)),\
#          (1.000, (0.23137254901960785, 0.1607843137254902, 0.11372549019607843)) ) )

cmap_energy \
    = LinearSegmentedColormap.from_list( 'energy',\
        ( (0.000, (1.000, 1.000, 1.000)),\
          (0.500, (0.6392156862745098, 0.45098039215686275, 0.3176470588235294)),\
          (1.000, (0.23137254901960785, 0.1607843137254902, 0.11372549019607843)) ) ) 


cmap_baryon \
    = LinearSegmentedColormap.from_list( 'baryon',\
        ( (0.000, (0.0, 0.5764705882352941, 0.5764705882352941)),\
          (0.500, (1.000, 1.000, 1.000)),\
          (1.000, (0.9568627450980393, 0.47843137254901963, 0.000)) ) )

cmap_strange \
    = LinearSegmentedColormap.from_list( 'strange',\
        ( (0.000, (0.5019607843137255, 0.0, 1.0)),\
          (0.500, (1.000, 1.000, 1.000)),\
          (1.000, (0.0, 0.5019607843137255, 0.000)) ) )

cmap_electric \
    = LinearSegmentedColormap.from_list( 'electric',\
        ( (0.000, (0.0, 0.0, 1.0)),\
          (0.500, (1.000, 1.000, 1.000)),\
          (1.000, (1.0, 0.0, 0.000)) ) )

filenameStem = sys.argv[1]
outfilenameStem = sys.argv[2]
insuffix = sys.argv[3]
outsuffix = sys.argv[4]
numberOfFrames = int(sys.argv[5])
plotLabel = sys.argv[6]
colToPlot = int(sys.argv[7])
mode = sys.argv[8]
use_log_scale = True if sys.argv[9] == "log" else False
colormap = {"energy_density": cmap_energy, \
            "baryon_density": cmap_baryon, \
            "strange_density": cmap_strange, \
            "electric_density": cmap_electric, \
            "temperature": plt.cm.inferno, \
            "baryon_chemical_potential": cmap_baryon, \
            "strange_chemical_potential": cmap_strange, \
            "electric_chemical_potential": cmap_electric}[mode]

fixed_maximum = True
minimum = 0.0
maximum = 0.0

for i in range(1, numberOfFrames+1):
    outfilename = outfilenameStem + f'{i:03}' + outsuffix
    if os.path.isfile(outfilename): # if file already exists, don't bother creating it again
        continue

    filename = filenameStem + str(i) + insuffix
    data = np.loadtxt(filename, usecols=(0,1,2,colToPlot))
    tau = data[0,0]
    data = data[:,1:]  # set tau and then discard that column
    data = data[np.where((np.abs(data[:,0])<15.0) & (np.abs(data[:,1])<15.0))]
    [x, y, f] = data.T
    
    if i==1 or not fixed_maximum:
        maximum = np.amax(np.abs(f))

    minimum = {"energy_density": 0.0, \
               "baryon_density": -maximum, \
               "strange_density": -maximum, \
               "electric_density": -maximum, \
               "temperature": 0.0, \
               "baryon_chemical_potential": -maximum, \
               "strange_chemical_potential": -maximum, \
               "electric_chemical_potential": -maximum}[mode]

    extent = np.min(x), np.max(x), np.min(y), np.max(y)
    fig, ax = plt.subplots( nrows=1, ncols=1 )

    length = int(np.sqrt(f.size))
    #psm = plt.imshow(f.reshape(length, length), cmap=plt.cm.inferno, interpolation='bicubic', extent=extent)
    if use_log_scale:
        psm = plt.imshow(f.reshape(length, length), cmap=colormap,\
                         norm=LogNorm(vmin=minimum, vmax=maximum),\
                         interpolation='bicubic', extent=extent)
    else:
        psm = plt.imshow(f.reshape(length, length), cmap=colormap,\
                         vmin=minimum, vmax=maximum,\
                         interpolation='bicubic', extent=extent)
        

    plt.text(0.075, 0.925, r'$\tau = %(t)5.2f$ fm$/c$'%{'t': tau}, \
            {'color': 'white', 'fontsize': 12}, transform=ax.transAxes,
            horizontalalignment='left', verticalalignment='top')

    ax.set_xlabel(r'$x$ (fm)', fontsize=16)
    ax.set_ylabel(r'$y$ (fm)', fontsize=16)
    cbar = fig.colorbar(psm, ax=ax)
    cbar.set_label(plotLabel, fontsize=16)

    #plt.show()
    fig.savefig(outfilename, bbox_inches='tight')
    plt.close(fig)
    
    print('Generated', outfilename, 'from', filename, flush=True)
