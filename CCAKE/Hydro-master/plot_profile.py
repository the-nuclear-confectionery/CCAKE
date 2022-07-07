import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os, sys


filenameStem = sys.argv[1]
outfilenameStem = sys.argv[2]
insuffix = sys.argv[3]
outsuffix = sys.argv[4]
numberOfFrames = int(sys.argv[5])
plotLabel = sys.argv[6]
colToPlot = int(sys.argv[7])


for i in range(1, numberOfFrames+1):
    outfilename = outfilenameStem + f'{i:03}' + outsuffix
    if os.path.isfile(outfilename): # if file already exists, don't bother creating it again
        continue

    filename = filenameStem + str(i) + insuffix
    data = np.loadtxt(filename, usecols=(0,1,2,colToPlot))
    tau = data[0,0]
    data = data[:,1:]  # set tau and then discard that column
    data = data[np.where((np.abs(data[:,0])<5.0) & (np.abs(data[:,1])<5.0))]
    [x, y, f] = data.T

    extent = np.min(x), np.max(x), np.min(y), np.max(y)
    fig, ax = plt.subplots( nrows=1, ncols=1 )

    length = int(np.sqrt(f.size))
    psm = plt.imshow(f.reshape(length, length), cmap=plt.cm.viridis, interpolation='bicubic', extent=extent)

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






########################################################################
### OLD VERSION BELOW THIS LINE
########################################################################

'''
filename = sys.argv[1]
outfilename = sys.argv[2]
plotLabel = sys.argv[3]
colToPlot = int(sys.argv[4])

data = np.loadtxt(filename, usecols=(0,1,2,colToPlot))
tau = data[0,0]
data = data[:,1:]  # set tau and then discard that column
data = data[np.where((np.abs(data[:,0])<5.0) & (np.abs(data[:,1])<5.0))]

[x, y, f] = data.T

extent = np.min(x), np.max(x), np.min(y), np.max(y)
fig, ax = plt.subplots( nrows=1, ncols=1 )

length = int(np.sqrt(f.size))
psm = plt.imshow(f.reshape(length, length), cmap=plt.cm.viridis, interpolation='bicubic', extent=extent)
#ax.set_facecolor('black')

plt.text(0.075, 0.925, r'$\tau = %(t)5.2f$ fm$/c$'%{'t': tau}, \
        {'color': 'white', 'fontsize': 12}, transform=ax.transAxes,
        horizontalalignment='left', verticalalignment='top')


ax.set_xlabel(r'$x$ (fm)', fontsize=16)
ax.set_ylabel(r'$y$ (fm)', fontsize=16)
cbar = fig.colorbar(psm, ax=ax)
cbar.set_label(plotLabel, fontsize=16)


#plt.show()
#dirname = os.path.dirname(filename)
#print('Saving to', outfilename)
fig.savefig(outfilename, bbox_inches='tight')
plt.close(fig)
'''

