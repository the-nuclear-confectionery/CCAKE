import h5py as h5
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap, LogNorm, SymLogNorm
from scipy.interpolate import CloughTocher2DInterpolator
import os, sys, time

fixed_maximum = True

# load arguments
f = h5.File(sys.argv[1], 'r')
event = f['Event']
event_keys = list(event.keys())
n_timesteps = min([len(event.keys()),10])

h = float(sys.argv[2])
knorm = 10.0/(7.0*np.pi*h*h)

#quantity = sys.argv[3]

fig = plt.figure(figsize=(12,12), dpi=125)
ax = fig.add_subplot(111)

xmin, xmax, ymin, ymax = -15, 15, -15, 15

#n = int(30/0.06)+1
n = 101
colormap = plt.cm.inferno

data = None
maximum = None
minimum = None
im = ax.imshow(np.zeros((n,n)))


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
def kernel(q):
    global knorm
    return np.piecewise(q, [q>=1, q<1], \
                        [lambda q: 0.25*knorm*(2.0-q)**3,\
                         lambda q: knorm*(1.0 - 1.5*q**2 + 0.75*q**3)])

#########################################################################################
def kernel1(q):  # q < 1
    global knorm
    return knorm*(1.0 - 1.5*q**2 + 0.75*q**3)

#########################################################################################
def kernel2(q):  # 2 > q >= 1
    global knorm
    return 0.25*knorm*(2.0-q)**3

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
    global data, maximum, minimum, im
    print('Plotting frame', i, flush=True)
    #fig.clear()
    frame = event[event_keys[i]]
    tau = frame.attrs['Time']
    x = np.array(frame['x'])
    y = np.array(frame['y'])
    T = np.array(frame['T'])
    data = np.c_[ x, y, T ]
    
    #print('max T =', np.amax(T))
    
    #print(data.shape)
        
    X, Y = np.meshgrid( np.linspace(xmin, xmax, n), np.linspace(ymin, ymax, n) )
    #print(np.c_[ X.flatten(), Y.flatten() ].shape)
    tic = time.perf_counter()

    #f = np.array([ evaluate_field(point) for point in np.c_[ X.flatten(), Y.flatten() ] ])

    interp = CloughTocher2DInterpolator(list(zip(x, y)), T, fill_value = 0.0)
    f = interp(X, Y)

    toc = time.perf_counter()
    print(f"Generated field grid in {toc - tic:0.4f} seconds", flush=True)

    
    if i==0 or not fixed_maximum:
        maximum = np.amax(np.abs(f))
        minimum = np.amin(f[np.abs(f)>0.0])
        
    #print(f.shape)
    
    extent = xmin, xmax, ymin, ymax
    #print(i)
    #print(f.reshape(n, n))
    tic = time.perf_counter()
    im = ax.imshow(f.reshape(n, n)+1e-10, cmap=colormap,\
                    norm=LogNorm(vmin=minimum+1e-15, vmax=maximum),\
                    interpolation='bicubic', extent=extent)
    toc = time.perf_counter()
    print(f"Generated frame in {toc - tic:0.4f} seconds")

    #im.set_data(f.reshape(n, n)+1e-15)

    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    
    if i==0:
        fig.savefig('frame' + str(i) + '.png', format='png')
    
    return [im]



#########################################################################################
def main():

    plt.margins(0, 0)
    
    ani = animation.FuncAnimation(fig, animate, np.arange(n_timesteps), \
                                  init_func=init, blit=True)
    #ani = animation.FuncAnimation(fig, animate, frames=10)

    out = "T_evo.gif"
    #out = sys.argv[2]
    print('Saving to', out)
    ani.save(out, writer='imagemagick', fps=25)
    print('Finished everything.')

    return 0

  
#########################################################################################
if __name__== "__main__":
  main()

