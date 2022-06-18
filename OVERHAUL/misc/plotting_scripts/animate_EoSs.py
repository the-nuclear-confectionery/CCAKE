import h5py
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import sys

f = h5py.File(sys.argv[1], 'r')
event = f['Event']
n_timesteps = min([len(event.keys()),5])
event_keys = list(event.keys())


def init():
    frame = event[event_keys[0]]
    x = np.array(frame['x'])
    y = np.array(frame['y'])
    #im = plt.scatter(x, y)
    im = plt.plot(x, y, 'bo', ms = 0.01)
    #return im,
    return im


def animate(i):
    print('Plotting frame', i, flush=True)
    frame = event[event_keys[i]]
    x = np.array(frame['x'])
    y = np.array(frame['y'])
    #im = plt.scatter(x, y, c = np.array(frame['e']), s = 0.000004,
    #                 cmap = cm.get_cmap('plasma') )
    im = plt.plot(x, y, 'bo', ms = 0.01)
    
    #plt.imsave(fname='old_animation_frames/frame' + str(i) + '.png', \
    #           arr=image, cmap=chosen_colormap, format='png')
    #return im,
    return im



def main():

    # Plot Volume Rendering
    fig = plt.figure(figsize=(10,10), dpi=150)
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
        
    # Do Volume Rendering at Different Viewing Angles
    ani = animation.FuncAnimation(fig, animate, np.arange(n_timesteps), \
                                  init_func=init, blit=True)

    out = "EoS_particle_evolution.gif"
    print('Saving to', out)
    #FFwriter = animation.FFMpegWriter(fps=2, extra_args=['-vcodec', 'libx264'])
    #ani.save(out, writer=FFwriter)
    ani.save(out, writer='imagemagick', fps=10)
    print('Finished everything.')

    return 0

  
if __name__== "__main__":
  main()

