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
n_timesteps = min([len(event.keys()),1000])
event_keys = list(event.keys())

fig = plt.figure(figsize=(12,12), dpi=125)


def init():
    im = plt.plot(np.array([0.0]), np.array([0.0]), alpha = 0.0)
    plt.xlim([-12, 12])
    plt.ylim([-12, 12])
    return im


def animate(i):
    print('Plotting frame', i, flush=True)
    fig.clear()
    frame = event[event_keys[i]]
    table_EOS              = np.where( np.array(frame['labels']) == 0 )
    tanh_conformal_EOS     = np.where( np.array(frame['labels']) == 1 )
    conformal_EOS          = np.where( np.array(frame['labels']) == 2 )
    conformal_diagonal_EOS = np.where( np.array(frame['labels']) == 3 )
    x = np.array(frame['x'])
    y = np.array(frame['y'])
    #im = plt.scatter(x, y, c = np.array(frame['e']), s = 0.000004,
    #                 cmap = cm.get_cmap('plasma') )
    im = plt.plot(x[table_EOS], y[table_EOS], 'o', color='blue', ms = 1)
    im = plt.plot(x[tanh_conformal_EOS], y[tanh_conformal_EOS], 'o', color='green', ms = 1)
    im = plt.plot(x[conformal_EOS], y[conformal_EOS], 'o', color='purple', ms = 1)
    im = plt.plot(x[conformal_diagonal_EOS], y[conformal_diagonal_EOS], 'o', color='red', ms = 1)
    plt.xlim([-12, 12])
    plt.ylim([-12, 12])
    
    if i==0:
        plt.savefig('frame' + str(i) + '.png', format='png')
    #return im,
    return im



def main():

    # Plot Volume Rendering
    plt.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    plt.margins(0, 0)
        
    # Do Volume Rendering at Different Viewing Angles
    ani = animation.FuncAnimation(fig, animate, np.arange(n_timesteps), \
                                  init_func=init, blit=True)

    #out = "EoS_particle_evolution.gif"
    out = sys.argv[2]
    print('Saving to', out)
    #FFwriter = animation.FFMpegWriter(fps=2, extra_args=['-vcodec', 'libx264'])
    #ani.save(out, writer=FFwriter)
    ani.save(out, writer='imagemagick', fps=10)
    print('Finished everything.')

    return 0

  
if __name__== "__main__":
  main()

