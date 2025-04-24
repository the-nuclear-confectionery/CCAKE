import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

use_log_scale = False

infilename = sys.argv[1]
outdirectory = sys.argv[2]
axisMode = sys.argv[3]

# HDF file structure
f = h5.File(infilename, 'r')
event = f['Event']
n_timesteps = min([len(event.keys()),1000])
event_keys = list(event.keys())

q  = 1.
e0 = 1.0

Nc = 3.
Nf = 2.5
cp = (2.*(Nc**2-1.) + 3.5*Nc*Nf)*np.pi**2/90.

#quantities = ['T','e','ux','uy']
quantities = ['T','e']
cols = dict(zip(quantities,range(2,len(quantities)+2)))
print('cols=',cols)


#===============================================================================
def eGubser(tau, r):
    return (e0/tau**4)*( (2.*q*tau)**(8./3.)
                        / ( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )**(4./3.)
                         )

#===============================================================================
def TGubser(tau, r):
    return ( eGubser(tau, r) / (3.*cp) )**0.25

#===============================================================================
def eFromT(T):
    return 3.*cp*T**4

#===============================================================================
def urGubser(tau, r):
    return 2.0*q**2*r*tau / np.sqrt( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )

#===============================================================================
#def get_selection(Tmin, Tmax):
#    frame    = event[event_keys[0]]
#    return np.where( (T0 >= Tmin) & (T0 <= Tmax) & (eos_tags == 0) )


#===============================================================================
#def frame_to_array(i, Tmin, Tmax):
#    frame = event[event_keys[i]]
#    tau = frame.attrs['Time']
#    T = np.array(frame['T'])
#    return np.c_[ T[selection], np.array(frame['muB'])[selection],\
#                  np.array(frame['muS'])[selection], np.array(frame['muQ'])[selection] ]


#===============================================================================
#def plot_all_particles():
#    # Set up
#    print('Setting up...', flush=True)
#    global selection
#    Tmin, Tmax = 400, 650
#    selection = get_selection(Tmin, Tmax)
#    data = np.stack([frame_to_array(i, Tmin, Tmax) for i in range(n_timesteps)])




#===============================================================================
def plot_slice(ax, hydroOutput, tau, axis, quantity):
    # c : column of quantity to plot in array
    c = cols[quantity]
    print('quantity=',quantity)
    print('c=',c)
    #cf   = [TGubser, eGubser, urGubser, urGubser][c]
    cf   = [None, None, TGubser, eGubser][c]
    if axis == '0':
        yEqAxisData = hydroOutput[np.where( np.abs(hydroOutput[:,1]) < 1e-6 )]
        ax.plot( yEqAxisData[:,0], yEqAxisData[:,c], 'r-' )
        xpts = np.linspace(np.amin(yEqAxisData[:,0]), np.amax(yEqAxisData[:,0]), 1001)
        ax.plot( xpts, cf(tau, xpts), 'b:' )
    elif axis == 'x':
        yeqxData = hydroOutput[np.where( np.isclose( hydroOutput[:,0], hydroOutput[:,1] ) )]
        rpts = np.sqrt(yeqxData[:,0]**2 + yeqxData[:,1]**2)
        ax.plot( rpts, yeqxData[:,c], 'r-' )
        rpts = np.linspace(0.0, np.amax(rpts), 1001)
        ax.plot( rpts, cf(tau, rpts), 'b:' )
    

#===============================================================================
if __name__ == "__main__":

    # set up figure
    toPlot = ['T', 'e']
    
    ncols = len(toPlot)
    nrows = 1
    fig, axs = plt.subplots( ncols=ncols, nrows=nrows, figsize=(5*ncols, 5*nrows) )

    # plot hydro output files
    for timestep in event_keys:
        # load Gubser check output files produced by hydro code
        # (eventually) use format: x [fm], y [fm], e [1/fm^4], u_x, u_y, ...
        frame = event[timestep]
        tau = frame.attrs['Time']
        print('tau=', tau)

        print(frame.keys())
        #exit(1)

        x = np.array(frame['x'])
        y = np.array(frame['y'])
        T = np.array(frame['T'])
        e = np.array(frame['e'])
        #ux = np.array(frame['ux'])
        #uy = np.array(frame['uy'])

        #exit(1)
        #hydroOutput = np.loadtxt( checkfile, skiprows=1 )
        hydroOutput = np.c_[ x, y, T, e ]

        # plot comparison along y==0 slice
        for i, ax in enumerate(axs.ravel()):
            if use_log_scale and ['T','e'].count(toPlot[i]) > 0:
                ax.set_yscale('log')
            plot_slice( ax, hydroOutput, tau, axisMode, toPlot[i] )
            if axisMode == '0':
                ax.set_xlim([-4.5, 4.5])
                ax.set_xlabel(r'$x$ (fm)')
            else:
                ax.set_xlim([0.0, 4.5])
                ax.set_xlabel(r'$r$ (fm)')
            if toPlot[i] == 'ux' or toPlot[i] == 'uy':
                ax.set_ylim([-2.0, 2.0])
            #ax.ylabel(r'$e$ (fm$^{-4}$)')
            
    #plt.show()
    plt.savefig('./yeq' + axisMode + '_slice_tau=' + tau + '.pdf')
    
