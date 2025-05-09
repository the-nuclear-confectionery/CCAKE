import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

use_log_scale = True

infilename = sys.argv[1]
outdirectory = sys.argv[2]
axisMode = sys.argv[3]
eta0 = float(sys.argv[4]) # slice in eta to compare


# HDF file structure
f = h5.File(infilename, 'r')
event = f['Event']
n_timesteps = min([len(event.keys()),1000])
event_keys = list(event.keys())

hbarc = 197.33 # MeV*fm
t0 = 0.5 # amount of temporal shift
q  = 1.
e0 = 80000.0 # units of MeV/fm^3 (cf. ideal_2+1d_gubser_ic_generator.py script)

Nc = 3.
Nf = 2.5
cp = (2.*(Nc**2-1.) + 3.5*Nc*Nf)*np.pi**2/90.

quantities = ['e','ur','ueta']
cols = dict(zip(quantities,range(3,len(quantities)+3)))
print('cols=',cols)


#===============================================================================
def taup(tau, eta):
    return np.sqrt(tau**2 + 2.*t0*tau*np.cosh(eta)+t0**2)
#==============================================================================
def etap(tau, eta):
    return np.arctanh( tau * np.sinh(eta) / (tau * np.cosh(eta) + t0) )
#==============================================================================
def eGubser(tau, r):
    return (e0/tau**4)*( (2.*q*tau)**(8./3.)
                        / ( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )**(4./3.)
                         )
#===============================================================================
def TGubser(tau, r):
    return hbarc*( eGubser(tau, r) / (3.*cp*hbarc) )**0.25
#===============================================================================
def urGubser(tau, r):
    return 2.0*q**2*r*tau / np.sqrt( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )
#===============================================================================
def utauGubser(tau, r):
    return np.sqrt(1.0 + urGubser(tau, r)**2)
#==============================================================================
def shifted_eGubser(tau, r, eta):
    return eGubser(taup(tau, eta), r)
#==============================================================================
def shifted_urGubser(tau, x, r, eta):
    return urGubser(taup(tau, eta), r)
#==============================================================================
def shifted_uetaGubser(tau, r, eta):
    return -utauGubser(taup(tau, eta), r) * ( t0 * np.sinh(eta) / (tau * taup(tau, eta)) )
#==============================================================================
#===============================================================================
def plot_slice(ax, hydroOutput, tau, axis, quantity):
    # c : column of quantity to plot in array
    c = cols[quantity]
    print('quantity =',quantity)
    print('c =',c)
    cf   = [None, None, None, shifted_eGubser, shifted_urGubser, shifted_uetaGubser][c]
    sliceData = hydroOutput[np.where( (np.isclose(hydroOutput[:,1], 0.0)) \
                                      & (np.isclose(hydroOutput[:,2], eta0)) )] # y == 0 ===>>> r == x
    ax.plot( sliceData[:,0], sliceData[:,c], 'r-' )
    xpts = np.linspace(np.amin(sliceData[:,0]), np.amax(sliceData[:,0]), 1001)
    ax.plot( xpts, cf(tau, xpts), 'b:' )
    

#===============================================================================
if __name__ == "__main__":

    # set up figure
    toPlot = ['e', 'ur', 'ueta']
    
    ncols = 3
    nrows = 1
    fig, axs = plt.subplots( ncols=ncols, nrows=nrows, figsize=(5*ncols, 5*nrows) )
    
    tau = 0.0
    
    # plot hydro output files
    for timestep in event_keys:
        # load Gubser check output files produced by hydro code
        # (eventually) use format: x [fm], y [fm], e [1/fm^4], u_x, u_y, ...
        frame = event[timestep]
        tau = frame.attrs['Time'][0]
        print('tau =', tau)

        print(frame.keys())
        #exit(1)

        x = np.array(frame['x'])
        y = np.array(frame['y'])
        eta = np.array(frame['eta'])
        e = np.array(frame['e'])
        ux = np.array(frame['ux'])
        uy = np.array(frame['uy'])
        ueta = np.array(frame['ueta'])

        #exit(1)
        #hydroOutput = np.loadtxt( checkfile, skiprows=1 )
        hydroOutput = np.c_[ x, y, eta, e, ux, ueta ] # ux == ur at y == 0

        # plot comparison along y==0 slice
        for i, ax in enumerate(axs.ravel()):
            if use_log_scale and ['T','e'].count(toPlot[i]) > 0:
                ax.set_yscale('log')
            plot_slice( ax, hydroOutput, tau, axisMode, toPlot[i] )
            ax.set_xlim([-4.5, 4.5])
            ax.set_xlabel(r'$r$ (fm)')
            if toPlot[i] == 'ur':
                ax.set_ylim([-3.0, 3.0])
            #ax.ylabel(r'$e$ (fm$^{-4}$)')
            
        #plt.show()
        print("tau = ", tau)
        print("float(int(tau)) = ", np.round(tau))
        print("np.isclose( tau, float(int(tau)), atol=1e-04 ) = ", np.isclose( tau, float(int(tau)), atol=1e-04 ))
        if np.isclose( tau, np.round(tau), atol=1e-04 ):
            plt.savefig('./yeq' + axisMode + '_slice_tau=' + f"{tau:.2f}" + '.pdf')
    
