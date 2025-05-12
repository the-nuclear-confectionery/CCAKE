import h5py as h5
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

use_log_scale = False

#infilename = sys.argv[1]
#outdirectory = sys.argv[2]
#axisMode = sys.argv[3]
#eta0 = float(sys.argv[4]) # slice in eta to compare

axisMode = sys.argv[1]
eta0 = float(sys.argv[2]) # slice in eta to compare
outdirectory = sys.argv[3]
infilenames = sys.argv[4:]



# HDF file structure
#f = h5.File(infilename, 'r')
#event = f['Event']
n_timesteps = min([len(infilenames),1000])
#event_keys = list(event.keys())

hbarc = 197.33 # MeV*fm
t0 = 0.5 # amount of temporal shift
q  = 1.
e0 = 80000.0 # units of MeV/fm^3 (cf. ideal_2+1d_gubser_ic_generator.py script)

Nc = 3.
Nf = 2.5
cp = (2.*(Nc**2-1.) + 3.5*Nc*Nf)*np.pi**2/90.

h = 0.1  # 3D cubic spline smoothing scale
knorm  = 1./(np.pi*h**3)

# points at which to compare exact and numerical solutions
xGrid = np.linspace(-5.0, 5.0, 1001)
yGrid = 0.0*xGrid
etaGrid = eta0 + 0.0*xGrid
grid3D = np.c_[ xGrid, yGrid, etaGrid ]

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
def shifted_urGubser(tau, r, eta):
    return urGubser(taup(tau, eta), r)
#==============================================================================
def shifted_uetaGubser(tau, r, eta):
    return -utauGubser(taup(tau, eta), r) * ( t0 * np.sinh(eta) / (tau * taup(tau, eta)) )
#==============================================================================
#===============================================================================
def get_time_step(filename):
    return float((open(filename)).readline()) # first line is header containing just timestep

#===============================================================================
def kernel(q):
    global knorm
    return np.piecewise(q, [q>=1, q<1], \
                        [lambda q: 0.25*knorm*(2.0-q)**3,\
                         lambda q: knorm*(1.0 - 1.5*q**2 + 0.75*q**3)])

#===============================================================================
def evaluate_field(r):
    global h, hydroOutput
    neighbors = hydroOutput[ ( r[0]-hydroOutput[:,0])**2 \
                             +(r[1]-hydroOutput[:,1])**2 \
                             +(r[2]-hydroOutput[:,2])**2 <= 4.0*h**2 ]
    weights = kernel( np.sqrt( ( r[0]-neighbors[:,0])**2 \
                               +(r[1]-neighbors[:,1])**2 \
                               +(r[2]-neighbors[:,2])**2 )/h )
    if np.sum(weights) < 1e-10:
        return 0
    else:
        return np.sum( neighbors*weights[:, np.newaxis], axis=0 ) / (np.sum(weights)+1e-10)

#===============================================================================
#def plot_slice(ax, hydroOutput, tau, axis, quantity):
#    # c : column of quantity to plot in array
#    # commented version plots particles directly
#    c = cols[quantity]
#    print('quantity =',quantity)
#    print('c =',c)
#    cf   = [None, None, None, shifted_eGubser, shifted_urGubser, shifted_uetaGubser][c]
#    
#    sliceData = hydroOutput[np.where( (np.isclose(hydroOutput[:,1], 0.0, atol=1e-2)) \
#                                      & (np.isclose(hydroOutput[:,2], eta0, atol=1e-2)) )] # y == 0 ===>>> r == x
#    if quantity == 'e':
#        sliceData[:,c] *= 1000. # GeV --> MeV
#    ax.plot( sliceData[:,0], sliceData[:,c], 'r-' )
#    xpts = np.linspace(np.amin(sliceData[:,0]), np.amax(sliceData[:,0]), 1001)
#    ax.plot( xpts, cf(tau, xpts, eta0), 'b:' )
    
#===============================================================================
def plot_slice(ax, f, tau, axis, quantity):
    # c : column of quantity to plot in array
    # version below plots interpolated fields
    c = cols[quantity]
    print('quantity =',quantity)
    print('c =',c)
    cf = [None, None, None, shifted_eGubser, shifted_urGubser, shifted_uetaGubser][c]    
    if quantity == 'e':
        f[:,c] *= 1000. # GeV --> MeV
    ax.plot( f[:,0], f[:,c], 'r-' )
    ax.plot( xGrid, cf(tau, xGrid, eta0), 'b:' )

    


#===============================================================================
if __name__ == "__main__":

    # set up figure
    toPlot = ['e', 'ur', 'ueta']
    
    ncols = 3
    nrows = 1
    fig, axs = plt.subplots( ncols=ncols, nrows=nrows, figsize=(5*ncols, 5*nrows) )
    
    tau = 0.0
    
    # plot hydro output files
    for infilename in infilenames:
        # load Gubser check output files produced by hydro code
        # (eventually) use format: x [fm], y [fm], e [1/fm^4], u_x, u_y, ...
        tau = get_time_step(infilename)
        print('tau =', tau)
        
        print('Loading', infilename)
        hydroOutput = np.loadtxt(infilename, skiprows=1, usecols=(2, 3, 4, 10, 33, 35))
        print('\t - finished.')
        
        print('Eliminating irrelevant particles to save time')
        hydroOutput = hydroOutput[np.where( (np.isclose(hydroOutput[:,1], 0.0, atol=3.0*h)) \
                                          & (np.isclose(hydroOutput[:,2], eta0, atol=3.0*h)) )] # y == 0 ===>>> r == x
        print('\t - finished.')

        
        # interpolate on regular grid
        print('Smoothing fields')
        f = np.array([ evaluate_field(point) for point in grid3D ])
        print('\t - finished.')
        print('\t - f.shape =', f.shape)

        # plot comparison along y==0 slice
        for i, ax in enumerate(axs.ravel()):
            if use_log_scale and ['T','e'].count(toPlot[i]) > 0:
                ax.set_yscale('log')
            plot_slice( ax, f, tau, axisMode, toPlot[i] )
            ax.set_xlim([-4.5, 4.5])
            ax.set_xlabel(r'$x$ (fm)')
            if toPlot[i] == 'ur':
                ax.set_ylim([-3.0, 3.0])
            #ax.ylabel(r'$e$ (fm$^{-4}$)')
            
    #plt.show()
    print("tau = ", tau)
    print("float(int(tau)) = ", np.round(tau))
    print("np.isclose( tau, float(int(tau)), atol=1e-04 ) = ", np.isclose( tau, float(int(tau)), atol=1e-04 ))
    #if np.isclose( tau, np.round(tau, decimals=1), atol=1e-04 ):
    plt.savefig(outdirectory + '/yeq' + axisMode \
                                 + '_slice_tau=' + f"{tau:.2f}" \
                                 + '_eta=' + f"{eta0:.2f}" + '.pdf')
    

