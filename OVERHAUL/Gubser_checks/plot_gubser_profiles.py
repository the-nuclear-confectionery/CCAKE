#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

use_semi_analytic = True

q  = 1.
#e0 = 1.0
e0 = 9126.*np.pi**2/3125. # normalization needed to get initial T0 = 1.2 1/fm
rhoB0, rhoS0, rhoQ0 = 0.5, 0.5, 0.5

Nc = 3.
Nf = 2.5
cp = (2.*(Nc**2-1.) + 3.5*Nc*Nf)*np.pi**2/90.

quantities = ['T','e','ux','uy','pixx','piyy','pixy','pizz']
cols = dict(zip(quantities,range(2,len(quantities)+2)))


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
def chargeGubser(tau, r):
    return (1./tau**3)*( (4.*q**2*tau**2)
                            / ( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )**2
                            )

#===============================================================================
def rhoBGubser(tau, r):
    return rhoB0 * chargeGubser(tau, r)

#===============================================================================
def rhoSGubser(tau, r):
    return rhoS0 * chargeGubser(tau, r)

#===============================================================================
def rhoQGubser(tau, r):
    return rhoQ0 * chargeGubser(tau, r)

#===============================================================================
def load_semi_analytic_files():
    return np.loadtxt('./ac/Initial_Profile_tau=1fm.dat'),  \
           np.loadtxt('./ac/y=0_tau=1.2_SemiAnalytic.dat'), \
           np.loadtxt('./ac/y=0_tau=1.5_SemiAnalytic.dat'), \
           np.loadtxt('./ac/y=0_tau=2.0_SemiAnalytic.dat'), \
           np.loadtxt('./ac/y=x_tau=1.2_SemiAnalytic.dat'), \
           np.loadtxt('./ac/y=x_tau=1.5_SemiAnalytic.dat'), \
           np.loadtxt('./ac/y=x_tau=2.0_SemiAnalytic.dat')

#===============================================================================
def plot_slice(ax, hydroOutput, tau, axis, quantity):
    # c : column of quantity to plot in array
    c = cols[quantity]
    if axis == '0':
        yeq0Data = hydroOutput[np.where( np.abs(hydroOutput[:,1]) < 1e-6 )]
        ax.plot( yeq0Data[:,0], yeq0Data[:,c], 'r-' )
        if not use_semi_analytic:
            cf   = [None, None, TGubser, eGubser, urGubser, urGubser, None, None, None, None][c]
            xpts = np.linspace(np.amin(yeq0Data[:,0]), np.amax(yeq0Data[:,0]), 1001)
            ax.plot( xpts, cf(tau, xpts), 'b--' )
    elif axis == 'x':
        yeqxData = hydroOutput[np.where( np.isclose( hydroOutput[:,0], hydroOutput[:,1] ) )]
        rpts = np.sqrt(hydroOutput[:,0]**2 + hydroOutput[:,1])
        ax.plot( rpts, yeqxData[:,c], 'r-' )
        if not use_semi_analytic:
            cf   = [None, None, TGubser, eGubser, urGubser, urGubser, None, None, None, None][c]
            rpts = np.linspace(0.0, np.amax(rpts), 1001)
            ax.plot( rpts, cf(tau, rpts), 'b--' )
    

#===============================================================================
def get_time_step(filename):
    return float((open(filename)).readline()) # first line is header containing just timestep

#===============================================================================
if __name__ == "__main__":
    # load files where semi-analytic calculations are stored that we can compare against (not used yet)
    ic, yEq0_tau1_2, yEq0_tau1_5, yEq0_tau2_0,  \
        yEqx_tau1_2, yEqx_tau1_5, yEqx_tau2_0 = \
        load_semi_analytic_files()

    # set up figure
    toPlot = ['e', 'ux', 'pixx']
    fig, axs = plt.subplots( ncols=len(toPlot), nrows=1, figsize=(5*len(toPlot), 5) )

    # plot hydro output files
    for checkfile in sys.argv[1:]:
        # load Gubser check output files produced by hydro code
        # (eventually) use format: x [fm], y [fm], e [1/fm^4], u_x, u_y, ...
        tau = get_time_step(checkfile)
        hydroOutput = np.loadtxt( checkfile, skiprows=1 )
        
        # plot comparison along y==0 slice
        for i, ax in enumerate(axs.ravel()):
            #ax.set_yscale('log')
            plot_slice( ax, hydroOutput, tau, '0', toPlot[i] )
            ax.set_xlim([-4.75, 4.75])
            ax.set_xlabel(r'$x$ (fm)')
            #ax.ylabel(r'$e$ (fm$^{-4}$)')
            
    # plot results of semi-analytic calculation if desired
    if use_semi_analytic:
        quantities = ['e','ux','uy','pixx','piyy','pixy','pizz']
        cols = dict(zip(quantities,range(2,len(quantities)+2)))
        for i, ax in enumerate(axs.ravel()):
            c = cols[toPlot[i]]
            if toPlot[i] == 'e':
                for data in [ic, yEq0_tau1_2, yEq0_tau1_5, yEq0_tau2_0]:
                    data[:,c] /= 0.1973
                    ax.plot( data[:,0], eFromT(data[:,c]), 'b--' )
            else:
                for data in [ic, yEq0_tau1_2, yEq0_tau1_5, yEq0_tau2_0]:
                    if ['pixx','piyy','pixy','pizz'].count(toPlot[i]) > 0:
                        data[:,c] /= 0.1973
                    ax.plot( data[:,0], data[:,c], 'b--' )
    
    #plt.show()
    plt.savefig('./yeq0_slice.pdf')
    
