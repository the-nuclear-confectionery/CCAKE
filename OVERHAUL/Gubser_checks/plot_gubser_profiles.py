#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

q  = 1.0
#e0 = 1.0
e0 = 9126.0*np.pi**2/3125.0
rhoB0, rhoS0, rhoQ0 = 0.5, 0.5, 0.5

#===============================================================================
def eGubser(tau, r):
    return (e0/tau**4)*( (2.*q*tau)**(8./3.)
                        / ( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )**(4./3.)
                         )

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
def plot_yeq0_slice(ax, hydroOutput, tau):
    yeq0Data = hydroOutput[:,np.where( np.abs(hydroOutput[:,1]) < 1e-6 )]
    ax.plot( yeq0Data[:,0], yeq0Data[:,2], 'ro' )
    xpts = np.linspace(np.amin(yeq0Data[:,0]),np.amax(yeq0Data[:,0]), 1001)
    ax.plot( xpts, eGubser(tau, xpts), 'b-' )
    

#===============================================================================
def get_time_step(filename):
    return float((open(filename)).readline()) # first line is header containing just timestep

#===============================================================================
if __name__ == "__main__":
    # load files where semi-analytic calculations are stored that we can compare against (not used yet)
    #ic, yEq0_tau1_2, yEq0_tau1_5, yEq0_tau2_0,  \
    #    yEqx_tau1_2, yEqx_tau1_5, yEqx_tau2_0 = \
    #    load_semi_analytic_files()
        
    # set up figure
    fig, ax = plt.subplots( nrows=1, ncols=1 )

    for checkfile in sys.argv[1:]:
        # load Gubser check output files produced by hydro code
        # (eventually) use format: x [fm], y [fm], e [1/fm^4], u_x, u_y, ...
        tau = get_time_step(checkfile)
        hydroOutput = np.loadtxt( checkfile, skiprows=1 )
        
        # plot comparison along y==0 slice
        plot_yeq0_slice( ax, hydroOutput, tau )
    
    #plt.show()
    plt.savefig('./yeq0_slice.pdf')
    
