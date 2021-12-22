import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys

use_semi_analytic = True
use_log_scale = False

axisMode = sys.argv[1]

q  = 1.
e0 = 1.0
if use_semi_analytic:
    e0 = 9126.*np.pi**2/3125. # normalization needed to get initial T0 = 1.2 1/fm
rhoB0, rhoS0, rhoQ0 = 0.5, 0.5, 0.5

Nc = 3.
Nf = 2.5
cp = (2.*(Nc**2-1.) + 3.5*Nc*Nf)*np.pi**2/90.

quantities = ['T','e','ux','uy','pixx','piyy','pixy','pizz','rhoB','rhoS','rhoQ']
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
                            / ( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )
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
    if axisMode == '0':
        return np.loadtxt('./ac/Initial_Profile_tau=1fm.dat'),  \
               np.loadtxt('./ac/y=0_tau=1.2_SemiAnalytic.dat'), \
               np.loadtxt('./ac/y=0_tau=1.5_SemiAnalytic.dat'), \
               np.loadtxt('./ac/y=0_tau=2.0_SemiAnalytic.dat')
    else:
        return np.loadtxt('./ac/Initial_Profile_tau=1fm.dat'),  \
               np.loadtxt('./ac/y=x_tau=1.2_SemiAnalytic.dat'), \
               np.loadtxt('./ac/y=x_tau=1.5_SemiAnalytic.dat'), \
               np.loadtxt('./ac/y=x_tau=2.0_SemiAnalytic.dat')

#===============================================================================
def plot_slice(ax, hydroOutput, tau, axis, quantity):
    # c : column of quantity to plot in array
    c = cols[quantity]
    if axis == '0':
        yEqAxisData = hydroOutput[np.where( np.abs(hydroOutput[:,1]) < 1e-6 )]
        ax.plot( yEqAxisData[:,0], yEqAxisData[:,c], 'r-' )
        if not use_semi_analytic:
            cf   = [None, None, TGubser, eGubser, urGubser, urGubser, \
                    None, None, None, None, rhoBGubser, rhoSGubser, rhoQGubser][c]
            xpts = np.linspace(np.amin(yEqAxisData[:,0]), np.amax(yEqAxisData[:,0]), 1001)
            ax.plot( xpts, cf(tau, xpts), 'b--' )
    elif axis == 'x':
        yeqxData = hydroOutput[np.where( np.isclose( hydroOutput[:,0], hydroOutput[:,1] ) )]
        rpts = np.sqrt(yeqxData[:,0]**2 + yeqxData[:,1]**2)
        ax.plot( rpts, yeqxData[:,c], 'r-' )
        if not use_semi_analytic:
            cf   = [None, None, TGubser, eGubser, urGubser, urGubser, \
                    None, None, None, None, rhoBGubser, rhoSGubser, rhoQGubser][c]
            rpts = np.linspace(0.0, np.amax(rpts), 1001)
            ax.plot( rpts, cf(tau, rpts), 'b--' )
    

#===============================================================================
def get_time_step(filename):
    return float((open(filename)).readline()) # first line is header containing just timestep

#===============================================================================
if __name__ == "__main__":
    # load files where semi-analytic calculations are stored that we can compare against (not used yet)
    ic, yEqAxis_tau1_2, yEqAxis_tau1_5, yEqAxis_tau2_0 = \
        load_semi_analytic_files()

    # set up figure
    toPlot = ['T', 'e', 'ux', 'rhoB', 'rhoS', 'rhoQ']
    if use_semi_analytic:
        toPlot = ['e', 'ux', 'pixx', 'piyy', 'pixy', 'pizz']
    
    ncols = len(toPlot)//2
    nrows = 2
    fig, axs = plt.subplots( ncols=ncols, nrows=nrows, figsize=(5*ncols, 5*nrows) )

    # plot hydro output files
    for checkfile in sys.argv[2:]:
        # load Gubser check output files produced by hydro code
        # (eventually) use format: x [fm], y [fm], e [1/fm^4], u_x, u_y, ...
        tau = get_time_step(checkfile)
        print('tau=', tau)
        hydroOutput = np.loadtxt( checkfile, skiprows=1 )
        
        # plot comparison along y==0 slice
        for i, ax in enumerate(axs.ravel()):
            if use_log_scale and ['T','e','rhoB','rhoS','rhoQ'].count(toPlot[i]) > 0:
                ax.set_yscale('log')
            plot_slice( ax, hydroOutput, tau, axisMode, toPlot[i] )
            ax.set_xlim([-4.75, 4.75])
            ax.set_xlabel(r'$r$ (fm)')
            #ax.ylabel(r'$e$ (fm$^{-4}$)')
            
    # plot results of semi-analytic calculation if desired
    if use_semi_analytic:
        quantities = ['e','ux','uy','pixx','piyy','pixy','pizz']
        cols = dict(zip(quantities,range(2,len(quantities)+2)))
        for i, ax in enumerate(axs.ravel()):
            if use_log_scale and ['T','e','rhoB','rhoS','rhoQ'].count(toPlot[i]) > 0:
                ax.set_yscale('log')
            c = cols[toPlot[i]]
            if toPlot[i] == 'e':
                for data in [ic[np.where(np.abs(ic[:,1])<1e-10)], \
                             yEqAxis_tau1_2, yEqAxis_tau1_5, yEqAxis_tau2_0]:
                    data[:,c] /= 0.1973
                    ax.plot( data[:,0], eFromT(data[:,c]), 'b--' )
            else:
                for data in [ic[np.where(np.abs(ic[:,1])<1e-10)], \
                             yEqAxis_tau1_2, yEqAxis_tau1_5, yEqAxis_tau2_0]:
                    if ['pixx','piyy','pixy','pizz'].count(toPlot[i]) > 0:
                        data[:,c] /= 0.1973
                    ax.plot( data[:,0], data[:,c], 'b--' )
    
    #plt.show()
    plt.savefig('./yeq' + axisMode + '_slice.pdf')
    
