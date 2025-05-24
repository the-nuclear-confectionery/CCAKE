import sys
import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import logging
from scipy.special import hyp2f1

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
#from variable_conversions import HBARC  # noqa
import my_plotting as myplt  # noqa

from matplotlib.colors import ListedColormap

custom_colors = ['#FF5F05', '#FF8C42', '#13294B', '#009FD4', '#8FC1DE', '#707372']
#['#FF5F05', '#13294B', '#009FD4', '#8FC1DE', '#707372']
cmap = ListedColormap(custom_colors)

logging.basicConfig(
    level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s"
)
logger = logging.getLogger(__name__)

plot_index = 1
# analytic_style = {'ls': '-', 'lw': 2}
# sim_style = {'ls':'-.','lw':3.,'alpha':.5}
# sim_style = {'facecolors': 'none'}

#time_list = np.arange(1.00, 1.50, 0.1)  # Use this to focus on before FO
#time_list=['1.00', '1.10', '1.20','1.30', '1.40', '1.50']# Use this to focus on before FO
time_list=['1.00', '1.10', '1.20']# Use this to focus on before FO
filter_criteria = 'abs(phi - 3.141592653589793/4.) < 1.e-2'

# cmap = myplt.get_cmap(len(time_list), 'cividis')

#mpl.rcParams['text.usetex'] = True
plt.rcParams.update({
    "text.usetex": True,
    "text.latex.preamble": r'\usepackage{amsmath}'
    })

dpi = 150
fig, ax = plt.subplot_mosaic([['e_at_eta0', 'ur_at_eta0', 'ueta_at_eta0', 'cbar'],
                              ['pixx_at_eta0', 'piyy_at_eta0', 'pietaeta_at_eta0', 'cbar']],
                             width_ratios=[1, 1, 1, 0.1],
                             figsize=np.array([7 * 3, 7 * 2]),
                             constrained_layout=True)




#============================================================
#========================
# my stuff below this line
use_log_scale = False


#axisMode = sys.argv[1]
#outdirectory = sys.argv[2]
#infilenames = sys.argv[3:]

#n_timesteps = min([len(infilenames),1000])
n_timesteps = len(time_list)


hbarc = 0.1973269804              # GeV*fm
t0 = 0.0 # amount of temporal shift
q  = 1.
e0 = 80000.0 # units of MeV*fm (cf. ideal_2+1d_gubser_ic_generator.py script)
e0 *= 0.001  # convert to units of GeV*fm

shearOVERs = 0.134  # i.e., specific shear viscosity (eta/s)
Nc = 3.
Nf = 2.5
cp = (2.*(Nc**2-1.) + 3.5*Nc*Nf)*np.pi**2/90.
#fs = 11.0                         # dimensionless
fs = cp                         # dimensionless

# derived parameters
H0 = 4.0*fs**0.25*shearOVERs/3.0
T0 = 0.25*e0**0.25*fs**0.75     # dimensionless

H0_by_9T0 = H0 / (9.0*T0)

h = 0.1  # 3D cubic spline smoothing scale
knorm  = 1./(np.pi*h**3)

# points at which to compare exact and numerical solutions
eta0, eta1 = 0.0, 1.0
xGrid = np.linspace(-5.0, 5.0, 1001)
yGrid = 0.0*xGrid
etaGrid = 0.0*xGrid
#grid3D = np.c_[ xGrid, yGrid, etaGrid ]

quantities = ['e','ur','ueta']
cols = dict(zip(quantities,range(3,len(quantities)+3)))
print('cols=',cols)


#==============================================================================
#==============================================================================
#==============================================================================
#==========================   MY STUFF   ======================================
#==============================================================================
#==============================================================================
#==============================================================================
def rho(tau,r):
    return np.arcsinh( (q**2 * (tau**2 - r**2) - 1.) / (2. * q * tau) )
#==============================================================================
#def eGubser(tau, r):
#    return (e0/tau**4)*( (2.*q*tau)**(8./3.)
#                        / ( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )**(4./3.)
#                         )
#===============================================================================
#def TGubser(tau, r):
#    return hbarc*( eGubser(tau, r) / (3.*cp*hbarc) )**0.25
#==============================================================================
def T_a(tau, r):
    s = np.sinh(rho(tau,r))
    return (hbarc/(tau * fs**0.25)) * ( T0 / np.cosh(rho(tau, r))**(2./3.) ) \
            * ( 1.0 + H0_by_9T0 * s**3 * hyp2f1(3./2., 7./6., 5./2., -s**2) )
#==============================================================================
def eps_a(tau, r):
    return fs*T_a(tau, r)**4/hbarc**3
#===============================================================================
def urGubser(tau, r):
    return 2.0*q**2*r*tau / np.sqrt( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )
#==============================================================================
def velocity_x(tau, x, r):
    return 2.0*q**2*x*tau / np.sqrt( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )
#==============================================================================
def velocity_y(tau, y, r):
    return 2.0*q**2*y*tau / np.sqrt( 1. + 2.*q**2*(tau**2 + r**2) + q**4*(tau**2 - r**2)**2 )
#===============================================================================
def utauGubser(tau, r):
    return np.sqrt(1.0 + urGubser(tau, r)**2)
#===============================================================================
def uetaGubser(tau, r):
    return 0.0
#==============================================================================
def pimunuGubser(tau, x, y, r):
    #shear = H0*eGubser(tau, r)**0.75
    shear = H0*(hbarc**3*eps_a(tau, r))**0.75
    prefactor = 4.*shear*np.tanh(rho(tau, r))/(3.*tau*hbarc**2) # N.B. - missing minus sign and fixed power of tau relative to 2503.XXXXX
    ux = velocity_x(tau, x, r)
    uy = velocity_y(tau, y, r)
    utau = utauGubser(tau, r)
    return prefactor * np.array([[ux**2+uy**2, ux*utau, uy*utau,           0],
                                 [ux*utau,     1+ux**2,   ux*uy,           0],
                                 [uy*utau,       ux*uy, 1+uy**2,           0],
                                 [0,                  0,       0, -2./tau**2]])
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
def evaluate_field(data, r):
    global h
    neighbors = data[ ( r[0]-data[:,0])**2 \
                      +(r[1]-data[:,1])**2 \
                      +(r[2]-data[:,2])**2 <= 4.0*h**2 ]
    weights = kernel( np.sqrt( ( r[0]-neighbors[:,0])**2 \
                               +(r[1]-neighbors[:,1])**2 \
                               +(r[2]-neighbors[:,2])**2 )/h )
    if np.sum(weights) < 1e-10:
        return np.pad(r, (0, neighbors.shape[1] - 3), 'constant')
    else:
        return np.sum( neighbors*weights[:, np.newaxis], axis=0 ) / (np.sum(weights)+1e-10)

#==============================================================================
#==============================================================================
#==============================================================================
#==========================   MY STUFF   ======================================
#==============================================================================
#==============================================================================
#==============================================================================



def plot_analytic_sol():
    for ii, t in enumerate(time_list):
        if ii % plot_index != 0:
            continue

        tau = float(t)
        cf = {'e_at_eta0': eps_a,
               'ur_at_eta0': urGubser,
               'ueta_at_eta0': uetaGubser,
               'pixx_at_eta0': eps_a,
               'piyy_at_eta0': eps_a,
               'pietaeta_at_eta0': eps_a}

        analytic_style = {'ls': '-', 'lw': 7, 'alpha': 0.4}

        for key in ['e_at_eta0', 'ur_at_eta0', 'ueta_at_eta0']:
            ax[key].plot( xGrid, cf[key](tau, xGrid), color=cmap(ii), **analytic_style)
            print(cf[key](tau, xGrid).shape)

        piM = np.array([pimunuGubser(tau, x, 0.0, x) for x in xGrid])
        print(piM.shape)
        #if True:
        #    exit(1)
        #for key in ['pixx_at_eta0', 'piyy_at_eta0', 'pietaeta_at_eta0']:
        #    ax[key].plot( xGrid, cf[key](tau, xGrid, eta0), color=cmap(ii), **analytic_style)
        ax['pixx_at_eta0'].plot( xGrid, piM[:,1,1], color=cmap(ii), **analytic_style)
        ax['piyy_at_eta0'].plot( xGrid, piM[:,2,2], color=cmap(ii), **analytic_style)
        ax['pietaeta_at_eta0'].plot( xGrid, piM[:,3,3], color=cmap(ii), **analytic_style)



def read_sim(sim_result_folder):
    dt = .001
    for ii, t in enumerate(time_list):
        if ii % plot_index != 0:
            continue
        col_names = [
            'id',
            't',
            'x',
            'y',
            'p',
            'T',
            'muB',
            'muS',
            'muQ',
            'e',
            'rhoB',
            'rhoS',
            'rhoQ',
            's',
            's_smoothed',
            's_specific',
            'sigma',
            'spec_s',
            'stauRelax',
            'bigTheta',
            '??',
            '??2',
            'pi00',
            'pixx',
            'piyy',
            'pixy',
            't2pi33',
            'v1',
            'v2',
            'gamma',
            'frz',
            'eos']
        idx = int(np.round((float(t) - 1) / dt) / 100)
        inp_path = os.path.join(sim_result_folder, f'system_state_{idx}.dat')
        print(inp_path)

        print('Loading', inp_path)
        hydroOutput = np.loadtxt(inp_path, skiprows=1, usecols=(2, 3, 4, 10, 33, 35, 25, 26, 32))
        print('\t - finished.')

        #print('\t - hydroOutput.shape =', hydroOutput.shape)
        print(hydroOutput[np.where(    (np.isclose(hydroOutput[:,0], 0.0, atol=1e-3))  \
                                     & (np.isclose(hydroOutput[:,1], 0.0, atol=1e-3)) \
                                     & (np.isclose(hydroOutput[:,2], 0.0, atol=1e-3)))])

        print('Eliminating irrelevant particles to save time')
        hydroOutput_at_eta0 = hydroOutput[np.where( (np.isclose(hydroOutput[:,1], 0.0, atol=3.0*h)) \
                                                    & (np.isclose(hydroOutput[:,2], eta0, atol=3.0*h)) )] # y == 0 ===>>> r == x
        #hydroOutput_at_eta1 = hydroOutput[np.where( (np.isclose(hydroOutput[:,1], 0.0, atol=3.0*h)) \
        #                                            & (np.isclose(hydroOutput[:,2], eta1, atol=3.0*h)) )] # y == 0 ===>>> r == x
        print('\t - finished.')
        #print('\t - hydroOutput.shape =', hydroOutput.shape)
        #print('\t - hydroOutput_at_eta0.shape =', hydroOutput_at_eta0.shape)
        #print('\t - hydroOutput_at_eta1.shape =', hydroOutput_at_eta1.shape)
        #print(hydroOutput_at_eta0)
        #print(hydroOutput_at_eta1)

        # interpolate on regular grid
        print('Smoothing fields')
        f0 = np.array([ evaluate_field(hydroOutput_at_eta0, point) for point in np.c_[ xGrid, yGrid, eta0 + etaGrid ] ])
        #f1 = np.array([ evaluate_field(hydroOutput_at_eta1, point) for point in np.c_[ xGrid, yGrid, eta1 + etaGrid ] ])
        print('\t - finished.')
        #print('\t - f0.shape =', f0.shape)
        #print('\t - f1.shape =', f1.shape)
        #print(f0)
        #print(f1)

        #if True:
        #    exit(1)

        #print(f0[f0[:,0] > 0, 3])
        #print(f1[f1[:,0] > 0, 3])

        sim_style = {'color': cmap(ii), 'ls': ':', 'lw': 2.5}

        ax['e_at_eta0'].plot(   f0[:,0], f0[:,3], **sim_style)
        ax['ur_at_eta0'].plot(  f0[:,0], f0[:,4], **sim_style)
        ax['ueta_at_eta0'].plot(f0[:,0], f0[:,5], **sim_style)
        ax['pixx_at_eta0'].plot(   f0[:,0], hbarc*f0[:,6], **sim_style)
        ax['piyy_at_eta0'].plot(  f0[:,0], hbarc*f0[:,7], **sim_style)
        ax['pietaeta_at_eta0'].plot(f0[:,0], hbarc*f0[:,8], **sim_style)
        print('Center:',f0[np.abs(f0[:,0])<1e-4])
        #ax['e_at_eta1'].plot(   f1[:,0], f1[:,3], **sim_style)
        #ax['ur_at_eta1'].plot(  f1[:,0], f1[:,4], **sim_style)
        #ax['ueta_at_eta1'].plot(f1[:,0], f1[:,5], **sim_style)



def beautify():
    # fig.set_tight_layout(True)
    #ylabels = {'e_at_eta0': r'$\mathcal E$ [GeV/fm$^3$]',
    #           'ur_at_eta0': r'$u^r$',
    #           'ueta_at_eta0': r'$u^\eta$ [1/fm]',
    #           'e_at_eta1': r'$\mathcal E$ [GeV/fm$^3$]',
    #           'ur_at_eta1': r'$u^r$',
    #           'ueta_at_eta1': r'$u^\eta$ [1/fm]'}
    ylabels = {'e_at_eta0': r'$\varepsilon_\text{sG}$ [GeV/fm$^3$]',
               'ur_at_eta0': r'$u^r_\text{sG}$',
               'ueta_at_eta0': r'$u^\eta_\text{sG}$ [1/fm]',
               'pixx_at_eta0': r'$\pi_{xx,\text{sG}}$ [GeV/fm$^3$]',
               'piyy_at_eta0': r'$\pi_{yy,_\text{sG}}$ [GeV/fm$^3$]',
               'pietaeta_at_eta0': r'$\pi_{\eta\eta,_\text{sG}}$ [GeV/fm$^5$]'}
    for key in ax.keys():
        if key == 'cbar':
            continue
        myplt.costumize_axis(ax=ax[key],
                             x_title=r'$r$ [fm]',
                             y_title=ylabels[key])
        ax[key].set_xlim(0, 4.5)

    myplt.costumize_axis(ax=ax['cbar'], x_title='', y_title='')

    tau_min = float(time_list[0])
    tau_max = float(time_list[-1])
    tau_list = [float(t) for t in time_list]

    # deal with the colorbar
    delta_time = tau_list[1] - tau_list[0]
    boundaries = np.arange(tau_min - delta_time / 2,
                           tau_max + 3 * delta_time / 2, delta_time)

    norm = mpl.colors.Normalize(vmin=tau_min, vmax=tau_max)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks=tau_list, label=r'$\tau$ (fm/c)',
                        cax=ax['cbar'],
                        boundaries=boundaries)
    style1 = {'ls': ':', 'lw': 3}
    style2 = {'ls': '-', 'lw': 3}
    analytic_style = {'ls': '-', 'lw': 7, 'alpha': 0.4}
    ax['ur_at_eta0'].plot(
        [],
        [],
        **analytic_style,
        label='Analytic')
        # edgecolors=cmap(0))
    ax['ur_at_eta0'].plot([], [], **style1, label='CCAKE', color=cmap(0))

    ax['e_at_eta0'].set_ylim(-0.5, 10.0)
    ax['ur_at_eta0'].set_ylim(0,2.0)
    ax['ueta_at_eta0'].set_ylim(-0.75, 0.05)
    ax['pixx_at_eta0'].set_ylim(-1.0, 1.0)
    ax['piyy_at_eta0'].set_ylim(-1.0, 1.0)
    ax['pietaeta_at_eta0'].set_ylim(-1.0, 1.0)

    for name, label in zip(ylabels.keys(),
                           ['a', 'b', 'c', 'd', 'e', 'f']):
        ax[name].text(
            0.93,
            0.93,
            f'({label})',
            transform=ax[name].transAxes,
            fontsize=20,
            bbox={'boxstyle': 'round', 'facecolor': 'white'},
            horizontalalignment='center',
        )

    ax['e_at_eta0'].text(
        0.8,
        0.5,
        r'$\eta = 0$',
        transform=ax['e_at_eta0'].transAxes,
        fontsize=30,
        #bbox={'boxstyle': 'round', 'facecolor': 'white'},
        horizontalalignment='center',
    )

    #ax['e_at_eta1'].text(
    #    0.8,
    #    0.5,
    #    r'$\eta = 1$',
    #    transform=ax['e_at_eta1'].transAxes,
    #    fontsize=30,
    #    #bbox={'boxstyle': 'round', 'facecolor': 'white'},
    #    horizontalalignment='center',
    #)

    ax['ur_at_eta0'].legend(loc='upper center', fontsize=18)


if __name__ == '__main__':
    #analytical_folder = sys.argv[1]
    #simulation_folder1 = sys.argv[2]
    plot_analytic_sol()
    read_sim('./V1_results_2+1d_viscous_Gubser/')
    beautify()

    #if (len(sys.argv) < 5):
    #    fig.savefig('inverter_comparison.pdf')
    #else:
    #    fig.savefig(sys.argv[3])
    fig.savefig('tmp_2_1d.pdf')
