#!/usr/bin/env python

import sys
import numpy as np
from scipy.special import hyp2f1

def check_input():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output>")
        sys.exit(1)

def write_header(stepx, stepy, stepEta, xmin, ymin, etamin):
    #TODO: parse from a trento input file
    f = open(sys.argv[1], "w")
    f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
    return f

hbarc = 0.19733
t0 = 0.5 # amount of temporal shift
tau0 = 1.0
eps0 = 80 # GeV/fm^3
shearOVERs = 0.134  # i.e., specific shear viscosity (eta/s)
q = 1.0
fs = 11.0

# derived parameters
H0 = 4.0*fs**0.25*shearOVERs/3.0
T0 = eps0**0.25/hbarc

H0_by_9T0 = H0 / (9.0*T0)

#==============================================================================
# Define Gubser model functions
#==============================================================================
def kappa(tau, r):
    return np.arctanh(2 * q**2 * tau * r / (1 + (q * tau)**2 + (q * r)**2))
#==============================================================================
def rho(tau,r):
    return np.arcsinh( (q**2 * (tau**2 - r**2) - 1.) / (2. * q * tau) )
#==============================================================================
def taup(tau, eta):
    return np.sqrt(tau**2 + 2.*t0*tau*np.cosh(eta)+t0**2)
#==============================================================================
def etap(tau, eta):
    return np.arctanh( tau * np.sinh(eta) / (tau * np.cosh(eta) + t0) )
#==============================================================================
def T_a(tau, r):
    s = np.sinh(rho(tau,r))
    return (hbarc/(tau * fs**0.25)) * ( T0 / np.cosh(rho(tau, r))**(2./3.) ) \
            * ( 1.0 + H0_by_9T0 * s**3 * hyp2f1(3./2., 7./6., 5./2., -s**2) )
#==============================================================================
def eps_a(tau, r):
    return fs*T_a(tau, r)**4
#==============================================================================
def velocity_tau(tau, r):
    return np.cosh(kappa(tau, r))
#==============================================================================
def velocity_x(tau, x, r):
    return np.sinh(kappa(tau, r)) * x / r if r != 0 else 0
#==============================================================================
def velocity_y(tau, y, r):
    return np.sinh(kappa(tau, r)) * y / r if r != 0 else 0
#==============================================================================
def shifted_eps(tau, r, eta):
    return eps_a(taup(tau, eta), r)
#==============================================================================
def shifted_velocity_x(tau, x, r, eta):
    return velocity_x(taup(tau, eta), x, r)
#==============================================================================
def shifted_velocity_y(tau, y, r, eta):
    return velocity_y(taup(tau, eta), y, r)
#==============================================================================
def shifted_velocity_eta(tau, r, eta):
    return -velocity_tau(taup(tau, eta), r) * ( t0 * np.sinh(eta) / (tau * taup(tau, eta)) )
#==============================================================================
#def jacobian(tau, r, eta):
#    return np.array([[(tau+t0*np.cosh(eta))/taup(tau,eta), t0*tau*np.sinh(eta)/taup(tau,eta)],
#                     [t0*tau*np.sinh(eta)/taup(tau,eta)**2, tau*(tau+t0*np.cosh(eta))/taup(tau,eta)**2]])
#==============================================================================
def jacobian(tau, eta):
    c, s, tp = np.cosh(eta), np.sinh(eta), taup(tau,eta)
    return np.array([[(tau+t0*c)/tp,  0, 0,      -t0*s],
                     [0,              1, 0,          0],
                     [0,              0, 1,          0],
                     [-t0*s/(tau*tp), 0, 0, 1+t0*c/tau]]).T
#==============================================================================
def pimunu(tau, x, y, r):
    shear = H0*eps_a(tau, r)**0.75
    prefactor = 2.*shear*np.tanh(rho(tau, r))/(3.*tau**4) # N.B. - missing minus sign relative to 2503.XXXXX
    ux = velocity_x(tau, x, r)
    uy = velocity_y(tau, y, r)
    utau = velocity_tau(tau, r)
    return prefactor * np.array([[ux**2+uy**2, ux*utau, uy*utau,           0],
                                 [ux*utau,     1+ux**2,   ux*uy,           0],
                                 [uy*utau,       ux*uy, 1+uy**2,           0],
                                 [0,                  0,       0, -2./tau**2]])
#==============================================================================
def shifted_pimunu(tau, x, y, r, eta):
    j = jacobian(tau, eta)
    return j.T @ pimunu(taup(tau,eta), x, y, r) @ j
#==============================================================================

def main():
    #print('eps_a(1,0) =', eps_a(1.0,0.0))
    #if True:
    #    exit(1)
       
    
    check_input()
    stepx = 0.05
    stepy = 0.05
    stepeta = 0.025
    xmax = 5.0
    ymax = 5.0
    etamax = 2.0
    xmin = -xmax
    ymin = -ymax
    etamin = -etamax
    hbarc = 0.1973269804 
    # scale for conformal EoS
    Nc    = 3.0 # three colors
    Nf    = 2.5 # u+d massless, s 'half massless'
    cpLoc = np.pi*np.pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0 #In units of fm^-4
    T0   = hbarc # T0 = 1 fm^-1 = .197 GeV

    #f = open(sys.argv[1], "r")
    #lines = f.readlines()
    #f.close()
    
    f = write_header(stepx, stepy, stepeta, xmin, ymin, etamin)

    for x in np.arange(xmin, xmax+stepx, stepx):
        for y in np.arange(ymin, ymax+stepy, stepy):
            r = np.sqrt(float(x)**2 + float(y)**2)
            for eta in np.arange(etamin, etamax+stepeta, stepeta):
                ux       = shifted_velocity_x(tau0, float(x), r, eta)
                uy       = shifted_velocity_y(tau0, float(y), r, eta)
                ueta     = shifted_velocity_eta(tau0, r, eta)
                eps      = shifted_eps(tau0, r, eta)
                piM      = shifted_pimunu(tau0, x, y, r, eta)
                pixx     = piM[1,1]
                piyy     = piM[2,2]
                pixy     = piM[1,2]
                pixeta   = piM[1,3]
                piyeta   = piM[2,3]
                pietaeta = piM[3,3]
                #utau = np.sqrt(1.0+ux**2+uy**2+(tau0*ueta)**2)
                #print('----------------------------------------')
                #np.set_printoptions(precision=8)
                #print('u =',utau,ux,uy,tau0**2*ueta)
                #print('pi =',piM)
                #if True:
                #    exit(1)
                #print('tr(pi) =',piM[0,0]-piM[1,1]-piM[2,2]-tau0**2*piM[3,3])
                #print('u*pi(0) =',utau*piM[0,0]-ux*piM[0,1]-uy*piM[0,2]-tau0**2*ueta*piM[0,3])
                #print('u*pi(1) =',utau*piM[1,0]-ux*piM[1,1]-uy*piM[1,2]-tau0**2*ueta*piM[1,3])
                #print('u*pi(2) =',utau*piM[2,0]-ux*piM[2,1]-uy*piM[2,2]-tau0**2*ueta*piM[2,3])
                #print('u*pi(3) =',utau*piM[3,0]-ux*piM[3,1]-uy*piM[3,2]-tau0**2*ueta*piM[3,3])
                #print(eps)
                
                f.write(f"{x} {y} {eta} {eps} 0 0 0 {ux} {uy} {ueta} 0 {pixx} {pixy} {pixeta} {piyy} {piyeta} {pietaeta}\n")

    f.close()

if __name__ == "__main__":
    main()
