#!/usr/bin/env python

import sys
import numpy as np

def check_input():
    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} <output>")
        sys.exit(1)

def write_header(stepx, stepy, stepEta, xmin, ymin, etamin):
    #TODO: parse from a trento input file
    f = open(sys.argv[1], "w")
    f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
    return f

t0 = 0.5 # amount of temporal shift
tau0 = 1.0
eps0 = 80
n0 = 1.0
q = 1.0

#==============================================================================
# Define Gubser model functions
def taup(tau, eta):
    return np.sqrt(tau**2 + 2.*t0*tau*np.cosh(eta)+t0**2)
#==============================================================================
def etap(tau, eta):
    return np.arctanh( tau * np.sinh(eta) / (tau * np.cosh(eta) + t0) )
#==============================================================================
def eps_a(tau, r):
    return (eps0 / tau**(4./3.)) * (2 * q)**(8./3.) / \
           (1 + 2 * q**2 * (tau**2 + r**2) + (q**2 * (tau**2 - r**2))**2)**(4./3.)
#==============================================================================
def kappa(tau, r):
    return np.arctanh(2 * q**2 * tau * r / (1 + (q * tau)**2 + (q * r)**2))
#==============================================================================
def rho(tau,r):
    return np.arcsinh( (q**2 * (tau**2 - r**2) - 1.) / (2. * q * tau) )
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

def main():
    check_input()
    stepx = 0.1
    stepy = 0.1
    stepeta = 0.025
    xmax = 5
    ymax = 5
    etamax = 2
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
                ux   = shifted_velocity_x(tau0, float(x), r, eta)
                uy   = shifted_velocity_y(tau0, float(y), r, eta)
                ueta = shifted_velocity_eta(tau0, r, eta)
                eps  = shifted_eps(tau0, r, eta)
                pixx = 0.
                piyy = 0.
                pixy = 0.
                pizz = 0.
                #print(eps)
                
                f.write(f"{x} {y} {eta} {eps} 0 0 0 {ux} {uy} {ueta} 0 {pixx} {pixy} 0 {piyy} 0 {pizz}\n")

    f.close()

if __name__ == "__main__":
    main()
