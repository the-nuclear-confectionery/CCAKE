#!/usr/bin/env python

import sys
import numpy as np

def check_input():
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <input> <output>")
        print("""The input is the semi-analytic solution at tau0 = 1 fm/c.
                 It is a tabular file with the columns
                 x  y  TLocal  ux  uy  pixx  piyy  pixy pizz;
              """)
        sys.exit(1)

def write_header(stepx, stepy, stepEta, xmin, ymin, etamin):
    #TODO: parse from a trento input file
    f = open(sys.argv[2], "w")
    f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
    return f


tau0 = 1.0
eps0 = 80
n0 = 1.0
q = 1.0

# Define Gubser model functions
def eps_a(r):
    return (eps0 / tau0**(4./3.)) * (2 * q)**(8./3.) / \
           (1 + 2 * q**2 * (tau0**2 + r**2) + (q**2 * (tau0**2 - r**2))**2)**(4./3.)
def kappa(r):
    return np.arctanh(2 * q**2 * tau0 * r / (1 + (q * tau0)**2 + (q * r)**2))
def velocity_x(x, r):
    return np.sinh(kappa(r)) * x / r if r != 0 else 0
def velocity_y(y, r):
    return np.sinh(kappa(r)) * y / r if r != 0 else 0

def main():
    check_input()
    stepx = .025
    stepy = .025
    xmax = 5
    ymax = 5
    xmin = -xmax
    ymin = -ymax
    hbarc = 0.1973269804 
    # scale for conformal EoS
    Nc    = 3.0 # three colors
    Nf    = 2.5 # u+d massless, s 'half massless'
    cpLoc = np.pi*np.pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0 #In units of fm^-4
    T0   = hbarc # T0 = 1 fm^-1 = .197 GeV

    f = open(sys.argv[1], "r")
    lines = f.readlines()
    f.close()

    f = write_header(stepx, stepy, .1, xmin, ymin, -.1)

    for ix, line in enumerate(lines):
        if line[0] == "#":
            continue
        else:
            data = line.split()
            x = data[0] #fm
            y = data[1] #fm
            r = np.sqrt(float(x)**2 + float(y)**2)
            ux = velocity_x(float(x), r)
            uy = velocity_y(float(y), r)
            eps = eps_a(r)
            pixx = 0.
            piyy = 0.
            pixy = 0.
            pizz = 0.
            #print(eps)

            f.write(f"{x} {y} 0 {eps} 0 0 0 {ux} {uy} 0 0 {pixx} {pixy} 0 {piyy} 0 {pizz}\n")
    f.close()

if __name__ == "__main__":
    main()