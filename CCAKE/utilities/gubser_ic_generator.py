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


def main():
    check_input()
    stepx = .02
    stepy = .02
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
            TLocal = float(data[2])# GeV - Multiply by 0.465259/0.159168 to match the old normalization
            ux = data[3]
            uy = data[4]
            pixx = data[5] #Gev/fm^3
            piyy = data[6] #Gev/fm^3
            pixy = data[7] #Gev/fm^3
            pizz = data[8] #Gev/fm^5

            eps = hbarc*3.0*cpLoc*TLocal*TLocal*TLocal*TLocal/T0**4

            f.write(f"{x} {y} 0 {eps} 0 0 0 {ux} {uy} 0 0 {pixx} {pixy} 0 {piyy} 0 {pizz}\n")
    f.close()

if __name__ == "__main__":
    main()