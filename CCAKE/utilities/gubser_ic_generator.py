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
    stepx = .05
    stepy = .05
    xmax = 5
    ymax = 5
    xmin = -xmax
    ymin = -ymax

    # Gubser parameters
    q     = 1.0 # 1/fm
    e0    = 9126.0*np.pi*np.pi/3125.0 # 1/fm^4
                                      # use this normalization to compare with semi-
                                      # analytic calculation in Phys. Rev. C 91, 014903
    rhoB0 = 0.0 # 1/fm^3
    rhoQ0 = 0.0 # 1/fm^3
    rhoS0 = 0.0 # 1/fm^3
    
    # scale for conformal EoS
    Nc    = 3.0 # three colors
    Nf    = 2.5 # u+d massless, s 'half massless'
    cpLoc = np.pi*np.pi*(2.0*(Nc*Nc-1.0)+(7.0/2.0)*Nc*Nf)/90.0

    f = open(sys.argv[1], "r")
    lines = f.readlines()
    f.close()

    f = write_header(stepx, stepy, .1, xmin, ymin, -.1)

    for ix, line in enumerate(lines):
        if line[0] == "#":
            continue
        else:
            data = line.split()
            x = data[0]
            y = data[1]
            TLocal = float(data[2])
            ux = data[3]
            uy = data[4]
            pixx = data[5]
            pixy = data[6]
            piyy = data[7]
            pizz = data[8]

            eps = 3.0*cpLoc*TLocal*TLocal*TLocal*TLocal

            f.write(f"{x} {y} 0 {eps} 0 0 0 {ux} {uy} 0 0 {pixx} {pixy} 0 {piyy} 0 {pizz}\n")
    f.close()

if __name__ == "__main__":
    main()