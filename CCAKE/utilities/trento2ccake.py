#!/usr/bin/env python

import sys
import os

def check_input():
    if len(sys.argv) != 5:
        print("Usage: trento2ccake.py <input file> <spacing> <xmax> <output>")
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print("File %s does not exist" % sys.argv[1])
        sys.exit(1)

def write_header(stepx, stepy, stepEta, xmin, ymin, etamin):
    #TODO: parse from a trento input file
    f = open(sys.argv[4], "w")
    f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
    return f

def main():
    check_input()
    f = open(sys.argv[1], "r")
    lines = f.readlines()
    f.close()

    stepx  = float(sys.argv[2])
    stepy  = float(sys.argv[2])
    stepEta  = .1
    xmin = -float(sys.argv[3])
    ymin = -float(sys.argv[3])
    etamin = -.5

    f = write_header(stepx, stepy, stepEta, xmin, ymin, etamin)

    for ix, line in enumerate(lines):
        if line[0] == "#":
            continue
        else:
            x = xmin + ix*stepx
            data = line.split()
            for iy, eps in enumerate(data):
                y = ymin + iy*stepy
                f.write(f"{x} {y} 0 {eps} 0 0 0 0 0 0 0 0 0 0 0 0 0\n")
    f.close()

if __name__ == "__main__":
    main()