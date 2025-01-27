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
    # Open the output file and write the header
    f = open(sys.argv[4], "w")
    f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
    return f

def main():
    check_input()
    
    # Read the input file and skip the first line
    with open(sys.argv[1], "r") as infile:
        lines = infile.readlines()[1:]

    stepx = float(sys.argv[2])
    stepy = float(sys.argv[2])
    stepEta = 0.1
    xmin = -float(sys.argv[3])
    ymin = -float(sys.argv[3])
    etamin = -0.5

    # Write the header to the output file
    outfile = write_header(stepx, stepy, stepEta, xmin, ymin, etamin)

    for line in lines:
        # Skip comment lines
        if line.startswith("#"):
            continue

        # Parse x, y, and eps from the line
        x, y, eps = map(float, line.split())

        # Write the formatted output
        outfile.write(f"{x} {y} 0 {eps} 0 0 0 0 0 0 0 0 0 0 0 0 0\n")

    outfile.close()

if __name__ == "__main__":
    main()
