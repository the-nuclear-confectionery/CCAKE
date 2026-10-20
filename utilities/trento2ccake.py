#!/usr/bin/env python3
import sys
import os

def check_input():
    if len(sys.argv) != 5:
        print("Usage: trento2ccake.py <input file> <spacing> <xmax> <output>")
        sys.exit(1)

    if not os.path.isfile(sys.argv[1]):
        print(f"File {sys.argv[1]} does not exist")
        sys.exit(1)

def write_header(outpath, stepx, stepy, stepEta, xmin, ymin, etamin):
    f = open(outpath, "w")
    # Keep your original header format
    f.write(f"#0 {stepx} {stepy} {stepEta} 0 {xmin} {ymin} {etamin}\n")
    return f

def read_trento_block(path):
    """
    Reads regular TRENTo write_text_file output:
      - header lines start with '#'
      - then each line is a row of floats (grid)
    Returns grid (list of rows), nx, ny
    """
    rows = []
    with open(path, "r") as infile:
        for line in infile:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                continue
            vals = [float(x) for x in line.split()]
            if vals:
                rows.append(vals)

    if not rows:
        raise RuntimeError("No grid rows found (after skipping '#').")

    ny = len(rows[0])
    for i, r in enumerate(rows):
        if len(r) != ny:
            raise RuntimeError(f"Inconsistent row length at row {i}: got {len(r)}, expected {ny}")

    nx = len(rows)
    return rows, nx, ny

def main():
    check_input()

    inpath = sys.argv[1]
    spacing = float(sys.argv[2])
    xmax = float(sys.argv[3])
    outpath = sys.argv[4]

    stepx = spacing
    stepy = spacing
    stepEta = 0.1
    xmin = -xmax
    ymin = -xmax
    etamin = -0.5

    grid, nx, ny = read_trento_block(inpath)

    # Sanity: grid extent implied by nx,ny and spacing
    implied_Lx = (nx - 1) * spacing
    implied_Ly = (ny - 1) * spacing
    implied_xmax = 0.5 * implied_Lx
    implied_ymax = 0.5 * implied_Ly

    # If user's xmax doesn't match implied xmax, we still follow the user's xmin/ymin,
    # but warn because it means coordinates won't line up with the grid spacing/count.
    if abs(xmax - implied_xmax) > 1e-6:
        print(f"WARNING: provided xmax={xmax} but grid implies xmax~{implied_xmax:.6f} "
              f"(nx={nx}, spacing={spacing}). Using provided xmax for coordinates.")

    outfile = write_header(outpath, stepx, stepy, stepEta, xmin, ymin, etamin)

    # Write one line per cell, reconstructing x,y from indices
    # We assume: row index -> x direction, column index -> y direction
    for ix in range(nx):
        x = xmin + ix * spacing
        row = grid[ix]
        for iy in range(ny):
            y = ymin + iy * spacing
            eps = row[iy]*75

            # Match your original 17-column output format:
            # x y eta eps ux uy ueta ... (all zeros)
            outfile.write(f"{x} {y} 0 {eps} 0 0 0 0 0 0 0 0 0 0 0 0 0\n")

    outfile.close()

    print("=== done ===")
    print(f"input : {inpath}")
    print(f"output: {outpath}")
    print(f"nx={nx} ny={ny} spacing={spacing}")
    print(f"xmin={xmin} ymin={ymin} etamin={etamin}")

if __name__ == "__main__":
    main()