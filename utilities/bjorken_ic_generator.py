import numpy as np
import pandas as pd
import sys

def init_cond(tau):
    stepeta = 0.1
    stepx = 0.1
    stepy = 0.1
    xmax = 0.1
    xmin = -xmax
    ymax = 0.1
    ymin = -ymax
    etamax = 10
    etamin = -etamax
    hbarc_G = 0.1973269804
    eps0 = 10 / hbarc_G**3
    etas = np.arange(-10,10.1,0.1).round(3)
    tau0 = 1
    u_eta = 0
    u_tau = 1
    fp = f"../bjorken/bjorken_tau_{tau}.dat"
    with open(fp, "w") as f:
        f.write(f"#0 {stepx} {stepy} {stepeta} 0 {xmin} {ymin} {etamin}\n")
        for eta in etas:
            eps = eps0 * (tau0 / tau)**(4/3)
            u_eta = 5
            f.write(f"0 0 {eta} {eps} 0 0 0 0 0 {u_eta} 0 0 0 0 0 0 0\n")
        print(f'Initial conditions were successfully generated for tau = {tau}')
        f.close()

def main():
    if len(sys.argv) != 2:
        print("No argument provided. Please provide a value for tau in the following format: "+sys.argv[0]+" 1.0")
    else:
        argument = float(sys.argv[1])
        init_cond(argument)

if __name__ == "__main__":
    main()
