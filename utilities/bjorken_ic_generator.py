import numpy as np
import pandas as pd
import sys

def init_cond(tau):
    eps0 = 1000
    etas = np.arange(-10,10.1,0.1).round(3)
    c_sq = 1/3
    a = 1
    tau0 = 1
    t0 = 0.01 * tau0
    c1 = (1 - c_sq**2) / (4 * c_sq)
    c2 = (1 + c_sq)**2 / (4 * c_sq)
    u_eta = 0
    u_tau = 0
    file_path = f"../bjorken/bjorken_tau_{tau}.dat"
    with open(file_path, "w") as f:
        for eta in etas:
            eps = eps0 * (tau0 / tau)**(4/3) 
#             eps = eps0 * (tau0 / tau)**2
            u_eta = 0 #np.sinh(eta)
            u_tau = 1 #np.cosh(eta) 
            f.write(f"0 0 {eta} {eps} 0 0 0 0 0 {u_eta} 0 0 0 0 0 0 0\n")
        f.close()
        f.close()

def main():
    if len(sys.argv) != 2:
        print("No argument provided. Please provide a value for tau in the following format: "+sys.argv[0]+" 1.0")
    else:
        argument = float(sys.argv[1])
        init_cond(argument)

if __name__ == "__main__":
    main()
