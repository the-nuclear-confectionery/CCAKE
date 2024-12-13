
import numpy as np
import sys

def init_cond(tau):
    stepeta = 0.01
    stepx = 0.1
    stepy = 0.1
    xmax = 0.1
    xmin = -xmax
    ymax = 0.1
    ymin = -ymax
    etamax = 10
    etamin = -etamax
    eps0 = 10
    etas = np.arange(-10,10.01,stepeta).round(3)
    #etas = np.array([num + 0.005 for num in etass])
    c_sq = 1/3
    a = 1
    tau0 = 1
    t0 = 0.01 * tau0
    c1 = (1 - c_sq**2) / (4 * c_sq)
    c2 = (1 + c_sq)**2 / (4 * c_sq)
    c = c1/a**2 - c2
    hbarc_G = 0.1973269804
    fname=f"../long_data/ic_long_tau_{tau}.dat"
    with open(fname, "w") as f:
        f.write(f"#0 {stepx} {stepy} {stepeta} 0 {xmin} {ymin} {etamin}\n")
        for eta in etas:
            eps = eps0 * (t0/tau0 + a*tau/tau0 * np.exp(eta))**c * (t0/tau0 + tau/(a*tau0) * np.exp(-eta))**c
            u = 0.5 / tau * (np.sqrt((t0*np.exp(-eta)+tau*a) / (t0*np.exp(eta)+(tau/a))) - np.sqrt((t0*np.exp(eta)+(tau/a)) / (t0*np.exp(-eta)+tau*a)))
            #eta_shift = eta + 2.0
            f.write(f"0 0 {eta} {eps} 0 0 0 0 0 {u} 0 0 0 0 0 0 0\n")
        f.close()
    print(f'Initial conditions were successfully generated for tau = {tau}')


def main():
    if len(sys.argv) != 2:
        print("No argument provided. Please provide a value for tau in the following format: "+sys.argv[0]+" 1.0")
    else:
        argument = float(sys.argv[1])
        init_cond(argument)

if __name__ == "__main__":
    main()
