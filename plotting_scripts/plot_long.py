#!/usr/bin/env python
import sys
import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

time_list=['1.0', '1.5', '2.0'] #, '1.2', '1.4',] '1.7', '2.', '2.2' ]
dpi=300
fig, ax = plt.subplots(2,figsize=np.array([1920,1920*2/3])/dpi, sharex=True,gridspec_kw={'wspace':.5})

def read_sol(analytic_sol_folder):
    for t in time_list:
        inp_path = os.path.join(analytic_sol_folder,'ic_long_tau_'+t+'.dat')
        col_names=['x', 'y', 'eta', 'eps', 'rhoB', 'rhoS', 'rhoQ', 'ux', 'uy', 'ueta','Bulk', 'Pixx', 'Pixy', 'Pixeta' 'Piyy', 'Piyeta', 'Pietaeta','extra']
        df = pd.read_table(inp_path, names=col_names, sep=" ", header=0)
        ax[0].plot(df['eta'],df['ueta'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[1].plot(df['eta'],df['eps'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[0].legend()
        ax[1].legend()

def read_sim(sim_result_folder):
    dt=.05
    for t in time_list:
        col_names=['id','tau', 'eta', 'p', 'T', 'muB', 'muS', 'muQ', 'eps', 'rhoB',
                   'rhoS', 'rhoQ', 's', 'smoothed_s', 'specific_s', 'sigma',
                   'norm_spec_s', 'stauRelax', 'bigtheta', "??", "??2",
                   'shv00', 'shv11', 'shv22', 'hydro_shv12', 't^2 shv33', 
                   'v_eta', 'gamma', 'freeze', 'eos_name']
        idx = int(np.round((float(t)-1)/dt)/10)
        inp_path = os.path.join(sim_result_folder,f'system_state_{idx}.dat')
        df = pd.read_table(inp_path,
                           names=col_names,sep=' ',header=1)
        df['ueta'] = df.loc[:,'v_eta']*df.loc[:,'gamma']
        ax[0].plot(df['eta'],df['ueta'],label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[1].plot(df['eta'],df['eps'],label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[0].legend()
        ax[1].legend()
    
    ax[0].set_ylabel(r"$u^{eta}$")
    ax[1].set_ylabel(r"$eps$")
    ax[0].set_xlabel("eta")
    ax[1].set_xlabel("eta")

#    ax[0].set_xlim(-4.8,4.8)
#    ax[1].set_xlim(-4.8,4.8)
#    ax[0].set_ylim(-5E-13,5E-13)


if __name__ == '__main__':
    analytical_folder = sys.argv[1]
    simulation_folder = sys.argv[2]
    read_sol(analytical_folder)
    read_sim(simulation_folder)
    if (len(sys.argv) < 4):
        fig.savefig('test.png',dpi=300)
    else:
        fig.savefig(sys.argv[3],dpi=300)
