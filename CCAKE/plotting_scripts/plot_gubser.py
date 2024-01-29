#!/usr/bin/env python
import sys
import os
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

time_list=['1.', '1.2', '1.4',]# '1.7', '2.', '2.2' ]
dpi=300
fig, ax = plt.subplots(2,3,figsize=np.array([1920,1920*2/3])/dpi, sharex=True,gridspec_kw={'wspace':.5})

def read_sol(analytic_sol_folder):
    for t in time_list:
        inp_path = os.path.join(analytic_sol_folder,'y=0_tau='+t+'_SemiAnalytic.dat')
        df = pd.read_table(inp_path,names=
                      ['x', 'y', 'T', 'ux', 'uy', 'Pixx', 'Piyy', 'Pixy', 'Pizz'], sep="     ", engine='python')
        ax[0][0].plot(df['x'],df['T'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[0][1].plot(df['x'],df['ux'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[0][2].plot(df['x'],df['uy'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[1][0].plot(df['x'],df['Pixx'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[1][1].plot(df['x'],df['Pixy'],label=r'$\tau = '+t+r'$ fm/c',color='blue')
        ax[1][2].plot(df['x'],df['Pizz']/float(t)**2,label=r'$\tau = '+t+r'$ fm/c',color='blue')
        

def read_sim(sim_result_folder):
    dt=.001
    for t in time_list:
        col_names=['id','t','x', 'y', 'p','T','muB','muS','muQ','e',
                   'rhoB','rhoS','rhoQ','s','s_smoothed','s_specific',
                   'sigma','spec_s','stauRelax','bigTheta','??',
                   '??2','pi00','pi11','pi22','pi12','t2pi33','v1','v2',
                   'gamma','frz','eos']
        idx = int(np.round((float(t)-1)/dt)/10)
        inp_path = os.path.join(sim_result_folder,f'system_state_{idx}.dat')
        df = pd.read_table(inp_path,
                           names=col_names,sep=' ',header=0)
        df['u1'] = df.loc[:,'v1']*df.loc[:,'gamma']
        df['u2'] = df.loc[:,'v2']*df.loc[:,'gamma']
        df_query = df.query(f"abs(y) < 1.5e-2")
        ax[0][0].plot(df_query['x'],df_query['T']/1000,label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[0][1].plot(df_query['x'],df_query['u1'],label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[0][2].plot(df_query['x'],df_query['u2'],label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[1][0].plot(df_query['x'],df_query['pi11']*.197,label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[1][1].plot(df_query['x'],df_query['pi12']*.197,label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')
        ax[1][2].plot(df_query['x'],df_query['t2pi33']*.197/float(t)**2,label=r'$\tau = '+t+r'$ fm/c',color='red',ls='--')


    ax[1][0].set_xlabel("x (fm)")
    ax[1][1].set_xlabel("x (fm)")
    ax[1][2].set_xlabel("x (fm)")
    ax[0][0].set_ylabel(r"$T$ (GeV)")
    ax[0][1].set_ylabel(r"$u^x$")
    ax[0][2].set_ylabel(r"$u^y$")
    ax[1][0].set_ylabel(r"$\pi^{xx}$")
    ax[1][1].set_ylabel(r"$\pi^{xy}$")
    ax[1][2].set_ylabel(r"$\pi^{\eta \eta}$")
    ax[0][0].set_xlim(-4.8,4.8)
    ax[0][1].set_xlim(-4.8,4.8)
    ax[0][2].set_xlim(-4.8,4.8)
    
    ax[0][0].set_ylim(0,.25)
    ax[0][1].set_ylim(-2.25,2.25)
    ax[0][2].set_ylim(-5E-13,5E-13)

    ax[1][0].set_ylim(-.45,.025)
    ax[1][1].set_ylim(-5E-15,5E-15)
    ax[1][2].set_ylim(-.025,.45)

if __name__ == '__main__':
    analytical_folder = sys.argv[1]
    simulation_folder = sys.argv[2]
    read_sol(analytical_folder)
    read_sim(simulation_folder)
    if (len(sys.argv) < 4):
        fig.savefig('test.png',dpi=300)
    else:
        fig.savefig(sys.argv[3],dpi=300)
