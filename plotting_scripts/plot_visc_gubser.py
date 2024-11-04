import sys
import os
import matplotlib as mpl
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from variable_conversions import HBARC  # noqa
import my_plotting as myplt  # noqa


plot_index = 1
analytic_style = {'ls': '-', 'lw': 2}
# sim_style = {'ls':'-.','lw':3.,'alpha':.5}
# sim_style = {'facecolors': 'none'}

time_list = np.arange(1.00, 1.50, 0.1)  # Use this to focus on before FO
# time_list=['1.00', '1.20',  '1.40',  '1.60']# Use this to focus on before FO
# time_list=['1.00', '1.50', '2.00', '3.00', '4.00', '5.00', '6.00', '7.00', '8.00', '9.00', '10.00' ]
filter_criteria = 'abs(phi - 3.141592653589793/4.) < 1.e-2'

cmap = myplt.get_cmap(len(time_list), 'cividis')

mpl.rcParams['text.usetex'] = True

dpi = 150
fig, ax = plt.subplot_mosaic([['e', 'rhoB', 'ux', 'cbar'],
                              ['pixx', 'Rey', 'pixy', 'cbar']],
                             width_ratios=[1, 1, 1, 0.1],
                             figsize=np.array([7 * 3, 7 * 2]),
                             constrained_layout=True)


def get_reynolds_number(df, t_squared):
    df['u0'] = np.sqrt(1 + df['ux']**2 + df['uy']**2)
    df['pitx'] = (df['pixx'] * df['ux'] + df['pixy'] * df['uy']) / df['u0']
    df['pity'] = (df['pixy'] * df['ux'] + df['piyy'] * df['uy']) / df['u0']
    df['pitt'] = (df['pitx'] * df['ux'] + df['pity'] * df['uy']) / df['u0']
    df['pizz'] = -(df['pixx'] + df['piyy'] - df['pitt']) / t_squared
    df['p'] = df['e'] / 3  # Pressure

    df['pi_norm'] = np.sqrt(df['pitt']**2 + df['pixx']**2 + df['piyy']**2
                            + (t_squared * df['pizz'])**2
                            + 2 * (df['pixy']**2 - df['pitx']**2
                                   - df['pity']**2))
    df['reynolds'] = df['pi_norm'] / df['p']
    return df


def read_sol(analytic_sol_folder):
    for ii, t in enumerate(time_list):
        if ii % plot_index != 0:
            continue
        inp_path = os.path.join(analytic_sol_folder, f'init_conditions.txt')
        df = pd.read_table(
            inp_path,
            names=[
                'x',
                'y',
                'eta',
                'e',
                'rhoB',
                'rhoS',
                'rhoQ',
                'ux',
                'uy',
                'ueta',
                'Bulk',
                'pixx',
                'pixy',
                'pixeta',
                'piyy',
                'pyeta',
                'pietaeta'],
            sep=" ",
            engine='python',
            header=1)

        df['r'] = np.sqrt(df['x']**2 + df['y']**2)
        df['phi'] = np.arctan2(df['y'], df['x'])
        df = get_reynolds_number(df, float(t)**2)

        df = df.query(filter_criteria)

        ax['e'].plot(df['r'], df['e'], color=cmap(ii), **analytic_style)
        ax['rhoB'].plot(df['r'], df['rhoB'], color=cmap(ii), **analytic_style)
        ax['ux'].plot(df['r'], df['ux'], color=cmap(ii), **analytic_style)
        ax['pixx'].plot(df['r'], df['pixx'], color=cmap(ii), **analytic_style)
        ax['Rey'].plot(
            df['r'],
            df['reynolds'],
            color=cmap(ii),
            **analytic_style)
        ax['pixy'].plot(
            df['r'],
            df['pixy'],
            color=cmap(ii),
            **analytic_style)


def read_sim(sim_result_folder):
    dt = .001
    for ii, t in enumerate(time_list):
        if ii % plot_index != 0:
            continue
        col_names = [
            'id',
            't',
            'x',
            'y',
            'p',
            'T',
            'muB',
            'muS',
            'muQ',
            'e',
            'rhoB',
            'rhoS',
            'rhoQ',
            's',
            's_smoothed',
            's_specific',
            'sigma',
            'spec_s',
            'stauRelax',
            'bigTheta',
            '??',
            '??2',
            'pi00',
            'pixx',
            'piyy',
            'pixy',
            't2pi33',
            'v1',
            'v2',
            'gamma',
            'frz',
            'eos']
        idx = int(np.round((float(t) - 1) / dt) / 100)
        inp_path = os.path.join(sim_result_folder, f'system_state_{idx}.dat')
        print(inp_path)
        df = pd.read_table(inp_path,
                           names=col_names, sep=' ', header=0)
        df['ux'] = df.loc[:, 'v1'] * df.loc[:, 'gamma']
        df['uy'] = df.loc[:, 'v2'] * df.loc[:, 'gamma']
        df['r'] = np.sqrt(df['x']**2 + df['y']**2)
        df['phi'] = np.arctan2(df['y'], df['x'])
        df['e'] = df['e'] / 1000  # convert to GeV/fm^3

        df['pixx'] = df['pixx'] * HBARC  # convert to GeV/fm^3
        df['pixy'] = df['pixy'] * HBARC  # convert to GeV/fm^3
        df['piyy'] = df['piyy'] * HBARC  # convert to GeV/fm^3

        df['t2pi33'] = df['t2pi33'] * HBARC  # convert to GeV/fm^3

        df = get_reynolds_number(df, float(t)**2)
        df_query = df.query(filter_criteria)

        sim_style = {'facecolor': cmap(ii), 's': 60, 'edgecolors': 'k'}

        stride = 1
        ax['e'].scatter(
            df_query['r'].to_numpy()[
                ::stride], df_query['e'].to_numpy()[
                ::stride], **sim_style)
        ax['rhoB'].scatter(
            df_query['r'].to_numpy()[
                ::stride], df_query['rhoB'].to_numpy()[
                ::stride], **sim_style)
        ax['ux'].scatter(
            df_query['r'].to_numpy()[
                ::stride], df_query['ux'].to_numpy()[
                ::stride], **sim_style)
        ax['pixx'].scatter(
            df_query['r'].to_numpy()[
                ::stride], df_query['pixx'].to_numpy()[
                ::stride], **sim_style)
        ax['Rey'].scatter(
            df_query['r'].to_numpy()[
                ::stride], df_query['reynolds'].to_numpy()[
                ::stride], **sim_style)
        ax['pixy'].scatter(
            df_query['r'].to_numpy()[
                ::stride], df_query['pixy'].to_numpy()[
                ::stride],  **sim_style)


def beautify():
    # fig.set_tight_layout(True)
    ylabels = {'e': r'$\mathcal E$ (GeV/fm$^3$)',
               'rhoB': r'$n_B$ (fm$^{-3}$)',
               'ux': r'$u^x$',
               'pixx': r'$\pi^{xx}$ (GeV/fm$^3$)',
               'Rey': r'$\mathcal{R}^{-1}$',
               # 'pixy':r'$\pi^{xy}$ (GeV/fm$^3$)',
               'pixy': r'$\pi^{xy}$ (GeV/fm$^3$)'}
    for key in ax.keys():
        if key == 'cbar':
            continue
        myplt.costumize_axis(ax=ax[key],
                             x_title=r'$r$ (fm)',
                             y_title=ylabels[key])
        ax[key].set_xlim(0, 4.5)

    myplt.costumize_axis(ax=ax['cbar'], x_title='', y_title='')

    tau_min = float(time_list[0])
    tau_max = float(time_list[-1])
    tau_list = [float(t) for t in time_list]

    # deal with the colorbar
    delta_time = tau_list[1] - tau_list[0]
    boundaries = np.arange(tau_min - delta_time / 2,
                           tau_max + 3 * delta_time / 2, delta_time)

    norm = mpl.colors.Normalize(vmin=tau_min, vmax=tau_max)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks=tau_list, label=r'$\tau$ (fm/c)',
                        cax=ax['cbar'],
                        boundaries=boundaries)
    ax['e'].plot([], [], **analytic_style, label='Analytic', color=cmap(0))
    ax['e'].scatter(
        [],
        [],
        label='CCAKE',
        edgecolors=cmap(0))
    ax['e'].text(
        0.90,
        0.85,
        "EoS 2",
        transform=ax['e'].transAxes,
        fontsize=20,
        bbox={'facecolor': 'white'},
        horizontalalignment='center',
    )

    # ax['ux'].set_ylim(0,3.3)
    # ax['ux'].set_ylim(0, 1.3)
    # ax['e'].set_ylim(0, 8)
    # ax['Rey'].set_ylim(0, 3.1)

    for name, label in zip(ylabels.keys(),
                           ['a', 'b', 'c', 'd', 'e', 'f']):
        ax[name].text(
            0.93,
            0.93,
            f'({label})',
            transform=ax[name].transAxes,
            fontsize=20,
            bbox={'boxstyle': 'round', 'facecolor': 'white'},
            horizontalalignment='center',
        )

    ax['e'].legend(loc='upper center', fontsize=18)


if __name__ == '__main__':
    analytical_folder = sys.argv[1]
    simulation_folder = sys.argv[2]
    read_sol(analytical_folder)
    read_sim(simulation_folder)
    beautify()

    if (len(sys.argv) < 4):
        fig.savefig('gubser_without_regulator.pdf')
    else:
        fig.savefig(sys.argv[3])
