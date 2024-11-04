from scipy.integrate import odeint
from scipy.interpolate import interp1d

from numpy import linspace
from numpy import ndarray
from numpy import array
from numpy import concatenate

import matplotlib.pyplot as plt

from my_plotting import costumize_axis

from typing import List

import sys
sys.path.append('../utilities')
from system_conformal_plasma import ConformalPlasma
from variable_conversions import HBARC
from variable_conversions import milne_T
from variable_conversions import milne_mu

T_PLOT = (0, 0)
MU_PLOT = (0, 1,)
PIXX_PLOT = (1, 0)
PIXY_PLOT = (1, 1)

E_PLOT = (0)
N_PLOT = (1,)
S_PLOT = (2,)

CONST_T0 = 1.0
CONST_MU0 = array([1.0, 1.0, 1.0]).reshape(-1, 1)

system = ConformalPlasma(temperature_0=CONST_T0, chem_potential_0=CONST_MU0)

def solve_and_plot(
        ax_1: plt.Axes,
        ax_2: plt.Axes,
        y0s: ndarray,
        rhos_1: ndarray,
        rhos_2: ndarray,
        xs: ndarray,
        taus: ndarray,
        color: List[str],
        linestyle: List[str],
        add_labels: bool = False,
) -> None:
    soln_1 = odeint(system.eom.for_scipy, y0s, rhos_1,)
    soln_2 = odeint(system.eom.for_scipy, y0s, rhos_2,)
    t_hat = concatenate((soln_1[:, 0][::-1], soln_2[:, 0]))
    mu_hat = [concatenate((soln_1[:, i][::-1], soln_2[:, i]))
              for i in [1, 2, 3]]
    pi_bar_hat = concatenate((soln_1[:, 4][::-1], soln_2[:, 4]))
    rhos = concatenate((rhos_1[::-1], rhos_2))

    t_interp = interp1d(rhos, t_hat)
    mu_interp = [interp1d(rhos, f) for f in mu_hat]
    pi_interp = interp1d(rhos, pi_bar_hat)

    for n, tau in enumerate(taus):
        t_evol = milne_T(tau, xs, 1, t_interp)
        mu_evol = milne_mu(tau, xs, 1, mu_interp[0])

        e_evol = system.milne_energy(tau, xs, 0.0, 1.0, t_interp, mu_interp)
        n_evol = system.milne_number(tau, xs, 0.0, 1.0, t_interp, mu_interp)[0]
        s_evol = system.milne_entropy(tau, xs, 0.0, 1.0, t_interp, mu_interp)

        ax_1[T_PLOT].plot(xs, t_evol,
                          color=color[n], lw=2, ls=linestyle[n],
                          label=r'$\mu_0/T_0=' + f'{y0s[1]/y0s[0]:.2f}$'
                          if n == 0 else None)
        ax_1[MU_PLOT].plot(xs, mu_evol,
                           color=color[n], lw=2, ls=linestyle[n],
                           label=r'$\tau = ' + f'{tau:.2f}' + r'$ [fm/$c$]'
                           if add_labels else None)

        pi_xx, pi_yy, pi_xy, pi_nn = system.milne_pi(
            tau,
            xs,
            0.0,
            1,
            t_interp,
            mu_interp,
            pi_interp,
            nonzero_xy=True,
        )

        ax_1[PIXX_PLOT].plot(xs,
                             pi_yy / (4.0 * e_evol / 3.0),
                             color=color[n], lw=2, ls=linestyle[n])

        ax_1[PIXY_PLOT].plot(xs,
                             pi_xy / (4.0 * e_evol / 3.0),
                             color=color[n], lw=2, ls=linestyle[n])

        ax_2[E_PLOT].plot(xs, e_evol,
                          color=color[n], lw=2, ls=linestyle[n],
                          label=r'$\mu_0/T_0=' + f'{y0s[1]/y0s[0]:.2f}$'
                          if n == 0 else None)
        ax_2[N_PLOT].plot(xs, n_evol,
                          color=color[n], lw=2, ls=linestyle[n],
                          label=r'$\tau=' + f'{tau:.2f}$ [fm/$c$]'
                          if add_labels else None)
        ax_2[S_PLOT].plot(xs, s_evol,
                          color=color[n], lw=2, ls=linestyle[n])


def main():
    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=(2 * 7, 2 * 7))
    fig.patch.set_facecolor('white')

    fig2, ax2 = plt.subplots(ncols=3, nrows=1, figsize=(3 * 7, 1 * 7))

    rhos_1 = linspace(-10, 0, 1000)[::-1]
    rhos_2 = linspace(0, 10, 1000)
    taus = array([1.2, 2.0, 3.0])
    xs = linspace(-6, 6, 200)

    t0 = 0.25 * 1.2 / HBARC
    colors = ['black', 'red', 'blue']
    for n, alpha in enumerate([0.025]):
        y0s = array([
            t0,
            *[alpha * t0 for _ in range(3)],
            0.0
        ])
        solve_and_plot(
            ax_1=ax,
            ax_2=ax2,
            y0s=y0s,
            rhos_1=rhos_1,
            rhos_2=rhos_2,
            xs=xs,
            taus=taus,
            color=3 * [colors[n]],
            linestyle=['solid', 'dashed', 'dotted'],
            add_labels=True if n == 0 else False,
        )

    costumize_axis(
        ax=ax[T_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$T(\tau, x)$ [GeV]'
    )
    costumize_axis(
        ax=ax[MU_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$\mu(\tau, x)$ [GeV]'
    )
    costumize_axis(
        ax=ax[PIXX_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$\pi^{yy}(\tau, x) / w(\tau, x)$'
    )
    costumize_axis(
        ax=ax[PIXY_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$\pi^{xy}(\tau, x) / w(\tau, x)$'
    )

    ax[T_PLOT].legend(loc='upper right', fontsize=20)
    ax[MU_PLOT].legend(loc='upper right', fontsize=20)
    fig.tight_layout()
    fig_name = './viscous-gubser-BSQ-CCAKE-1.pdf'
    print(f'saving figure to {fig_name}')
    fig.savefig(fig_name)

    costumize_axis(
        ax=ax2[E_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$\mathcal E(\tau, x)$ [GeV/fm$^{-3}$]'
    )
    ax2[E_PLOT].set_yscale('log')
    costumize_axis(
        ax=ax2[N_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$n(\tau, x)$ [fm$^{-3}$]'
    )
    ax2[N_PLOT].set_yscale('log')
    ax2[N_PLOT].set_ylim(bottom=1e-2)
    costumize_axis(
        ax=ax2[S_PLOT],
        x_title=r'$x$ [fm]',
        y_title=r'$s(\tau, x)$ [fm$^{-3}$]'
    )
    ax2[S_PLOT].set_yscale('log')

    ax2[E_PLOT].legend(loc='upper right', fontsize=20)
    ax2[N_PLOT].legend(loc='upper right', fontsize=20)
    fig2.tight_layout()
    fig_name = './viscous-gubser-BSQ-CCAKE-2.pdf'
    print(f'saving figure to {fig_name}')
    fig2.savefig(fig_name)


if __name__ == "__main__":
    main()
