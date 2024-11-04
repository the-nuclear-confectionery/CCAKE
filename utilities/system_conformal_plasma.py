
from numpy import ndarray
from numpy import pi
from numpy import array
from numpy import sum as npsum
from numpy import sinh
from numpy import tanh
from numpy import arcsinh
from numpy import arctanh
from numpy import sqrt
from numpy import fabs
from numpy import where

from scipy.interpolate import interp1d

from typing import Callable
from typing import List
from typing import Union

import variable_conversions as conv

# Global constants

ALPHA = 3 * pi ** 2 * (16 + 105 / 4) / 90   # degeneracy of QGP
CTAUR = 5                                   # relaxation time constant
ETA_S = 0.2                                 # shear viscosity to entory ratio

# Equation of State class
class EoS:
    '''
    Corresponds to EoS 1 in the paper and describes a masslass quark gluon
    plasma
    '''
    def __init__(
            self,
            temperature_0: float,
            chem_potential_0: ndarray,
    ):
        r'''
        parameters:
        ----------
        temperature_0 float
            The constant $T_\ast$ in that sets the scale of the temperature
        chem_potential_0 ndarray
            The constants $\mu_\ast$ for baryon, strange and electric charge
            chemical potentials that sets their relative size
        '''
        self.temperature_0 = temperature_0
        self.chem_potential_0 = chem_potential_0

    def pressure(
            self,
            temperature: Union[float, ndarray],
            chem_potential: Union[float, ndarray],
    ) -> Union[float, ndarray]:
        '''
        parametes:
        ----------
        temperature Union[float, np.ndarray]
        chem_potential Union[float, np.ndarray]
            Assumed to be a single baryon chemical potential

        returns the pressure
        '''
        return_value = (temperature / self.temperature_0) ** 2
        return_value += npsum((chem_potential / self.chem_potential_0) ** 2)
        return_value = self.temperature_0 ** 4 * return_value ** 2
        return ALPHA * return_value

    def energy(
            self,
            temperature: Union[float, ndarray],
            chem_potential: Union[float, ndarray],
    ) -> Union[float, ndarray]:
        '''
        parametes:
        ----------
        temperature Union[float, np.ndarray]
        chem_potential Union[float, np.ndarray]
            Assumed to be a single baryon chemical potential

        returns the energy density
        '''
        return_value = (temperature / self.temperature_0) ** 2
        return_value += npsum((array(chem_potential) / self.chem_potential_0) ** 2)
        return_value = self.temperature_0 ** 4 * return_value ** 2
        return 3 * ALPHA * return_value

    def entropy(
            self,
            temperature: Union[float, ndarray],
            chem_potential: Union[float, ndarray],
    ) -> Union[float, ndarray]:
        '''
        parametes:
        ----------
        temperature Union[float, np.ndarray]
        chem_potential Union[float, np.ndarray]
            Assumed to be a single baryon chemical potential

        returns the entropy density
        '''
        temp_value = (temperature / self.temperature_0) ** 2
        temp_value += npsum((chem_potential / self.chem_potential_0) ** 2)
        return_value = 4 * temperature * temp_value \
            * self.temperature_0 ** 2
        return ALPHA * return_value

    def number(
            self,
            temperature: Union[float, ndarray],
            chem_potential: Union[float, ndarray],
    ) -> Union[float, ndarray]:
        '''
        parametes:
        ----------
        temperature Union[float, np.ndarray]
        chem_potential Union[float, np.ndarray]
            Assumed to be a single baryon chemical potential

        returns the number density
        '''
        # Note that we assume the order is [B, S, Q]
        temp_value = (temperature / self.temperature_0) ** 2
        temp_value += npsum((chem_potential / self.chem_potential_0) ** 2)
        return_value = 4 * self.temperature_0 ** 4 * chem_potential
        return_value *= temp_value / self.chem_potential_0 ** 2
        return ALPHA * return_value

# The equations of motion class
class EoM:
    '''
    This class implementst the equations of motion for temperature, chemcical
    potentails and shear pressure for the viscous Gubser with BSQ currents
    '''
    def __init__(
            self,
            eos: EoS,
            temperature_0: float,
            chem_potential_0: ndarray,
    ):
        r'''
        parameters:
        ----------
        tau_R Callable
            A function that returns the relaxation time
        temperature_0 float
            The constant $T_\ast$ in that sets the scale of the temperature
        chem_potential_0 ndarray
            The constants $\mu_\ast$ for baryon, strange and electric charge
            chemical potentials that sets their relative size
        '''
        self.eos = eos
        self.temperature_0 = temperature_0
        self.chem_potential_0 = chem_potential_0

    def tau_R(
            self,
            temperature: float,
            chem_potential: ndarray,
    ) -> float:
        r'''
        The relaxation time calculate by calling the equation of state with the
        natural equation of state variables

        parameters:
        ----------
        temperature float
        chem_potential float
        temperature_0 float
            The constant $T_\ast$ in that sets the scale of the temperature
        chem_potential_0 ndarray
            The constants $\mu_\ast$ for baryon, strange and electric charge
            chemical potentials that sets their relative size

        return the relaxation time
        '''
        e = self.eos.energy(
                temperature=temperature, chem_potential=chem_potential)
        p = self.eos.pressure(
            temperature=temperature, chem_potential=chem_potential)
        s = self.eos.entropy(
                temperature=temperature, chem_potential=chem_potential)
        return CTAUR * ETA_S * s / (e + p)

    def dT_drho(
            self,
            ys: ndarray,
            rho: Union[float, ndarray],
    ) -> Union[float, ndarray]:
        '''
        Derivative of temperature with respect to rho (de Sitter time)

        parameters:
        -----------
        ys ndarray
            The quantities being evolved.
            Expected to follow the order:
                [temperature, chemical potential 1, chemcical potential 2, ...]
        rho Union[float, ndarray]
            de Sitter time for time stepping

        returns dT/drho evaluated at time rho
        '''
        temperature, mu_B, mu_S, mu_Q, pi_hat = ys
        chem_potential = array([mu_B, mu_S, mu_Q])
        return_value = npsum(pi_hat * (chem_potential * self.temperature_0
                                       / (self.chem_potential_0 * temperature)) ** 2)
        return_value += (1 / 3) * (-2 + pi_hat)
        return temperature * return_value * tanh(rho)

    def dmu_drho(
            self,
            ys: ndarray,
            rho: Union[float, ndarray],
    ) -> Union[float, ndarray]:
        '''
        Derivative of chemical potential(s) with respect to rho (de Sitter
        time)

        parameters:
        -----------
        ys ndarray
            The quantities being evolved.
            Expected to follow the order:
                [temperature, chemical potential 1, chemcical potential 2, ...]
        rho Union[float, ndarray]
            de Sitter time for time stepping

        returns dmu/drho evaluated at time rho
        '''
        _, mu_B, mu_S, mu_Q, pi_hat = ys
        chem_potential = array([mu_B, mu_S, mu_Q])
        return_value = 1 + pi_hat
        return -(2 / 3) * chem_potential * return_value * tanh(rho)

    def dpi_drho(
            self,
            ys: ndarray,
            rho: Union[float, ndarray]
    ) -> Union[float, ndarray]:
        '''
        Derivative of shear pressure with respect to rho (de Sitter time)

        parameters:
        -----------
        ys ndarray
            The quantities being evolved.
            Expected to follow the order:
                [temperature, chemical potential 1, chemcical potential 2, ...]
        rho Union[float, ndarray]
            de Sitter time for time stepping

        return dpi/drho evaluated at time rho
        '''
        temperature, mu_B, mu_S, mu_Q, pi_hat = ys
        chem_potential = array([mu_B, mu_S, mu_Q])
        tau_r = self.tau_R(temperature, chem_potential)
        return_value = (4 / 3 / CTAUR) * tanh(rho)
        return_value -= pi_hat / tau_r
        return_value -= (4 / 3) * pi_hat ** 2 * tanh(rho)
        return return_value

    def for_scipy(
            self,
            ys: ndarray,
            rho: Union[float, ndarray],
            ideal: bool = False
    ) -> ndarray:
        '''
        Packs equations of motion into single array for
        scipy.integrate.odeint to integrate

        parameters:
        -----------
        ys ndarray
            The quantities being evolved.
            Expected to follow the order:
                [temperature, chemical potential 1, chemcical potential 2, ...]
        rho Union[float, ndarray]
            de Sitter time for time stepping
        ideal bool
            If true, sets the viscous corrections to zero

        return dpi/drho evaluated at time rho
        '''
        dTdrho = self.dT_drho(ys, rho)
        dmudrho = self.dmu_drho(ys, rho)
        dpidrho = 0 if ideal else self.dpi_drho(ys, rho)
        return array([dTdrho, *dmudrho.reshape(-1,), dpidrho])


class ConformalPlasma:
    def __init__(
            self,
            temperature_0: float,
            chem_potential_0: ndarray,
        ):
        r'''
        parameters:
        ----------
        temperature_0 float
            The constant $T_\ast$ in that sets the scale of the temperature
        chem_potential_0 ndarray
            The constants $\mu_\ast$ for baryon, strange and electric charge
            chemical potentials that sets their relative size
        '''
        self.eos = EoS(temperature_0, chem_potential_0)
        self.eom = EoM(self.eos, temperature_0, chem_potential_0)
        self.temperature_0 = temperature_0
        self.chem_potential_0 = chem_potential_0


    def denergy_drho(
            self,
            ys: ndarray,
            rho: Union[float, ndarray],
    ) -> ndarray:
        '''
        Function used to help calculate freeze-out surface and normal vectors


        parameters:
        -----------
        ys ndarray
            The quantities being evolved.
            Expected to follow the order:
                [temperature, chemical potential 1, chemcical potential 2, ...]
        rho Union[float, ndarray]
            de Sitter time for time stepping

        returns d(energy density)/drho
        '''
        temperature, mu_B, mu_S, mu_Q, _ = ys
        chem_potential = array([mu_B, mu_S, mu_Q])
        temp_value = (temperature / self.temperature_0) ** 2
        temp_value += npsum((chem_potential / self.chem_potential_0) ** 2)
        return_value = (
            temperature *
            self.eom.dT_drho(ys, rho) /
            self.temperature_0 ** 2 +
            npsum(
                chem_potential *
                self.eom.dmu_drho(ys, rho) /
                self.chem_potential_0 ** 2
            )
        )
        return 12 * ALPHA * self.temperature_0 ** 4 * temp_value * return_value

    # Utility functions for converting to Milne Coordinates
    def milne_energy(
            self,
            tau: float,
            x: Union[float, ndarray],
            y: Union[float, ndarray],
            q: float,
            ads_T: interp1d,
            ads_mu: List[interp1d],
            tol: float = 1e-20,
    ) -> Union[float, ndarray]:
        r'''
        The de Sitter energy density converted in Milne spacetime

        parameters:
        -----------
        tau float
            Milne time
        x float
            Milne x coordinate
        y float
            Milne y coordinate
        q float
            Dimensional parameter introduce by Gubser to transform for Milne
            coordinates to de Sitter coordinates
        ads_T interp1d
            An interpolating function theat represents the numerical solution
            to the equations of motion temperature
        ads_mu interp1d
            An interpolating function that represents the numerical solution
            to the equations of motion for the chemical potential
        tol float (default 1e-20)
            A lower bound on the evaluations of the interpolating functions
            which are necessary to have physically correct thermodynamic
            variables

        return the Milne energy density evaluated at a point in spacetime
        '''
        r = sqrt(x ** 2 + y ** 2)
        rh = conv.rho(tau, r, q)
        temp = ads_T(rh)
        mu = array([f(rh) for f in ads_mu])

        if isinstance(temp, ndarray):
            return conv.HBARC * self.eos.energy(
                temperature=temp,
                chem_potential=mu,
            ) / tau ** 4
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol
        e = conv.HBARC * self.eos.energy(
            temperature=temp,
            chem_potential=mu,
        ) / tau ** 4
        return tol if e < tol else e

    def milne_number(
            self,
            tau: float,
            x: Union[float, ndarray],
            y: Union[float, ndarray],
            q: float,
            ads_T: interp1d,
            ads_mu: Union[List[interp1d], interp1d],
            tol: float = 1e-20,
    ) -> Union[float, ndarray]:
        r'''
        The de Sitter number density converted in Milne spacetime

        parameters:
        -----------
        tau float
            Milne time
        x float
            Milne x coordinate
        y float
            Milne y coordinate
        q float
            Dimensional parameter introduce by Gubser to transform for Milne
            coordinates to de Sitter coordinates
        ads_T interp1d
            An interpolating function theat represents the numerical solution
            to the equations of motion temperature
        ads_mu interp1d
            An interpolating function that represents the numerical solution
            to the equations of motion for the chemical potential
        tol float (default 1e-20)
            A lower bound on the evaluations of the interpolating functions
            which are necessary to have physically correct thermodynamic
            variables

        return the Milne number density evaluated at a point in spacetime
        '''
        r = sqrt(x ** 2 + y ** 2)
        rh = conv.rho(tau, r, q)
        temp = ads_T(rh)
        if isinstance(ads_mu, List):
            mu = array([f(rh) for f in ads_mu])
        else:
            mu = ads_mu(rh)

        if isinstance(temp, ndarray):
            return self.eos.number(
                temperature=temp,
                chem_potential=mu,
            ) / tau ** 3
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol
        n = self.eos.number(
            temperature=temp,
            chem_potential=mu,
        ) / tau ** 3
        n[where(n < tol)] = tol
        return n

    def milne_entropy(
            self,
            tau: Union[float, ndarray],
            x: Union[float, ndarray],
            y: Union[float, ndarray],
            q: float,
            ads_T: interp1d,
            ads_mu: List[interp1d],
            tol: float = 1e-20,
    ) -> Union[float, ndarray]:
        r'''
        The de Sitter entropy density converted in Milne spacetime

        parameters:
        -----------
        tau float
            Milne time
        x float
            Milne x coordinate
        y float
            Milne y coordinate
        q float
            Dimensional parameter introduce by Gubser to transform for Milne
            coordinates to de Sitter coordinates
        ads_T interp1d
            An interpolating function theat represents the numerical solution
            to the equations of motion temperature
        ads_mu interp1d
            An interpolating function that represents the numerical solution
            to the equations of motion for the chemical potential
        tol float (default 1e-20)
            A lower bound on the evaluations of the interpolating functions
            which are necessary to have physically correct thermodynamic
            variables

        return the Milne entropy density evaluated at a point in spacetime
        '''
        r = sqrt(x ** 2 + y ** 2)
        rh = conv.rho(tau, r, q)
        temp = ads_T(rh)
        mu = array([f(rh) for f in ads_mu])

        if isinstance(temp, ndarray):
            return self.eos.entropy(
                temperature=temp,
                chem_potential=mu,
            ) / tau ** 3
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol
        s = self.eos.entropy(
            temperature=temp,
            chem_potential=mu,
        ) / tau ** 3
        return tol if s < tol else s

    def milne_pi(
            self,
            tau: float,
            x: Union[float, ndarray],
            y: Union[float, ndarray],
            q: float,
            ads_T: interp1d,
            ads_mu: List[interp1d],
            ads_pi_bar_hat: interp1d,
            tol: float = 1e-20,
            nonzero_xy: bool = False,
    ) -> float:
        r'''
        The de Sitter shear pressure converted to the full shear stress tensor
        in Milne spacetime

        parameters:
        -----------
        tau float
            Milne time
        x float
            Milne x coordinate
        y float
            Milne y coordinate
        q float
            Dimensional parameter introduce by Gubser to transform for Milne
            coordinates to de Sitter coordinates
        ads_T interp1d
            An interpolating function theat represents the numerical solution
            to the equations of motion temperature
        ads_mu interp1d
            An interpolating function that represents the numerical solution
            to the equations of motion for the chemical potential
        tol float (default 1e-20)
            A lower bound on the evaluations of the interpolating functions
            which are necessary to have physically correct thermodynamic
            variables

        return the Milne shear stress tensor evaluated at a point in spacetime
        '''
        r = sqrt(x ** 2 + y ** 2)
        temp = ads_T(conv.rho(tau, r, q))
        rh = conv.rho(tau, r, q)
        mu = array([f(rh) for f in ads_mu])

        if isinstance(temp, ndarray):
            pass
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol

        e = self.eos.energy(temp, mu)
        p = self.eos.pressure(temp, mu)

        pi_hat = conv.HBARC * (e + p) * ads_pi_bar_hat(conv.rho(tau, r, q))
        pi_nn = pi_hat / tau ** 6
        pi_xx = -0.5 * (1 + conv.u_x(tau, x, y, q) ** 2) * pi_hat / tau ** 4
        pi_yy = -0.5 * (1 + conv.u_y(tau, x, y, q) ** 2) * pi_hat / tau ** 4
        if nonzero_xy:
            y = x
        pi_xy = -0.5 * conv.u_x(tau, x, y, q) * conv.u_y(tau, x, y, q) * pi_hat / tau ** 4

        if isinstance(temp, ndarray):
            pass
        else:
            pi_nn = tol if fabs(pi_nn) < tol else pi_nn
            pi_xx = tol if fabs(pi_xx) < tol else pi_xx
            pi_yy = tol if fabs(pi_yy) < tol else pi_yy
            pi_xy = tol if fabs(pi_xy) < tol else pi_xy

        return [pi_xx, pi_yy, pi_xy, pi_nn]
