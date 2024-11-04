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
NC = 3                                          # number of colors
NF = 2.5                                        # number of flavors
ALPHA = (2 * (NC **2 - 1) + 7 * NC * NF / 2)    # qgp degeneracy
BETA = 4 * NC * NF                              # fermi gas degeneracy
CTAUR = 5                                       # relaxation time constant
ETA_S = 0.2                                     # shear viscosity to entory ratio


# Equation of State class
class EoS:
    '''
    Corresponds to EoS 1 in the paper and describes a masslass quark gluon
    plasma
    '''
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
        return_value = BETA * chem_potential ** 2 * temperature ** 2 / 216
        return_value += BETA * chem_potential ** 4 / 324 / pi ** 2
        return_value += ALPHA * pi ** 2 * temperature ** 4 / 90
        return return_value

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
        return 3 * self.pressure(temperature, chem_potential)

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
        return_value = 2 * pi ** 2 * ALPHA * temperature ** 3 / 45
        return_value += BETA * temperature * chem_potential ** 2 / 108
        return return_value

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
        return_value = BETA * temperature ** 2 * chem_potential / 108
        return_value += BETA * chem_potential ** 3 / 81 / pi ** 2
        return ALPHA * return_value

# The equations of motion class
class EoM:
    '''
    This class implementst the equations of motion for temperature, chemcical
    potentails and shear pressure for the viscous Gubser with BSQ currents
    '''
    def __init__(
            self,
            eos: EoS
    ):
        self.eos = eos

    def tau_R(
            self,
            temperature: float,
            chem_potential: float,
    ) -> float:
        '''
        The relaxation time calculate by calling the equation of state with the
        natural equation of state variables

        parameters:
        ----------
        temperature float
        chem_potential float
        '''
        e = self.eos.energy(temperature=temperature, chem_potential=chem_potential)
        p = self.eos.pressure(temperature=temperature, chem_potential=chem_potential)
        s = self.eos.entropy(temperature=temperature, chem_potential=chem_potential)
        return CTAUR * ETA_S * s / (e + p)

    def f(
            self,
            temperature: float,
            chem_potential: float,
    ) -> float:
        numerator = 5 * BETA * 3 * pi ** 2 * chem_potential ** 2 * temperature ** 2
        numerator += 5 * BETA * 2 * chem_potential ** 4
        numerator += 36 * pi ** 2 * ALPHA * temperature ** 4

        denominator = 72 * pi ** 4 * ALPHA * temperature ** 4
        denominator += 288 * pi ** 2 * ALPHA * temperature ** 2 * chem_potential ** 2
        denominator += -15 * pi ** 2 * BETA * chem_potential ** 2 * temperature ** 2
        denominator += 20 * BETA * chem_potential ** 4
        return numerator / denominator

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
        temperature, chem_potential, pi_hat = ys
        return_value = 1 + (4 * chem_potential ** 2) / (pi ** 2 * temperature ** 2)
        return_value *= -self.f(temperature, chem_potential) * pi_hat
        return_value += 1
        return -(2 / 3) * temperature * return_value * tanh(rho)

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
        temperature, chem_potential, pi_hat = ys
        return_value = 1 + 2 * self.f(temperature, chem_potential) * pi_hat
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
        temperature, chem_potenial, pi_hat = ys
        tau_r = self.tau_R(temperature=temperature, chem_potential=chem_potenial)
        return_value = (4 / 3 / CTAUR) * tanh(rho)
        return_value -= pi_hat / tau_r
        return_value -= (4 / 3) * pi_hat ** 2 * tanh(rho)
        return return_value

    def for_scipy(
            self,
            ys: ndarray,
            rho: Union[float, ndarray],
            ideal: bool = False,
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
            If true, sets the viscous contributions to 0

        return dpi/drho evaluated at time rho
        '''
        dTdrho = self.dT_drho(ys, rho)
        dmudrho = self.dmu_drho(ys, rho)
        dpidrho = 0 if ideal else self.dpi_drho(ys, rho)
        return array([dTdrho, *dmudrho.reshape(-1,), dpidrho])


class MasslessQGP:
    '''
    Class that contains the implementation details for a massless QGP equation
    of state
    '''

    def __init__(self):
        self.eos = EoS()
        self.eom = EoM(self.eos)

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
        temperature, chem_potential, _ = ys
        return_value_1 = 24 * ALPHA * pi ** 2 * temperature ** 2
        return_value_1 += 5 * BETA * chem_potential ** 2
        return_value_1 *= 3 * temperature * self.eom.dT_drho(ys, rho)
        return_value_2 = 3 * pi **2 * temperature ** 2
        return_value_2 += 4 * chem_potential ** 2
        return_value_2 *= 5 * BETA * chem_potential * self.eom.dmu_drho(ys, rho) / pi ** 2
        return (return_value_1 + return_value_2) / 540

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
        mu = ads_mu(rh)

        if isinstance(temp, ndarray):
            return conv.HBARC * self.eos.energy(temperature=temp, chem_potential=mu) / tau ** 4
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol
        e = conv.HBARC * self.eos.energy(temperature=temp, chem_potential=mu) / tau ** 4
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
        mu = ads_mu(rh)

        if isinstance(temp, ndarray):
            return self.eos.number(temperature=temp, chem_potential=mu) / tau ** 3
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol
        n = self.eos.number(temperature=temp, chem_potential=mu) / tau ** 3
        return tol if n < tol else n

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
        mu = ads_mu(rh)

        if isinstance(temp, ndarray):
            return self.eos.entropy(temperature=temp, chem_potential=mu) / tau ** 3
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol
        s = self.eos.entropy(temperature=temp, chem_potential=mu) / tau ** 3
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
        mu = ads_mu(conv.rho(tau, r, q))

        if isinstance(temp, ndarray):
            pass
        else:
            if temp <= tol:
                temp = tol
            if mu <= tol:
                temp = tol

        e = self.eos.energy(temperature=temp, chem_potential=mu)
        p = self.eos.pressure(temperature=temp, chem_potential=mu)

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