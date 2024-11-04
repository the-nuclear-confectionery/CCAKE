from numpy import ndarray
from numpy import sinh
from numpy import arcsinh
from numpy import arctanh
from numpy import sqrt
from numpy import fabs

from scipy.interpolate import interp1d

from typing import Union

HBARC = 0.19733


def rho(
        tau: float,
        r: Union[float, ndarray],
        q: float,
) -> float:
    '''
    Calculate de Sitter time from cylindral Milne coordinates

    parameters:
    -----------
    tau float
        Milne time
    r float
        Milne transverse radius
    q float
        Dimensional parameter introduce by Gubser to transform for Milne
        coordinates to de Sitter coordinates

    returns the de Sitter time
    '''
    return arcsinh(-(1 - (q * tau) ** 2 + (q * r) ** 2) / (2 * q * tau))


def kappa(
        tau: float,
        r: Union[float, ndarray],
        q: float,
) -> float:
    '''
    A conveience function that allows us to represent Gubser solution for
    the fluid's velocity field in Milne coordinates in concise manner.

    parameters:
    -----------
    tau float
        Milne time
    r float
        Milne transverse radius
    q float
        Dimensional parameter introduce by Gubser to transform for Milne
        coordinates to de Sitter coordinates

    returns the conveince function evaluated at the point in Milne
    spacetime
    '''
    return arctanh((2 * q ** 2 * r * tau) /
                   (1 + (q * tau) ** 2 + (q * r) ** 2))


def u_r(
        tau: float,
        r: Union[float, ndarray],
        q: float,
) -> float:
    '''
    The radial component of the Gubser solution for the fluid's velocity
    field in Milne coordinates

    parameters:
    -----------
    tau float
        Milne time
    r float
        Milne transverse radius
    q float
        Dimensional parameter introduce by Gubser to transform for Milne
        coordinates to de Sitter coordinates

    returns the radial velocity field evaluated at the point in Milne
    spacetime
    '''
    return sinh(kappa(tau, r, q))


def u_x(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
) -> float:
    '''
    The x-component of the Gubser solution for the fluid's velocity
    field in Milne coordinates

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

    returns the x-component of the velocity field evaluated at the point
    in Milne spacetime
    '''
    r = sqrt(x ** 2 + y ** 2)
    return (x / r) * sinh(kappa(tau, r, q))


def u_y(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
) -> float:
    '''
    The y-component of the Gubser solution for the fluid's velocity
    field in Milne coordinates

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

    returns the y-component of the velocity field evaluated at the point
    in Milne spacetime
    '''
    r = sqrt(x ** 2 + y ** 2)
    return (y / r) * sinh(kappa(tau, r, q))


def milne_T(
        tau: float,
        r: Union[float, ndarray],
        q: float,
        ads_T: interp1d,
) -> float:
    '''
    The de Sitter temperature converted in Milne spacetime

    parameters:
    -----------
    tau float
        Milne time
    r float
        Milne transverse radius
    q float
        Dimensional parameter introduce by Gubser to transform for Milne
        coordinates to de Sitter coordinates
    ads_T interp1d
        An interpolating function that represents the numerical solution
        to the equations of motion for the temperature

    return the Milne temperature evaluated at a point in spacetime
    '''
    return HBARC * ads_T(rho(tau, r, q)) / tau


def milne_mu(
        tau: float,
        r: Union[float, ndarray],
        q: float,
        ads_mu: interp1d
) -> float:
    '''
    The de Sitter chemical potentials converted in Milne spacetime

    parameters:
    -----------
    tau float
        Milne time
    r float
        Milne transverse radius
    q float
        Dimensional parameter introduce by Gubser to transform for Milne
        coordinates to de Sitter coordinates
    ads_mu interp1d
        An interpolating function that represents the numerical solution
        to the equations of motion for the chemical potential

    return the Milne chemical potential evaluated at a point in spacetime
    '''
    return HBARC * ads_mu(rho(tau, r, q)) / tau