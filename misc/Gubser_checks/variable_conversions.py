from numpy import ndarray
from numpy import sinh
from numpy import arcsinh
from numpy import arctanh
from numpy import sqrt
from numpy import fabs
from numpy import where
from numpy import array

from scipy.interpolate import interp1d

from typing import Union
from typing import List

from equations_of_motion import energy
from equations_of_motion import number
from equations_of_motion import pressure
from equations_of_motion import entropy

HBARC = 0.19733


def rho(
        tau: float,
        r: Union[float, ndarray],
        q: float,
) -> float:
    return arcsinh(-(1 - (q * tau) ** 2 + (q * r) ** 2) / (2 * q * tau))


def kappa(
        tau: float,
        r: Union[float, ndarray],
        q: float,
) -> float:
    return arctanh((2 * q ** 2 * r * tau) /
                   (1 + (q * tau) ** 2 + (q * r) ** 2))


def u_r(
        tau: float,
        r: Union[float, ndarray],
        q: float,
) -> float:
    return sinh(kappa(tau, r, q))


def u_x(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
) -> float:
    r = sqrt(x ** 2 + y ** 2)
    return (x / r) * sinh(kappa(tau, r, q))


def u_y(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
) -> float:
    r = sqrt(x ** 2 + y ** 2)
    return (y / r) * sinh(kappa(tau, r, q))


def milne_T(
        tau: float,
        r: Union[float, ndarray],
        q: float,
        ads_T: interp1d,
) -> float:
    return HBARC * ads_T(rho(tau, r, q)) / tau


def milne_mu(
        tau: float,
        r: Union[float, ndarray],
        q: float,
        ads_mu: interp1d
) -> float:
    return HBARC * ads_mu(rho(tau, r, q)) / tau


def milne_energy(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
        ads_T: interp1d,
        ads_mu: List[interp1d],
        temperature_0: float,
        chem_potential_0: ndarray,
        tol: float = 1e-20,
) -> Union[float, ndarray]:
    r = sqrt(x ** 2 + y ** 2)
    rh = rho(tau, r, q)
    temp = ads_T(rh)
    mu = array([f(rh) for f in ads_mu])

    if isinstance(temp, ndarray):
        return HBARC * energy(
            temperature=temp,
            chem_potential=mu,
            temperature_0=temperature_0,
            chem_potential_0=chem_potential_0
        ) / tau ** 4
    else:
        if temp <= tol:
            temp = tol
        if mu <= tol:
            temp = tol
    e = HBARC * energy(
        temperature=temp,
        chem_potential=mu,
        temperature_0=temperature_0,
        chem_potential_0=chem_potential_0
    ) / tau ** 4
    return tol if e < tol else e


def milne_number(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
        ads_T: interp1d,
        ads_mu: Union[List[interp1d], interp1d],
        temperature_0: float,
        chem_potential_0: ndarray,
        tol: float = 1e-20,
) -> Union[float, ndarray]:
    r = sqrt(x ** 2 + y ** 2)
    rh = rho(tau, r, q)
    temp = ads_T(rh)
    if isinstance(ads_mu, List):
        mu = array([f(rh) for f in ads_mu])
    else:
        mu = ads_mu(rh)

    if isinstance(temp, ndarray):
        return number(
            temperature=temp,
            chem_potential=mu,
            temperature_0=temperature_0,
            chem_potential_0=chem_potential_0
        ) / tau ** 3
    else:
        if temp <= tol:
            temp = tol
        if mu <= tol:
            temp = tol
    n = number(
        temperature=temp,
        chem_potential=mu,
        temperature_0=temperature_0,
        chem_potential_0=chem_potential_0
    ) / tau ** 3
    n[where(n < tol)] = tol
    return n


def milne_entropy(
        tau: Union[float, ndarray],
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
        ads_T: interp1d,
        ads_mu: List[interp1d],
        temperature_0: float,
        chem_potential_0: ndarray,
        tol: float = 1e-20,
) -> Union[float, ndarray]:
    r = sqrt(x ** 2 + y ** 2)
    rh = rho(tau, r, q)
    temp = ads_T(rh)
    mu = array([f(rh) for f in ads_mu])

    if isinstance(temp, ndarray):
        return entropy(
            temperature=temp,
            chem_potential=mu,
            temperature_0=temperature_0,
            chem_potential_0=chem_potential_0
        ) / tau ** 3
    else:
        if temp <= tol:
            temp = tol
        if mu <= tol:
            temp = tol
    s = entropy(
        temperature=temp,
        chem_potential=mu,
        temperature_0=temperature_0,
        chem_potential_0=chem_potential_0
    ) / tau ** 3
    return tol if s < tol else s


def milne_pi(
        tau: float,
        x: Union[float, ndarray],
        y: Union[float, ndarray],
        q: float,
        ads_T: interp1d,
        ads_mu: List[interp1d],
        ads_pi_bar_hat: interp1d,
        temperature_0: float,
        chem_potential_0: ndarray,
        tol: float = 1e-20,
        nonzero_xy: bool = False,
) -> float:
    r = sqrt(x ** 2 + y ** 2)
    temp = ads_T(rho(tau, r, q))
    rh = rho(tau, r, q)
    mu = array([f(rh) for f in ads_mu])

    if isinstance(temp, ndarray):
        pass
    else:
        if temp <= tol:
            temp = tol
        if mu <= tol:
            temp = tol

    e = energy(temp, mu, temperature_0, chem_potential_0)
    p = pressure(temp, mu, temperature_0, chem_potential_0)

    pi_hat = HBARC * (e + p) * ads_pi_bar_hat(rho(tau, r, q))
    pi_nn = pi_hat / tau ** 6
    pi_xx = -0.5 * (1 + u_x(tau, x, y, q) ** 2) * pi_hat / tau ** 4
    pi_yy = -0.5 * (1 + u_y(tau, x, y, q) ** 2) * pi_hat / tau ** 4
    if nonzero_xy:
        y = x
    pi_xy = -0.5 * u_x(tau, x, y, q) * u_y(tau, x, y, q) * pi_hat / tau ** 4

    if isinstance(temp, ndarray):
        pass
    else:
        pi_nn = tol if fabs(pi_nn) < tol else pi_nn
        pi_xx = tol if fabs(pi_xx) < tol else pi_xx
        pi_yy = tol if fabs(pi_yy) < tol else pi_yy
        pi_xy = tol if fabs(pi_xy) < tol else pi_xy

    return [pi_xx, pi_yy, pi_xy, pi_nn]
