"""
Distances in a flat FRW universe with radiation, matter and a cosmological constant.

    H(z)^2 = H0^2 [ Omega_r (1+z)^4 + Omega_m (1+z)^3 + Omega_L ]

luminosity_distance:  d_L(z) = (1+z) c/H0  int_0^z dz' / E(z')
graviton_distance:    D(z)   = (1+z) c/H0  int_0^z dz' / ((1+z')^2 E(z'))

D(z) is the distance that controls how much a massive graviton lags behind a massless one
(C. M. Will, Phys. Rev. D 57, 2061 (1998)).  Defaults are the Planck 2015 values
(arXiv:1502.01589) with Omega_L fixed by flatness.
"""
from dataclasses import dataclass

import numpy as np

from .constants import C, MPC


@dataclass(frozen=True)
class Cosmology:
    H0: float = 67.74          # km/s/Mpc
    omega_m: float = 0.3089
    omega_r: float = 9.0e-5

    @property
    def omega_l(self):
        return 1.0 - self.omega_m - self.omega_r

    @property
    def hubble_distance_mpc(self):
        """c / H0 in Mpc."""
        return C / 1e3 / self.H0

    def E(self, z):
        z = np.asarray(z, dtype=float)
        return np.sqrt(self.omega_r * (1 + z)**4 + self.omega_m * (1 + z)**3 + self.omega_l)


PLANCK15 = Cosmology()


def _integrate(f, z, n=4001):
    zz = np.linspace(0.0, z, n)
    return np.trapz(f(zz), zz)


def luminosity_distance(z, cosmo=PLANCK15):
    """Luminosity distance d_L(z) in Mpc."""
    return (1 + z) * cosmo.hubble_distance_mpc * _integrate(lambda x: 1 / cosmo.E(x), z)


def graviton_distance(z, cosmo=PLANCK15):
    """The distance D(z) entering the massive-graviton phase, in Mpc."""
    return (1 + z) * cosmo.hubble_distance_mpc * _integrate(
        lambda x: 1 / ((1 + x)**2 * cosmo.E(x)), z)


def mpc_to_m(d_mpc):
    return d_mpc * MPC
