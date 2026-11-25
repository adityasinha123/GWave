"""
A simplified LISA response, following Cutler, Phys. Rev. D 57, 7089 (1998) and
Barack & Cutler, Phys. Rev. D 69, 082005 (2004).

LISA's three arms act as two L-shaped detectors I and II with 60-degree opening angle,
h_A = (sqrt(3)/2) [F^A_+ h_+ + F^A_x h_x].  Detector II is detector I rotated by 45 degrees,
F^II(theta, phi, psi) = F^I(theta, phi - pi/4, psi).

Simplifications: the constellation's orientation is held fixed (fine for signals much
shorter than a year), but the Doppler phase from its orbit around the Sun is kept, since
it matters already for signals lasting weeks:

    Phi_D(t) = 2 pi f(t) (R/c) sin(theta_s) cos(2 pi t / 1 yr - phi_s)

with (theta_s, phi_s) the source position in ecliptic coordinates and R = 1 AU.
"""
import numpy as np

from .constants import AU, C, YEAR
from .inspiral import antenna_pattern, gw_frequency, polarizations


def lisa_antenna_patterns(theta, phi, psi):
    """(F^I_+, F^I_x, F^II_+, F^II_x) in the detector frame."""
    f1p, f1x = antenna_pattern(theta, phi, psi)
    f2p, f2x = antenna_pattern(theta, phi - np.pi / 4, psi)
    return f1p, f1x, f2p, f2x


def doppler_phase(f_gw, t, theta_s, phi_s, t0=0.0):
    """Phase shift from LISA's heliocentric orbit; t0 sets the orbital phase at t = 0."""
    return 2 * np.pi * f_gw * (AU / C) * np.sin(theta_s) * np.cos(2 * np.pi * (t + t0) / YEAR - phi_s)


def lisa_strain(binary, t, distance_mpc, iota, theta=0.0, phi=0.0, psi=0.0,
                theta_s=None, phi_s=None, doppler=True, pn_order=2):
    """
    Strains (h_I, h_II) seen by LISA. (theta, phi) is the source direction in the detector
    frame, (theta_s, phi_s) in ecliptic coordinates (default: the same angles).
    """
    theta_s = theta if theta_s is None else theta_s
    phi_s = phi if phi_s is None else phi_s
    extra = doppler_phase(gw_frequency(binary, t, pn_order=pn_order), t, theta_s, phi_s) if doppler else 0.0
    hp, hx = polarizations(binary, t, distance_mpc, iota, pn_order=pn_order, extra_phase=extra)
    f1p, f1x, f2p, f2x = lisa_antenna_patterns(theta, phi, psi)
    k = np.sqrt(3) / 2
    return k * (f1p * hp + f1x * hx), k * (f2p * hp + f2x * hx)
