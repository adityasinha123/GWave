"""
Post-Newtonian inspiral of a compact binary and its signal in an L-shaped detector.

Time-domain "restricted" waveforms: the phase is evolved to 0PN (Newtonian) or 2PN order
in the dimensionless time variable Theta of L. Blanchet, Living Rev. Relativ. 17, 2 (2014);
the amplitude is kept at leading (quadrupole) order.

All masses are in solar masses, distances in Mpc, times in seconds, frequencies in Hz.
Internally masses are converted to seconds via T_sun = G M_sun / c^3.

For a source at redshift z, every mass scale in the signal is redshifted, M -> (1+z) M,
and the distance is the luminosity distance: use Binary.redshifted(z) together with
gwave.cosmology.luminosity_distance(z).
"""
from dataclasses import dataclass, replace

import numpy as np

from .constants import C, MPC, T_SUN


# --- the binary ------------------------------------------------------------------------
@dataclass(frozen=True)
class Binary:
    """Two point masses m1, m2 in solar masses."""
    m1: float
    m2: float

    def redshifted(self, z):
        """The same binary as seen from redshift z: detector-frame masses (1+z) m."""
        return replace(self, m1=self.m1 * (1 + z), m2=self.m2 * (1 + z))

    @property
    def total_mass(self):
        return self.m1 + self.m2

    @property
    def eta(self):
        """Symmetric mass ratio  nu = m1 m2 / M^2  (1/4 for equal masses)."""
        return self.m1 * self.m2 / self.total_mass**2

    @property
    def chirp_mass(self):
        """Chirp mass  M_c = mu^(3/5) M^(2/5) = eta^(3/5) M  [M_sun]."""
        return self.eta**0.6 * self.total_mass

    @property
    def M_s(self):
        """Total mass in seconds, G M / c^3."""
        return self.total_mass * T_SUN

    @property
    def Mc_s(self):
        """Chirp mass in seconds, G M_c / c^3."""
        return self.chirp_mass * T_SUN

    @property
    def f_isco(self):
        """GW frequency at the Schwarzschild ISCO, where the inspiral approximation ends."""
        return 1.0 / (6**1.5 * np.pi * self.M_s)

    def time_to_coalescence(self, f_gw):
        """Newtonian time from GW frequency f_gw to coalescence [s]."""
        return 5.0 / 256.0 * self.Mc_s**(-5 / 3) * (np.pi * f_gw)**(-8 / 3)


# --- phase evolution ---------------------------------------------------------------------
def theta(binary, t, tc=0.0):
    """Dimensionless time to coalescence  Theta = nu (tc - t) / (5 G M / c^3)."""
    return binary.eta * (tc - np.asarray(t, dtype=float)) / (5.0 * binary.M_s)


def _pn_coefficients(eta):
    a1 = 743 / 2688 + 11 / 32 * eta                                     # 1PN
    a2 = 1855099 / 14450688 + 56975 / 258048 * eta + 371 / 2048 * eta**2  # 2PN
    b1 = 3715 / 8064 + 55 / 96 * eta
    b2 = 9275495 / 14450688 + 284875 / 258048 * eta + 1855 / 2048 * eta**2
    return a1, a2, b1, b2


def gw_frequency(binary, t, tc=0.0, pn_order=2):
    """GW frequency f(t) = 2 x orbital frequency [Hz], to 0PN or 2PN order."""
    th = theta(binary, t, tc)
    series = th**(-3 / 8)
    if pn_order == 2:
        a1, a2, _, _ = _pn_coefficients(binary.eta)
        series = series + a1 * th**(-5 / 8) - 3 * np.pi / 10 * th**(-3 / 4) + a2 * th**(-7 / 8)
    elif pn_order != 0:
        raise ValueError("pn_order must be 0 or 2")
    return series / (8 * np.pi * binary.M_s)


def gw_phase(binary, t, tc=0.0, phi_c=0.0, pn_order=2):
    """GW phase Phi(t) = 2 x orbital phase, to 0PN or 2PN order.  Phi(tc) = phi_c."""
    th = theta(binary, t, tc)
    series = th**(5 / 8)
    if pn_order == 2:
        _, _, b1, b2 = _pn_coefficients(binary.eta)
        series = series + b1 * th**(3 / 8) - 3 * np.pi / 4 * th**(1 / 4) + b2 * th**(1 / 8)
    elif pn_order != 0:
        raise ValueError("pn_order must be 0 or 2")
    return phi_c - 2.0 / binary.eta * series


# --- waveform ------------------------------------------------------------------------------
def polarizations(binary, t, distance_mpc, iota, tc=0.0, phi_c=0.0, pn_order=2,
                  extra_phase=0.0):
    """
    Restricted-PN polarizations h_+(t), h_x(t) for inclination iota (radians).
    extra_phase is added to the GW phase (e.g. the Doppler phase of a moving detector).
    """
    f = gw_frequency(binary, t, tc, pn_order)
    phase = gw_phase(binary, t, tc, phi_c, pn_order) + extra_phase
    r_s = distance_mpc * MPC / C                                  # distance in seconds
    h0 = 4.0 * binary.Mc_s**(5 / 3) * (np.pi * f)**(2 / 3) / r_s
    h_plus = h0 * 0.5 * (1 + np.cos(iota)**2) * np.cos(phase)
    h_cross = h0 * np.cos(iota) * np.sin(phase)
    return h_plus, h_cross


def antenna_pattern(theta_sky, phi_sky, psi):
    """
    Response F_+, F_x of an L-shaped interferometer with arms along x and y, for a source
    at polar angle theta_sky and azimuth phi_sky (detector frame), polarization angle psi.
    """
    a = 0.5 * (1 + np.cos(theta_sky)**2) * np.cos(2 * phi_sky)
    b = np.cos(theta_sky) * np.sin(2 * phi_sky)
    f_plus = a * np.cos(2 * psi) - b * np.sin(2 * psi)
    f_cross = a * np.sin(2 * psi) + b * np.cos(2 * psi)
    return f_plus, f_cross


def strain(binary, t, distance_mpc, iota, theta_sky=0.0, phi_sky=0.0, psi=0.0,
           tc=0.0, phi_c=0.0, pn_order=2):
    """Detector strain h(t) = F_+ h_+ + F_x h_x."""
    hp, hx = polarizations(binary, t, distance_mpc, iota, tc, phi_c, pn_order)
    fp, fx = antenna_pattern(theta_sky, phi_sky, psi)
    return fp * hp + fx * hx


def gw_cycles(binary, f_low, f_high=None, pn_order=2, n=400_000):
    """
    Number of GW cycles N(f) accumulated from f_low, evaluated on a log grid up to f_high
    (default: the ISCO).  Returns (f_grid, N).
    """
    f_high = binary.f_isco if f_high is None else f_high
    t = -np.geomspace(1.5 * binary.time_to_coalescence(f_low), 1e-6, n)   # dense near merger
    f = gw_frequency(binary, t, pn_order=pn_order)
    phase = gw_phase(binary, t, pn_order=pn_order)
    increasing = np.concatenate([[True], np.diff(f) > 0])
    stop = np.argmin(increasing) if not increasing.all() else len(f)
    f, phase = f[:stop], phase[:stop]
    grid = np.geomspace(f_low, f_high, 300)
    return grid, (np.interp(grid, f, phase) - np.interp(f_low, f, phase)) / (2 * np.pi)


def inspiral_times(binary, f_low, sample_rate, pn_order=2):
    """
    Sample times from GW frequency f_low up to the ISCO (or up to where the truncated PN
    series stops increasing), with coalescence at t = 0.
    """
    t = np.arange(-binary.time_to_coalescence(f_low), 0.0, 1.0 / sample_rate)
    f = gw_frequency(binary, t, pn_order=pn_order)
    keep = (f < binary.f_isco) & np.concatenate([[True], np.diff(f) > 0])
    stop = np.argmin(keep) if not keep.all() else len(t)
    return t[:stop]


def inspiral(binary, distance_mpc=None, iota=0.0, f_low=20.0, sample_rate=8192.0,
             theta_sky=0.0, phi_sky=0.0, psi=0.0, pn_order=2, z=None):
    """
    Convenience wrapper: the strain of an L-shaped detector from f_low up to the ISCO.
    Give either distance_mpc (a nearby source) or z (masses are taken as source-frame,
    redshifted, and the luminosity distance is computed). Returns (t, h), coalescence at t = 0.
    """
    if z is not None:
        from .cosmology import luminosity_distance
        binary = binary.redshifted(z)
        distance_mpc = luminosity_distance(z) if distance_mpc is None else distance_mpc
    if distance_mpc is None:
        raise ValueError("give distance_mpc or z")
    t = inspiral_times(binary, f_low, sample_rate, pn_order)
    h = strain(binary, t, distance_mpc, iota, theta_sky, phi_sky, psi, pn_order=pn_order)
    return t, h
