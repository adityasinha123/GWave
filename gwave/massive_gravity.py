"""
Propagation of gravitational waves with a massive graviton, E^2 = p^2 c^2 + m_g^2 c^4.

Over a cosmological distance each frequency picks up a phase (Will 1998)

    Psi_g(f) = - pi c D(z) / ((1+z) lambda_g^2 f),          lambda_g = h / (m_g c),

with D(z) from gwave.cosmology.graviton_distance.  Lower frequencies travel more slowly, so
the frequency f arrives late by the group delay

    dt(f) = (1/2pi) dPsi_g/df = c D / (2 (1+z) lambda_g^2 f^2),

and its phase at arrival is shifted by -2 Psi_g(f) (the phase and group velocities differ).
For a chirping inspiral this compresses the waveform: its early, low-frequency part is
delayed more than its late, high-frequency part. Compression makes each frequency sweep past
faster, so the amplitude grows by 1/sqrt(dt_arrival/dt) to keep |h(f)| unchanged.

`propagate` applies this directly in the time domain, sample by sample. That avoids taking
an FFT of a finite stretch of signal, whose circular wrap-around and edge effects otherwise
have to be tamed with windowing.
"""
import numpy as np

from .constants import C, HC_EV_M, MPC
from .cosmology import PLANCK15, graviton_distance


def compton_wavelength(m_g_ev):
    """Graviton Compton wavelength lambda_g = h / (m_g c) in metres, for m_g in eV/c^2."""
    return HC_EV_M / m_g_ev


def graviton_mass_ev(lambda_g):
    """Graviton mass in eV/c^2 for a Compton wavelength in metres."""
    return HC_EV_M / lambda_g


def _scale(z, lambda_g, cosmo):
    """c D / ((1+z) lambda_g^2), in Hz (so that Psi_g = -pi * scale / f)."""
    D = graviton_distance(z, cosmo) * MPC
    return C * D / ((1 + z) * lambda_g**2)


def phase_shift(f, z, lambda_g, cosmo=PLANCK15):
    """Psi_g(f) in radians, for observed frequency f [Hz]."""
    return -np.pi * _scale(z, lambda_g, cosmo) / np.asarray(f, dtype=float)


def group_delay(f, z, lambda_g, cosmo=PLANCK15):
    """Arrival delay of frequency f relative to a massless graviton [s]."""
    return 0.5 * _scale(z, lambda_g, cosmo) / np.asarray(f, dtype=float)**2


def propagate(t, f, phase, z, lambda_g, f_ref=None, cosmo=PLANCK15):
    """
    Map a GR signal, given by its sample times t, frequency f(t) and phase Phi(t), to the
    signal after massive-graviton propagation.
    Returns (t_arrival, phase_arrival, amplitude_factor).

    The constant delay and phase of the reference frequency f_ref (default: the last
    sample) are removed: they are absorbed into the coalescence time and phase.
    """
    f = np.asarray(f, dtype=float)
    f_ref = f[-1] if f_ref is None else f_ref
    t_arr = t + group_delay(f, z, lambda_g, cosmo) - group_delay(f_ref, z, lambda_g, cosmo)
    ph_arr = phase - 2 * (phase_shift(f, z, lambda_g, cosmo) - phase_shift(f_ref, z, lambda_g, cosmo))
    if np.any(np.diff(t_arr) <= 0):
        raise ValueError("arrival times not monotonic: graviton mass too large for this signal")
    amp = 1.0 / np.sqrt(np.gradient(t_arr, t))
    return t_arr, ph_arr, amp


def massive_inspiral(binary, z, lambda_g, iota=0.0, f_low=20.0, sample_rate=8192.0,
                     theta_sky=0.0, phi_sky=0.0, psi=0.0, pn_order=2, cosmo=PLANCK15):
    """
    Strain of an L-shaped detector for a source at redshift z, in GR and with a massive
    graviton, on a common time grid (GR coalescence at t = 0). Returns (t, h_gr, h_massive).
    """
    from .cosmology import luminosity_distance
    from .inspiral import antenna_pattern, gw_frequency, gw_phase, inspiral_times

    b = binary.redshifted(z)
    dL = luminosity_distance(z, cosmo)
    t = inspiral_times(b, f_low, sample_rate, pn_order)
    f = gw_frequency(b, t, pn_order=pn_order)
    phase = gw_phase(b, t, pn_order=pn_order)
    t_arr, ph_arr, amp = propagate(t, f, phase, z, lambda_g, cosmo=cosmo)

    fp, fx = antenna_pattern(theta_sky, phi_sky, psi)
    r_s = dL * MPC / C

    def h_of(f_, ph_, a_=1.0):
        h0 = a_ * 4.0 * b.Mc_s**(5 / 3) * (np.pi * f_)**(2 / 3) / r_s
        return h0 * (fp * 0.5 * (1 + np.cos(iota)**2) * np.cos(ph_) + fx * np.cos(iota) * np.sin(ph_))

    h_gr = h_of(f, phase)
    inside = (t >= t_arr[0]) & (t <= t_arr[-1])
    h_m = np.zeros_like(t)
    h_m[inside] = h_of(np.interp(t[inside], t_arr, f), np.interp(t[inside], t_arr, ph_arr),
                       np.interp(t[inside], t_arr, amp))
    return t, h_gr, h_m
