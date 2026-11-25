"""Tests for cosmology, the LISA response and massive-graviton propagation."""
import numpy as np
import pytest

import gwave
from gwave import Binary, cosmology, lisa
from gwave import massive_gravity as mg
from gwave.constants import AU, C


# --- cosmology -------------------------------------------------------------------------------
def test_low_redshift_hubble_law():
    z = 1e-3
    assert cosmology.luminosity_distance(z) == pytest.approx(z * cosmology.PLANCK15.hubble_distance_mpc, rel=2e-3)


def test_matter_only_closed_forms():
    """For Omega_m = 1 both distances have closed forms (Will 1998)."""
    eds = cosmology.Cosmology(H0=70.0, omega_m=1.0, omega_r=0.0)
    dh = eds.hubble_distance_mpc
    for z in [0.1, 1.0, 3.0]:
        assert cosmology.luminosity_distance(z, eds) == pytest.approx(2 * dh * (1 + z) * (1 - (1 + z)**-0.5), rel=1e-6)
        assert cosmology.graviton_distance(z, eds) == pytest.approx(0.4 * dh * (1 + z) * (1 - (1 + z)**-2.5), rel=1e-6)


def test_redshift_scales_masses_not_mass_ratio():
    b, bz = Binary(36, 29), Binary(36, 29).redshifted(0.09)
    assert bz.chirp_mass == pytest.approx(1.09 * b.chirp_mass)
    assert bz.eta == pytest.approx(b.eta)


def test_inspiral_with_redshift():
    b, z = Binary(36, 29), 0.09
    t1, h1 = gwave.inspiral(b, z=z, f_low=30)
    t2, h2 = gwave.inspiral(b.redshifted(z), distance_mpc=cosmology.luminosity_distance(z), f_low=30)
    np.testing.assert_allclose(h1, h2)


# --- LISA ----------------------------------------------------------------------------------
def test_lisa_detector_two_is_rotated_detector_one():
    """Matches the explicit F^II expressions in Cutler (1998)."""
    th, ph, ps = 0.7, 1.1, 0.3
    _, _, f2p, f2x = lisa.lisa_antenna_patterns(th, ph, ps)
    a = 0.5 * (1 + np.cos(th)**2)
    assert f2p == pytest.approx(a * np.sin(2 * ph) * np.cos(2 * ps) + np.cos(th) * np.cos(2 * ph) * np.sin(2 * ps))
    assert f2x == pytest.approx(a * np.sin(2 * ph) * np.sin(2 * ps) - np.cos(th) * np.cos(2 * ph) * np.cos(2 * ps))


def test_doppler_phase_amplitude():
    t = np.linspace(0, gwave.YEAR, 10001)
    d = lisa.doppler_phase(0.01, t, np.pi / 2, 0.0)
    assert np.max(np.abs(d)) == pytest.approx(2 * np.pi * 0.01 * AU / C, rel=1e-6)


# --- massive graviton ----------------------------------------------------------------------
def test_compton_wavelength_of_ligo_bound():
    # m_g < 1.27e-23 eV  <->  lambda_g > ~9.8e16 m
    assert mg.compton_wavelength(1.27e-23) == pytest.approx(9.76e16, rel=1e-2)


def test_group_delay_is_phase_derivative():
    f = np.linspace(10, 500, 200001)
    psi = mg.phase_shift(f, 0.2, 1e16)
    np.testing.assert_allclose(np.gradient(psi, f)[1:-1] / (2 * np.pi),
                               mg.group_delay(f, 0.2, 1e16)[1:-1], rtol=1e-6)


def test_massless_limit():
    t, h_gr, h_m = mg.massive_inspiral(Binary(30, 30), 0.1, 1e24, f_low=30)
    assert np.max(np.abs(h_m - h_gr)) < 1e-6 * np.max(np.abs(h_gr))


def test_matches_frequency_domain_propagation():
    """
    Independent check: multiply the (zero-padded, tapered) FFT of a GR chirp by the
    propagation phase and transform back; compare with the time-domain mapping.
    """
    b, z, fs, lam = Binary(5, 5), 0.1, 2048.0, 5e15
    bz = b.redshifted(z)
    t = gwave.inspiral_times(bz, 30, fs)
    f = gwave.gw_frequency(bz, t)
    t, f = t[f < 150], f[f < 150]
    ph = gwave.gw_phase(bz, t)
    n = len(t)
    taper = np.ones(n)
    k = n // 10
    ramp = 0.5 - 0.5 * np.cos(np.pi * np.arange(k) / k)
    taper[:k], taper[-k:] = ramp, ramp[::-1]
    a = (np.pi * f)**(2 / 3) * taper

    # frequency domain (numpy's FFT sign convention: multiply by exp(-i Psi_g))
    pad, off = 8 * n, 3 * n
    x = np.zeros(pad)
    x[off:off + n] = a * np.cos(ph)
    fr = np.fft.rfftfreq(pad, 1 / fs)
    fref = f[-1]
    psi = np.zeros_like(fr)
    psi[1:] = -(mg.phase_shift(fr[1:], z, lam) - mg.phase_shift(fref, z, lam)
                - 2 * np.pi * mg.group_delay(fref, z, lam) * (fr[1:] - fref))
    y = np.fft.irfft(np.fft.rfft(x) * np.exp(1j * psi), pad)[off:off + n]

    # time domain
    t_arr, ph_arr, amp = mg.propagate(t, f, ph, z, lam)
    h_m = np.interp(t, t_arr, a * amp * np.cos(ph_arr), left=0, right=0)

    gr_vs_massive = np.max(np.abs(y - a * np.cos(ph))) / np.max(np.abs(y))
    assert gr_vs_massive > 1.0                                   # the effect is large here
    assert np.max(np.abs(y - h_m)) / np.max(np.abs(y)) < 0.01    # yet the two methods agree
