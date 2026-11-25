import numpy as np
import pytest

import gwave
from gwave import Binary


def test_mass_parameters():
    b = Binary(1.4, 1.4)
    assert b.total_mass == pytest.approx(2.8)
    assert b.eta == pytest.approx(0.25)
    assert b.chirp_mass == pytest.approx(1.2188, rel=1e-4)       # standard BNS value
    assert Binary(36, 29).chirp_mass == pytest.approx(28.1, rel=1e-2)  # GW150914-like


def test_isco_frequency():
    # f_ISCO ~ 4.4 kHz / (M / M_sun)
    assert Binary(1, 1).f_isco == pytest.approx(4397 / 2, rel=1e-3)


@pytest.mark.parametrize("pn_order", [0, 2])
def test_frequency_is_derivative_of_phase(pn_order):
    """2 pi f(t) must equal dPhi/dt -- checks the PN coefficients against each other."""
    b = Binary(1.4, 1.2)
    t = np.linspace(-10.0, -0.5, 200001)
    phase = gwave.gw_phase(b, t, pn_order=pn_order)
    f = gwave.gw_frequency(b, t, pn_order=pn_order)
    dphase = np.gradient(phase, t)
    np.testing.assert_allclose(dphase[1:-1], 2 * np.pi * f[1:-1], rtol=1e-6)


def test_newtonian_time_to_coalescence_inverts_frequency():
    b = Binary(10, 10)
    for f in [10.0, 30.0, 100.0]:
        t = -b.time_to_coalescence(f)
        assert gwave.gw_frequency(b, t, pn_order=0) == pytest.approx(f, rel=1e-12)


@pytest.mark.parametrize("m1, m2, newtonian, total_2pn", [
    # cycles from 10 Hz to the ISCO; literature values (Blanchet 2014), 2PN = sum of terms to 2PN
    (1.4, 1.4, 16031, 16031 + 439 - 208 + 9.9),
    (1.4, 10, 3576, 3576 + 213 - 181 + 9.8),
    (10, 10, 601, 601 + 59 - 51 + 4.1),
])
def test_cycle_counts_match_literature(m1, m2, newtonian, total_2pn):
    b = Binary(m1, m2)
    assert gwave.gw_cycles(b, 10.0, pn_order=0)[1][-1] == pytest.approx(newtonian, rel=2e-3)
    assert gwave.gw_cycles(b, 10.0, pn_order=2)[1][-1] == pytest.approx(total_2pn, rel=2e-3)


def test_antenna_pattern_rotation_invariant():
    """F_+^2 + F_x^2 cannot depend on the polarization angle psi."""
    th, ph = 0.6, 0.4
    norms = [sum(F**2 for F in gwave.antenna_pattern(th, ph, psi)) for psi in np.linspace(0, np.pi, 7)]
    np.testing.assert_allclose(norms, norms[0], rtol=1e-12)


def test_antenna_pattern_overhead():
    fp, fx = gwave.antenna_pattern(0.0, 0.0, 0.0)
    assert (fp, fx) == pytest.approx((1.0, 0.0))


def test_strain_amplitude_scale():
    """GW170817-like binary at 40 Mpc should give a peak strain of order 1e-22 -- 1e-21."""
    t, h = gwave.inspiral(Binary(1.46, 1.27), distance_mpc=40.0, f_low=30.0)
    assert 1e-22 < np.abs(h).max() < 1e-20
    assert t[-1] < 0
