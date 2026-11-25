"""Tests for the EMRI analytic-kludge model and the LISA frame it uses."""
import numpy as np
import pytest

import gwave
from gwave import emri as em
from gwave.constants import T_SUN, YEAR
from gwave.lisa import lisa_frame, unit_vector


def test_lisa_frame_matches_cutler():
    """Source position in the detector frame, Cutler (1998) eqs. for theta(t), phi(t)."""
    t = np.linspace(0, YEAR, 1001)
    ths, phs = 0.7, 1.9
    n = unit_vector(ths, phs)
    x, y, zd = lisa_frame(t)
    pb = 2 * np.pi * t / YEAR
    np.testing.assert_allclose(zd @ n, 0.5 * np.cos(ths) - np.sqrt(3) / 2 * np.sin(ths) * np.cos(pb - phs), atol=1e-12)
    with np.errstate(divide="ignore"):
        phi_c = pb + np.arctan((np.sqrt(3) * np.cos(ths) + np.sin(ths) * np.cos(pb - phs))
                               / (2 * np.sin(ths) * np.sin(pb - phs)))
    d = np.mod(np.arctan2(y @ n, x @ n) - phi_c + 0.1, np.pi) - 0.1      # arctan is defined mod pi
    np.testing.assert_allclose(d, 0, atol=1e-9)


def test_newtonian_limit_of_frequency_evolution():
    """With no spin and v -> 0, dnu/dt is the quadrupole chirp (f_gw = 2 nu)."""
    e = em.EMRI(mu=10, M=1e6)
    nu = 1e-5
    dnu = em.rates(nu, e)[0]
    b = gwave.Binary(10, 1e6)
    dfgw = 96 / 5 * np.pi**(8 / 3) * b.Mc_s**(5 / 3) * (2 * nu)**(11 / 3)
    assert dnu == pytest.approx(dfgw / 2, rel=2e-3)        # (1PN correction ~ 1e-3 here)


def test_precession_rates():
    e = em.EMRI(mu=1, M=1e6, chi=0.1, lam=np.pi / 6)
    nu = 5e-4
    M = 1e6 * T_SUN
    _, dphi, dgamma, dalpha = em.rates(nu, e)
    assert dalpha == pytest.approx(8 * np.pi**2 * nu**2 * e.chi * M)                       # eq. 4.47
    x = (2 * np.pi * M * nu)**(2 / 3)
    expected = 6 * np.pi * nu * x * (1 + 6.5 * x) - 24 * np.pi**2 * nu**2 * np.cos(e.lam) * e.chi * M
    assert dgamma == pytest.approx(expected)                                                # eq. 4.46
    assert em.rates(nu, em.EMRI(mu=1, M=1e6))[3] == 0


def test_face_on_amplitude():
    """Seen face-on by a detector with arms along x and y: amplitude 4 mu (2 pi M nu)^(2/3) / D."""
    n = np.array([0.0, 0.0, 1.0])
    e1 = np.tile([1.0, 0.0, 0.0], (200, 1))
    e2 = np.tile([0.0, 1.0, 0.0], (200, 1))
    phase = np.linspace(0, 2 * np.pi, 200)
    A = np.full(200, 1.0)
    h = em.quadrupole_wave(n, e1, e2, phase, A)
    D = 0.5 * (np.diag([1, 0, 0]) - np.diag([0, 1, 0]))
    s = np.einsum("ij,nij->n", D, h)
    assert np.max(np.abs(s)) == pytest.approx(2 * A[0])       # A here = 2 mu (..)/D


def test_tensor_response_matches_antenna_pattern():
    """For any geometry the tensor projection equals F_+ h_+ + F_x h_x for some psi."""
    rng = np.random.default_rng(1)
    th, ph = 1.0, 0.6
    n = unit_vector(th, ph)
    L = unit_vector(2.0, -1.3)
    e1 = np.cross(n, L); e1 /= np.linalg.norm(e1)            # line of nodes
    e2 = np.cross(L, e1)
    phase = rng.uniform(0, 2 * np.pi, 50)
    A = np.ones(50)
    h = em.quadrupole_wave(n, np.tile(e1, (50, 1)), np.tile(e2, (50, 1)), phase, A)
    D = 0.5 * (np.diag([1.0, 0, 0]) - np.diag([0, 1.0, 0]))
    s = np.einsum("ij,nij->n", D, h)
    # textbook polarizations for inclination iota, phase measured from the line of nodes
    ci = L @ n
    hp = -A * (1 + ci**2) * np.cos(2 * phase)
    hx = -2 * A * ci * np.sin(2 * phase)
    psis = np.linspace(0, np.pi, 3601)
    errs = [min(np.max(np.abs(s - (fp * hp + sg * fx * hx))) for sg in (1, -1))
            for fp, fx in zip(*gwave.antenna_pattern(th, ph, psis))]
    assert min(errs) < 1e-3


def test_spin_precession_modulates_amplitude():
    """
    Thesis Figs. 4.18-4.22: Lense-Thirring precession of the orbital plane modulates the
    amplitude, increasingly with spin. Dividing by the spinless envelope removes the slow
    modulation from LISA's own motion.
    """
    t = np.arange(0, 5 * 86400, 20.0)
    kw = dict(nu0=5e-4, distance_mpc=gwave.cosmology.luminosity_distance(1.0), theta_s=np.pi / 4,
              phi_s=np.pi / 4, theta_k=np.pi / 3, phi_k=np.pi / 4, z=1.0, envelope=True)
    env = {chi: em.lisa_emri(em.EMRI(1, 1e6, chi=chi, lam=np.pi / 6), t, **kw)[3]
           for chi in (0.0, 1e-4, 0.1)}
    spread = lambda e: e.max() - e.min()
    assert spread(env[1e-4] / env[0.0]) < 0.01
    assert spread(env[0.1] / env[0.0]) > 0.1


def test_plunge_stops_evolution():
    e = em.EMRI(mu=10, M=1e5)
    M = 1e5 * T_SUN
    nu_isco = 6**-1.5 / (2 * np.pi * M)
    orbit = em.evolve(e, np.linspace(0, 1e7, 20001), 0.8 * nu_isco)
    assert len(orbit["t"]) < 20001
    assert orbit["nu"][-1] <= nu_isco * 1.001
