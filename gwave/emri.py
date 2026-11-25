"""
Extreme-mass-ratio inspirals (EMRIs): a compact object of mass mu on a quasi-circular orbit
around a spinning supermassive black hole of mass M >> mu, in the "analytic kludge" model of
L. Barack and C. Cutler, Phys. Rev. D 69, 082005 (2004), restricted to circular orbits as in
the dissertation (eqs. 4.46-4.48). The BH spin is chi = S / M^2 (G = c = 1), and lambda is
the fixed angle between the orbital angular momentum L and the spin S.

Evolved quantities (all functions of time):
    nu       orbital frequency, with post-Newtonian radiation reaction and spin-orbit coupling
    Phi      orbital phase,                         dPhi/dt = 2 pi nu
    gamma    periapsis angle, measured in the orbital plane from L x S
             (relativistic periapsis advance plus Lense-Thirring term)
    alpha    azimuth of L around S (Lense-Thirring precession of the orbital plane)

The wave is the quadrupole field of the orbit, h_ij = (2/D) d^2(mu r_i r_j)/dt^2, built from
vectors and projected transverse-traceless to the line of sight. Projecting onto LISA's
detector tensors then needs no polarization-angle formulas.
"""
from dataclasses import dataclass

import numpy as np
from scipy.integrate import solve_ivp

from .constants import AU, C, MPC, T_SUN, YEAR
from .lisa import lisa_detector_tensors, unit_vector


@dataclass(frozen=True)
class EMRI:
    mu: float            # compact object mass [M_sun]
    M: float             # black-hole mass [M_sun]
    chi: float = 0.0     # dimensionless spin S / M^2
    lam: float = 0.0     # angle between L and S [rad]


def rates(nu, emri, z=0.0):
    """(dnu/dt, dPhi/dt, dgamma/dt, dalpha/dt) at orbital frequency nu (observer frame)."""
    M = emri.M * (1 + z) * T_SUN
    mu = emri.mu * (1 + z) * T_SUN
    v3 = 2 * np.pi * M * nu                     # (v/c)^3
    x = v3**(2 / 3)                             # (v/c)^2
    cl = np.cos(emri.lam)
    dnu = 96 / (10 * np.pi) * mu / M**3 * v3**(11 / 3) * (1 + 1273 / 336 * x - 73 / 12 * v3 * emri.chi * cl)
    dphi = 2 * np.pi * nu
    dgamma = 6 * np.pi * nu * x * (1 + 6.5 * x) - 12 * np.pi * nu * emri.chi * cl * v3
    dalpha = 4 * np.pi * nu * emri.chi * v3
    return dnu, dphi, dgamma, dalpha


def evolve(emri, t, nu0, z=0.0, phi0=0.0, gamma0=0.0, alpha0=0.0):
    """
    Integrate the orbit on the sample times t (starting at t[0] with orbital frequency nu0).
    Stops at the Schwarzschild ISCO, (2 pi M nu)^(2/3) = 1/6. Returns a dict of arrays,
    truncated to the samples before plunge.
    """
    M = emri.M * (1 + z) * T_SUN

    def rhs(_, y):
        return rates(y[0], emri, z)

    def plunge(_, y):
        return (2 * np.pi * M * y[0])**(2 / 3) - 1 / 6
    plunge.terminal = True

    sol = solve_ivp(rhs, (t[0], t[-1]), [nu0, phi0, gamma0, alpha0], t_eval=t, events=plunge,
                    method="DOP853", rtol=1e-11, atol=1e-12)
    n = sol.y.shape[1]
    return {"t": t[:n], "nu": sol.y[0], "Phi": sol.y[1], "gamma": sol.y[2], "alpha": sol.y[3]}


def _perp_basis(v):
    """Two unit vectors orthogonal to unit vector v (and to each other)."""
    a = np.array([1.0, 0.0, 0.0]) if abs(v[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    e1 = np.cross(v, a)
    e1 /= np.linalg.norm(e1)
    return e1, np.cross(v, e1)


def orbital_frame(S_hat, lam, alpha):
    """
    L_hat(t) precessing about S_hat at fixed angle lam, and an orthonormal basis (e1, e2)
    of the orbital plane with e1 along L x S (the reference direction for gamma).
    """
    s1, s2 = _perp_basis(S_hat)
    L = np.cos(lam) * S_hat + np.sin(lam) * (np.cos(alpha)[:, None] * s1 + np.sin(alpha)[:, None] * s2)
    if np.sin(lam) > 1e-12:
        e1 = np.cross(L, S_hat)
        e1 /= np.linalg.norm(e1, axis=-1, keepdims=True)
    else:                                   # orbit in the spin's equatorial plane: no node
        e1 = np.broadcast_to(s1, L.shape)
    e2 = np.cross(L, e1)
    return L, e1, e2


def quadrupole_wave(n_hat, e1, e2, phase, amplitude):
    """
    Transverse-traceless h_ij (shape (N, 3, 3)) of a circular orbit with in-plane basis
    (e1, e2) and orbital angle `phase` measured from e1, seen along n_hat:
        h_ij = -2 A [cos 2phase (e1 e1 - e2 e2) + sin 2phase (e1 e2 + e2 e1)],  then TT.
    `amplitude` is A = 2 mu (2 pi M nu)^(2/3) / D.
    """
    c, s = np.cos(2 * phase)[:, None, None], np.sin(2 * phase)[:, None, None]
    o = lambda a, b: np.einsum("ni,nj->nij", a, b)
    h = -2 * amplitude[:, None, None] * (c * (o(e1, e1) - o(e2, e2)) + s * (o(e1, e2) + o(e2, e1)))
    P = np.eye(3) - np.outer(n_hat, n_hat)
    hP = P @ h @ P
    return hP - 0.5 * P * np.trace(hP, axis1=1, axis2=2)[:, None, None]


def lisa_emri(emri, t, nu0, distance_mpc, theta_s, phi_s, theta_k, phi_k, z=0.0,
              doppler=True, envelope=False):
    """
    LISA strains (h_I, h_II) for an EMRI. (theta_s, phi_s) is the source position and
    (theta_k, phi_k) the BH spin direction, both in ecliptic coordinates. Masses are source
    frame; give the luminosity distance for redshift z. With envelope=True also return the
    slowly varying amplitudes |h_I|, |h_II| (useful for plotting long signals).
    """
    orbit = evolve(emri, np.asarray(t, float), nu0, z)
    t = orbit["t"]
    n_hat = unit_vector(theta_s, phi_s)
    S_hat = unit_vector(theta_k, phi_k)
    L, e1, e2 = orbital_frame(S_hat, emri.lam, orbit["alpha"])

    M = emri.M * (1 + z) * T_SUN
    mu = emri.mu * (1 + z) * T_SUN
    D = distance_mpc * MPC / C
    A = 2 * mu * (2 * np.pi * M * orbit["nu"])**(2 / 3) / D
    phase = orbit["Phi"] + orbit["gamma"]
    if doppler:   # GW phase 2*phase gains 2 pi f_gw (R/c) sin(theta_s) cos(orbital phase - phi_s)
        phase = phase + 2 * np.pi * orbit["nu"] * (AU / C) * np.sin(theta_s) * np.cos(2 * np.pi * t / YEAR - phi_s)

    d1, d2 = lisa_detector_tensors(t)
    h = quadrupole_wave(n_hat, e1, e2, phase, A)
    out = [np.einsum("nij,nij->n", d1, h), np.einsum("nij,nij->n", d2, h)]
    if envelope:
        hq = quadrupole_wave(n_hat, e1, e2, phase + np.pi / 4, A)      # quadrature partner
        out += [np.hypot(out[0], np.einsum("nij,nij->n", d1, hq)),
                np.hypot(out[1], np.einsum("nij,nij->n", d2, hq))]
    return (t, *out)
