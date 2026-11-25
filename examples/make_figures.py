"""Generate the figures in figures/ used by the README.  Run from the repository root:

    python examples/make_figures.py
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import gwave
from gwave import Binary

OUT = os.path.join(os.path.dirname(__file__), "..", "figures")
os.makedirs(OUT, exist_ok=True)

# palette shared with asinha5's website
PAPER, INK, MUTED, LINE = "#f3eee3", "#1c2130", "#6a6556", "#d8cfbd"
ORANGE, TEAL, OCHRE, BRICK = "#c4561d", "#2b7a6f", "#b5862a", "#b8452f"

plt.rcParams.update({
    "figure.facecolor": PAPER, "axes.facecolor": PAPER, "savefig.facecolor": PAPER,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.grid": True,
    "grid.color": LINE, "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
    "font.family": "DejaVu Sans", "font.size": 10.5, "axes.titlesize": 11.5, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "legend.frameon": False, "lines.linewidth": 1.4,
    "mathtext.fontset": "dejavusans",
})


def save(fig, name):
    path = os.path.join(OUT, name)
    fig.savefig(path, dpi=160, bbox_inches="tight")
    plt.close(fig)
    print("wrote", os.path.relpath(path))


# 1. the chirp: Newtonian vs 2PN for a heavy binary black hole -----------------------------
def fig_chirp():
    b = Binary(30, 30)
    f_low, fs = 20.0, 16384.0
    t, h2 = gwave.inspiral(b, distance_mpc=400, f_low=f_low, sample_rate=fs, pn_order=2)
    t0, h0 = gwave.inspiral(b, distance_mpc=400, f_low=f_low, sample_rate=fs, pn_order=0)

    fig, (ax, axf) = plt.subplots(2, 1, figsize=(8.2, 5.4), sharex=True,
                                  gridspec_kw={"height_ratios": [2.2, 1], "hspace": 0.12})
    ax.plot(t0, h0 * 1e21, color=TEAL, alpha=0.75, label="Newtonian (0PN) phase")
    ax.plot(t, h2 * 1e21, color=ORANGE, label="2PN phase")
    ax.set_ylabel(r"strain $h\ [10^{-21}]$")
    ax.set_title("Inspiral of a 30 + 30 $M_\\odot$ binary black hole at 400 Mpc, from 20 Hz to the ISCO")
    ax.legend(loc="upper left", ncol=2)

    axf.plot(t0, gwave.gw_frequency(b, t0, pn_order=0), color=TEAL, alpha=0.75)
    axf.plot(t, gwave.gw_frequency(b, t, pn_order=2), color=ORANGE)
    axf.axhline(b.f_isco, color=MUTED, ls=(0, (3, 3)), lw=1)
    axf.text(t0[0], b.f_isco * 1.04, "  ISCO", color=MUTED, fontsize=9, va="bottom")
    axf.set_ylabel("GW frequency [Hz]")
    axf.set_xlabel("time before coalescence [s]")
    save(fig, "chirp.png")


# 2. how wrong is the Newtonian phase? ------------------------------------------------------
def fig_dephasing():
    systems = [("BNS  1.4 + 1.4", Binary(1.4, 1.4), ORANGE),
               ("NSBH 1.4 + 10", Binary(1.4, 10), OCHRE),
               ("BBH  10 + 10", Binary(10, 10), TEAL)]
    fig, ax = plt.subplots(figsize=(8.2, 4.2))
    for label, b, col in systems:
        f, n0 = gwave.gw_cycles(b, 10.0, pn_order=0)
        _, n2 = gwave.gw_cycles(b, 10.0, pn_order=2)
        ax.plot(f, n0 - n2, color=col, label=f"{label} $M_\\odot$   (N$_{{2PN}}$ = {n2[-1]:,.0f} cycles)")
    ax.set_xscale("log")
    ax.set_xlabel("GW frequency [Hz]")
    ax.set_ylabel(r"$N_{\rm Newt} - N_{\rm 2PN}$  [cycles]")
    ax.set_title("Phase error of the Newtonian approximation, in GW cycles counted from 10 Hz")
    ax.legend(loc="upper left")
    save(fig, "dephasing.png")


# 3. where on the sky an L-shaped detector is deaf -----------------------------------------
def fig_antenna():
    lon = np.linspace(-np.pi, np.pi, 361)
    lat = np.linspace(-np.pi / 2, np.pi / 2, 181)
    LON, LAT = np.meshgrid(lon, lat)
    fp, fx = gwave.antenna_pattern(np.pi / 2 - LAT, LON, 0.0)
    maps = [(fp, r"$F_+$  ($\psi = 0$)", "RdBu_r"), (fx, r"$F_\times$  ($\psi = 0$)", "RdBu_r"),
            (np.hypot(fp, fx), r"$\sqrt{F_+^2 + F_\times^2}$  (any $\psi$)", "magma")]
    fig = plt.figure(figsize=(12, 3.0))
    for i, (z, title, cmap) in enumerate(maps):
        ax = fig.add_subplot(1, 3, i + 1, projection="mollweide")
        ax.set_facecolor(PAPER)
        vmin = 0 if cmap == "magma" else -1
        im = ax.pcolormesh(LON, LAT, z, cmap=cmap, vmin=vmin, vmax=1, shading="auto", rasterized=True)
        ax.set_title(title, loc="center", pad=10)
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.grid(color="white", alpha=0.35, lw=0.5)
        fig.colorbar(im, ax=ax, orientation="horizontal", shrink=0.6, pad=0.06)
    fig.suptitle("Antenna pattern of an interferometer with arms along x and y (detector frame)",
                 x=0.06, y=1.02, ha="left", fontweight="bold", fontsize=11.5)
    fig.subplots_adjust(top=0.98, wspace=0.08)
    save(fig, "antenna_pattern.png")


# 4. a massive graviton compresses the chirp ----------------------------------------------
def fig_massive_graviton():
    from gwave import massive_gravity as mg

    fig, (ax, axp) = plt.subplots(1, 2, figsize=(12, 4.2), gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.22})

    b, z, lam = Binary(10, 10), 0.1, 1e16
    t, h_gr, h_m = mg.massive_inspiral(b, z, lam, f_low=30, sample_rate=8192)
    sel = (t > -1.2) & (t < -0.02)
    ax.plot(t[sel], h_gr[sel] * 1e22, color=TEAL, alpha=0.8, label="GR (massless graviton)")
    ax.plot(t[sel], h_m[sel] * 1e22, color=ORANGE, label=f"$\\lambda_g$ = {lam:.0e} m")
    ax.set_xlabel("time [s]")
    ax.set_ylabel(r"strain $h\ [10^{-22}]$")
    ax.set_title("A massive graviton compresses the chirp  (10 + 10 $M_\\odot$, z = 0.1)")
    ax.legend(loc="upper left", ncol=2)
    ax.text(0.99, 0.03, f"$m_g$ = {mg.graviton_mass_ev(lam):.1e} eV, already excluded by LIGO",
            transform=ax.transAxes, ha="right", fontsize=9, color=MUTED)

    f = np.geomspace(1e-5, 3e3, 400)
    for lam_, col in [(1e16, BRICK), (1e17, ORANGE), (1e18, OCHRE)]:
        axp.loglog(f, np.abs(mg.phase_shift(f, 1.0, lam_)), color=col,
                   label=f"$\\lambda_g$ = $10^{{{int(np.log10(lam_))}}}$ m")
    axp.axvspan(1e-4, 1e-1, color=TEAL, alpha=0.10)
    axp.axvspan(10, 1e3, color=ORANGE, alpha=0.10)
    axp.text(3e-3, 3e-2, "LISA", color=TEAL, ha="center")
    axp.text(100, 3e-2, "LIGO", color=ORANGE, ha="center")
    axp.axhline(1, color=MUTED, ls=(0, (3, 3)), lw=1)
    axp.set_ylim(1e-2, 1e16)
    axp.set_xlabel("GW frequency [Hz]")
    axp.set_ylabel(r"$|\Psi_g(f)|$ [rad]")
    axp.set_title("Phase accumulated from z = 1")
    axp.legend(loc="upper right")
    save(fig, "massive_graviton.png")


if __name__ == "__main__":
    fig_chirp()
    fig_dephasing()
    fig_antenna()
    fig_massive_graviton()
