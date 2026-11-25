"""
gwave: gravitational-wave signals from inspiralling compact binaries.

    gwave.inspiral         post-Newtonian (0PN / 2PN) waveforms and the L-shaped detector response
    gwave.cosmology        luminosity distance and the massive-graviton distance D(z)
    gwave.lisa             simplified LISA response with the orbital Doppler phase
    gwave.massive_gravity  propagation with a massive graviton
"""
from .constants import AU, C, MPC, T_SUN, YEAR
from .inspiral import (Binary, antenna_pattern, gw_cycles, gw_frequency, gw_phase, inspiral,
                       inspiral_times, polarizations, strain, theta)
from . import cosmology, lisa, massive_gravity
