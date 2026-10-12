"""Propagation of deep-space orbits the AIAA verification set does not reach.

The article's satellites leave parts of the deep-space model untried. None of
its five twelve-hour satellites has an eccentricity between 0.65 and 0.66, where
one of the fits of the resonance to eccentricity changes over; none has a mean
motion just below the band where that resonance is held to apply; and all six of
its geosynchronous satellites are so nearly circular that the eccentricity terms
of their own fit contribute a millionth of what they could. Altering any of
those three places in the propagator leaves the AIAA tests passing.

A seventh day is asked for as well as the first, because the beat between a
resonant orbit and the Earth's rotation is integrated forward from the epoch and
so builds up with time. A difference in the constants that describe it can be
too small to show within a day and plain within a week.

The orbits below sit in each of those gaps. Their expected positions come from
the ``sgp4`` package, a port of the reference implementation the article
describes, which reproduces every state of ``aiaa_results`` to within five
nanometres. That package is not needed to run these tests, only to have written
them.
"""

import numpy as np
import pytest

from pyorbital.orbital import Orbital

POSITION_TOLERANCE = 5e-6  # km, as for the AIAA cases

CASES = {
    "eccentricity above the fit boundary at 0.65": (
        "1 40001U 06001A   06176.00000000  .00000000  00000+0  00000+0 0  9990",
        "2 40001  63.4000  45.0000 6550000 120.0000  30.0000  1.99500000    16",
        (
            (0.0, (-4663.844135811, -13491.263180781, -12429.628762728)),
            (360.0, (29673.513675195, 8151.069052328, -30392.205431881)),
            (1440.0, (-5159.103677060, -13242.779329675, -11423.226602012)),
            (-720.0, (-4411.288916872, -13600.330410202, -12919.198839071)),
            (10000.0, (8065.030788060, 12257.077967201, 6233.620791191)),
        ),
    ),
    "mean motion just below the twelve hour resonance band": (
        "1 40002U 06001A   06176.00000000  .00000000  00000+0  00000+0 0  9991",
        "2 40002  63.4000  45.0000 6000000 120.0000  30.0000  1.88570000    12",
        (
            (0.0, (-7026.383636164, -14550.621581078, -10590.893711190)),
            (360.0, (29517.980694216, 7079.922133772, -31691.360946876)),
            (1440.0, (-693.801125596, 6913.594467648, 10737.220017617)),
            (-720.0, (-1159.549983466, -15907.214392140, -20775.504005053)),
            (10000.0, (2850.484294481, -15296.203768858, -25754.090168744)),
        ),
    ),
    "geosynchronous and eccentric enough to feel its own fit": (
        "1 40003U 06001A   06176.00000000  .00000000  00000+0  00000+0 0  9992",
        "2 40003  10.0000  45.0000 6000000 120.0000  30.0000  1.00270000    12",
        (
            (0.0, (-3908.635579786, -29105.396115816, -3167.620339585)),
            (360.0, (50157.215190843, -34792.143566207, -10620.591900083)),
            (1440.0, (-3138.053160242, -29721.818498590, -3341.469673313)),
            (-720.0, (65153.254340756, -5998.684734843, -8873.889192963)),
            (10000.0, (-13068.505240716, -18046.712453215, -640.852992890)),
        ),
    ),
}


@pytest.mark.parametrize("description", sorted(CASES))
def test_deep_space_orbit_in_a_gap_of_the_verification_set(description):
    """The propagated states match the reference implementation."""
    line1, line2, expected_states = CASES[description]
    orbital = Orbital("unknown", line1=line1, line2=line2)

    for delay, expected_position in expected_states:
        utc_time = np.timedelta64(round(delay * 60 * 1e6), "us") + orbital.tle.epoch
        position, _ = orbital.get_position(utc_time, normalize=False)

        np.testing.assert_allclose(position, expected_position, rtol=0,
                                   atol=POSITION_TOLERANCE,
                                   err_msg=f"{description} at {delay} min")
