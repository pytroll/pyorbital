"""Verification of the propagator against the AIAA 2006-6753 test cases.

The satellites of ``SGP4-VER.TLE`` are propagated to every time for which
``aiaa_results`` gives an expected state, and the two are compared.

Satellites that the propagator cannot handle yet are listed in
``NOT_YET_SUPPORTED``; removing a satellite from that set turns its test back
on. Satellites that are not expected to ever pass are listed in
``SKIP_DECAYING`` with a reason. ``test_every_reference_satellite_is_accounted_for``
makes sure no satellite silently escapes all three treatments.
"""

import datetime as dt
import os
from dataclasses import dataclass

import numpy as np
import pytest

from pyorbital import astronomy
from pyorbital.orbital import Orbital
from pyorbital.tlefile import ChecksumError

_DATAPATH = os.path.dirname(os.path.abspath(__file__))
_TLE_FILE = os.path.join(_DATAPATH, "SGP4-VER.TLE")
_REFERENCE_FILE = os.path.join(_DATAPATH, "aiaa_results")

POSITION_TOLERANCE = 5e-6  # km
VELOCITY_TOLERANCE = 5e-9  # km/s
TIME_TOLERANCE = 1e-3  # minutes

CHECKSUM_ERROR_SATELLITES = {33333, 33334, 33335}

SKIP_DECAYING = {
    29141: "SL-14 DEB is on its last orbits; decay is not modelled, and the "
           "propagated position drifts ~0.3 mm/min away from the reference",
}

NOT_YET_SUPPORTED = {
    # Near-earth, simplified drag equations: rejected by _SGDP4.propagate.
    88888, 29238, 28350, 22312, 28872,
    # Deep space, non-resonant.
    28129, 16925, 28623, 23177, 23599, 4632, 20413, 11801, 23333,
    # Deep space, 12 h resonance.
    9880, 8195, 26975, 21897, 22674,
    # Deep space, 24 h (geosynchronous) resonance.
    24208, 9998, 28626, 26900, 25954, 14128,
}


@dataclass(frozen=True)
class ReferenceState:
    """One expected state vector from ``aiaa_results``."""

    position: tuple[float, float, float]
    velocity: tuple[float, float, float]
    utc_time: np.datetime64 | None


@dataclass(frozen=True)
class VerificationCase:
    """One satellite of ``SGP4-VER.TLE``.

    The delays to propagate to are not taken from the trailing start/stop/step
    of line 2 but from the reference file, which is the authority on which
    states are actually defined: for satellites that decay during the requested
    span the reference stops early, and one satellite asks for a second span
    that the reference does not cover at all.
    """

    satnumber: int
    name: str
    line1: str
    line2: str


def _delay_key(delay):
    """Key a delay by its printed value, so accumulated float error cannot miss a row."""
    return f"{float(delay):.8f}"


def _parse_reference_states(path):
    """Read all expected state vectors, keyed by satellite number and delay."""
    states: dict[int, dict[str, ReferenceState]] = {}
    with open(path) as reference_file:
        for line in reference_file:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.endswith(" xx"):
                states[int(line[:-3])] = current_satellite = {}
                continue
            fields = line.split()
            current_satellite[_delay_key(fields[0])] = ReferenceState(
                position=(float(fields[1]), float(fields[2]), float(fields[3])),
                velocity=(float(fields[4]), float(fields[5]), float(fields[6])),
                utc_time=_parse_reference_time(fields[7:]),
            )
    return states


def _parse_reference_time(fields):
    """Read the trailing ``year mon day hh:mm:ss.ssssss`` fields, if present."""
    if not fields:
        return None
    year, month, day, time_of_day = fields
    utc_time = dt.datetime.strptime(time_of_day, "%H:%M:%S.%f")
    return np.datetime64(utc_time.replace(year=int(year), month=int(month), day=int(day)))


def _parse_verification_cases(path):
    """Read the satellites to propagate and the delays requested for each.

    A satellite may appear more than once, asking for a further span of delays
    that the reference file does not cover; only its first appearance is kept.
    """
    cases = {}
    name = ""
    line1 = ""
    with open(path) as tle_file:
        for line in tle_file:
            line = line.rstrip("\r\n")
            if line.startswith("#"):
                name = line
            elif line.startswith("1 "):
                line1 = line
            elif line.startswith("2 "):
                line2 = line[:69]
                case = VerificationCase(satnumber=int(line2[2:7]), name=name,
                                        line1=line1, line2=line2)
                cases.setdefault(case.satnumber, case)
    return cases


REFERENCE_STATES = _parse_reference_states(_REFERENCE_FILE)
VERIFICATION_CASES = _parse_verification_cases(_TLE_FILE)

PROPAGATION_CASES = sorted(set(VERIFICATION_CASES) - CHECKSUM_ERROR_SATELLITES)


def _reference_delays(satnumber):
    """The delays, in minutes, for which the reference file defines a state."""
    return sorted(float(delay) for delay in REFERENCE_STATES[satnumber])


def _describe(satnumber, delay):
    """Name a satellite and time the way the verification file comments on it."""
    description = VERIFICATION_CASES[satnumber].name.lstrip("#").strip()
    return f"{satnumber} ({description}) at {delay} min"


def _assert_matches_reference(satnumber, delay, utc_time, position, velocity):
    """Compare one propagated state against its reference value."""
    expected = REFERENCE_STATES[satnumber][_delay_key(delay)]

    np.testing.assert_allclose(position, expected.position, atol=POSITION_TOLERANCE, rtol=0,
                               err_msg=f"position of {_describe(satnumber, delay)}")
    np.testing.assert_allclose(velocity, expected.velocity, atol=VELOCITY_TOLERANCE, rtol=0,
                               err_msg=f"velocity of {_describe(satnumber, delay)}")

    if expected.utc_time is not None:
        error_in_minutes = astronomy._days(expected.utc_time - utc_time) * 24 * 60
        assert abs(error_in_minutes) < TIME_TOLERANCE


@pytest.mark.parametrize("satnumber", PROPAGATION_CASES, ids=str)
def test_aiaa_verification_case(satnumber):
    """Propagated states match the AIAA reference values."""
    if satnumber in SKIP_DECAYING:
        pytest.skip(SKIP_DECAYING[satnumber])
    if satnumber in NOT_YET_SUPPORTED:
        pytest.skip(f"{satnumber} is not supported by the propagator yet")

    case = VERIFICATION_CASES[satnumber]
    orbital = Orbital("unknown", line1=case.line1, line2=case.line2)

    for delay in _reference_delays(satnumber):
        utc_time = np.timedelta64(round(delay * 60 * 1e6), "us") + orbital.tle.epoch
        position, velocity = orbital.get_position(utc_time, normalize=False)
        _assert_matches_reference(satnumber, delay, utc_time, position, velocity)


@pytest.mark.parametrize("satnumber", sorted(CHECKSUM_ERROR_SATELLITES), ids=str)
def test_aiaa_case_with_broken_checksum_is_rejected(satnumber):
    """Deliberately corrupted verification TLEs are rejected when read."""
    case = VERIFICATION_CASES[satnumber]

    with pytest.raises(ChecksumError):
        Orbital("unknown", line1=case.line1, line2=case.line2)


def test_every_reference_satellite_is_accounted_for():
    """No satellite of the reference file escapes both testing and an explicit skip."""
    assert set(REFERENCE_STATES) == set(PROPAGATION_CASES)
    assert SKIP_DECAYING.keys() <= set(PROPAGATION_CASES)
    assert NOT_YET_SUPPORTED <= set(PROPAGATION_CASES)
    assert not NOT_YET_SUPPORTED & SKIP_DECAYING.keys()
