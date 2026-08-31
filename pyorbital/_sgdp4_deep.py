"""Deep-space corrections to the SGP4 propagator, following AIAA 2006-6753.

A satellite whose orbit takes 225 minutes or more is far enough from the Earth
that the Sun and the Moon move it as much as the Earth's own oblateness does,
and that its orbit can beat against the Earth's rotation. SGP4 leaves both out;
this module adds them, in the three steps the article describes:

* the geometry of the Sun and the Moon relative to the orbit at the epoch, and
  the coefficients derived from it, computed once (``dscom``, ``dsinit``);
* the slow drift those two bodies impose, added to the secular rates before
  each propagation, together with the beat against the Earth's rotation for the
  orbits that resonate with it (``dspace``);
* the periodic wobble they impose, added to the orbital elements themselves
  (``dpper``).

The names of the coefficients are the ones used in the article, so that the two
can be read side by side.
"""

import numpy as np

TWO_PI = 2.0 * np.pi

# Amplitude of the solar and lunar perturbations, and the eccentricity of the
# Earth's orbit about the Sun and of the Moon's about the Earth.
SOLAR_PERTURBATION = 2.9864797e-6
LUNAR_PERTURBATION = 4.7968065e-7
SOLAR_ECCENTRICITY = 0.01675
LUNAR_ECCENTRICITY = 0.05490

# The Sun's position in the orbital plane of the Earth, at the reference epoch.
SOLAR_SIN_INCLINATION = 0.39785416
SOLAR_COS_INCLINATION = 0.91744867
SOLAR_COS_ARG_PERIGEE = 0.1945905
SOLAR_SIN_ARG_PERIGEE = -0.98088458

# Days between 1949 December 31 and 1900 January 0, the epoch the lunar
# ephemeris below is written against.
DAYS_FROM_1900 = 18261.5


class DeepSpace:
    """The Sun and the Moon, as they act on one satellite's orbit.

    Built once from the orbit at its epoch, then asked at each propagation for
    the drift and the wobble to add to the mean elements.
    """

    def __init__(self, days_since_1949, eccentricity, arg_perigee, inclination,
                 right_ascension, mean_motion, mean_anomaly, sidereal_time_at_epoch,
                 mean_anomaly_rate, arg_perigee_rate, right_ascension_rate,
                 semi_major_axis):
        """Work out the geometry and the coefficients that follow from it."""
        self._geometry = compute_lunisolar_geometry(
            days_since_1949, eccentricity, arg_perigee, inclination, right_ascension,
            mean_motion)
        self._secular_rates = compute_secular_rates(self._geometry)
        self._mean_motion_at_epoch = mean_motion
        self._arg_perigee_at_epoch = arg_perigee
        self._arg_perigee_rate = arg_perigee_rate
        self.resonance = classify_resonance(mean_motion, eccentricity)
        if self.resonance != NOT_RESONANT:
            self._prepare_resonance(eccentricity, arg_perigee, right_ascension,
                                    mean_anomaly, sidereal_time_at_epoch, mean_motion,
                                    mean_anomaly_rate, arg_perigee_rate,
                                    right_ascension_rate, semi_major_axis)

    def _prepare_resonance(self, eccentricity, arg_perigee, right_ascension, mean_anomaly,
                           sidereal_time_at_epoch, mean_motion, mean_anomaly_rate,
                           arg_perigee_rate, right_ascension_rate, semi_major_axis):
        """Work out how strongly the orbit beats against the Earth, and from where."""
        rates = self._secular_rates
        aonv = 1.0 / semi_major_axis
        theta = sidereal_time_at_epoch

        if self.resonance != RESONANT_TWICE_A_DAY:
            raise NotImplementedError(
                "A satellite that circles once a day is not supported yet")

        self._amplitudes = _twice_daily_resonance(self._geometry, eccentricity,
                                                  mean_motion, aonv)
        self._longitude_at_epoch = np.fmod(
            mean_anomaly + right_ascension + right_ascension - theta - theta, TWO_PI)
        self._longitude_rate_offset = (
            mean_anomaly_rate + rates["dmdt"]
            + 2.0 * (right_ascension_rate + rates["dnodt"] - EARTH_ROTATION)
            - mean_motion)

    def resonant_correction(self, minutes_since_epoch, right_ascension, arg_perigee,
                            sidereal_time_at_epoch):
        """Integrate the beat against the Earth up to this time.

        Returns the mean motion there and the mean anomaly it implies. The
        integration always restarts from the epoch, so that asking for one time
        never depends on which times were asked for before.
        """
        longitude, mean_motion = _integrate_resonance(
            minutes_since_epoch, self._longitude_at_epoch, self._mean_motion_at_epoch,
            self._longitude_rate_offset, self._arg_perigee_at_epoch,
            self._arg_perigee_rate, self.resonance, self._amplitudes)

        # The bulges of the equator the orbit beats against turn with the Earth.
        sidereal_time = np.fmod(
            sidereal_time_at_epoch + minutes_since_epoch * EARTH_ROTATION, TWO_PI)

        mean_anomaly = longitude - 2.0 * right_ascension + 2.0 * sidereal_time
        return mean_motion, mean_anomaly

    def secular_drift(self, minutes_since_epoch):
        """Give the slow change of each element since the epoch."""
        rates = self._secular_rates
        return dict(eccentricity=rates["dedt"] * minutes_since_epoch,
                    inclination=rates["didt"] * minutes_since_epoch,
                    mean_anomaly=rates["dmdt"] * minutes_since_epoch,
                    arg_perigee=rates["domdt"] * minutes_since_epoch,
                    right_ascension=rates["dnodt"] * minutes_since_epoch)

    def periodic_shift(self, minutes_since_epoch):
        """Give the periodic displacement of each element at this time."""
        return compute_periodics(self._geometry, minutes_since_epoch)


# Rate at which the Earth turns under the orbit, radians per minute.
EARTH_ROTATION = 4.37526908801129966e-3

# Amplitudes of the tesseral harmonics that a resonant orbit beats against, and
# the longitudes of the corresponding bulges.
ROOT22, ROOT32, ROOT44, ROOT52, ROOT54 = (1.7891679e-6, 3.7393792e-7, 7.3636953e-9,
                                          1.1428639e-7, 2.1765803e-9)
G22, G32, G44, G52, G54 = 5.7686396, 0.95240898, 1.8014998, 1.0508330, 4.4108898

# The step the resonance is integrated with, in minutes.
RESONANCE_STEP = 720.0

NOT_RESONANT = 0
RESONANT_ONCE_A_DAY = 1
RESONANT_TWICE_A_DAY = 2


def classify_resonance(mean_motion, eccentricity):
    """Say whether the orbit beats against the Earth's rotation, and how.

    An orbit that circles once or twice while the Earth turns once keeps meeting
    the same bulges of the equator, so their small pull accumulates instead of
    averaging away.
    """
    if 0.0034906585 < mean_motion < 0.0052359877:
        return RESONANT_ONCE_A_DAY
    if 8.26e-3 <= mean_motion <= 9.24e-3 and eccentricity >= 0.5:
        return RESONANT_TWICE_A_DAY
    return NOT_RESONANT


def _twice_daily_eccentricity_terms(eccentricity):
    """Fit the tesseral amplitudes to the eccentricity, for a twelve hour orbit."""
    e = eccentricity
    esq = e * e
    ecube = e * esq

    g201 = -0.306 - (e - 0.64) * 0.440
    if e <= 0.65:
        g211 = 3.616 - 13.2470 * e + 16.2900 * esq
        g310 = -19.302 + 117.3900 * e - 228.4190 * esq + 156.5910 * ecube
        g322 = -18.9068 + 109.7927 * e - 214.6334 * esq + 146.5816 * ecube
        g410 = -41.122 + 242.6940 * e - 471.0940 * esq + 313.9530 * ecube
        g422 = -146.407 + 841.8800 * e - 1629.014 * esq + 1083.4350 * ecube
        g520 = -532.114 + 3017.977 * e - 5740.032 * esq + 3708.2760 * ecube
    else:
        g211 = -72.099 + 331.819 * e - 508.738 * esq + 266.724 * ecube
        g310 = -346.844 + 1582.851 * e - 2415.925 * esq + 1246.113 * ecube
        g322 = -342.585 + 1554.908 * e - 2366.899 * esq + 1215.972 * ecube
        g410 = -1052.797 + 4758.686 * e - 7193.992 * esq + 3651.957 * ecube
        g422 = -3581.690 + 16178.110 * e - 24462.770 * esq + 12422.520 * ecube
        if e > 0.715:
            g520 = -5149.66 + 29936.92 * e - 54087.36 * esq + 31324.56 * ecube
        else:
            g520 = 1464.74 - 4664.75 * e + 3763.64 * esq

    if e < 0.7:
        g533 = -919.22770 + 4988.6100 * e - 9064.7700 * esq + 5542.21 * ecube
        g521 = -822.71072 + 4568.6173 * e - 8491.4146 * esq + 5337.524 * ecube
        g532 = -853.66600 + 4690.2500 * e - 8624.7700 * esq + 5341.4 * ecube
    else:
        g533 = -37995.780 + 161616.52 * e - 229838.20 * esq + 109377.94 * ecube
        g521 = -51752.104 + 218913.95 * e - 309468.16 * esq + 146349.42 * ecube
        g532 = -40023.880 + 170470.89 * e - 242699.48 * esq + 115605.82 * ecube

    return g201, g211, g310, g322, g410, g422, g520, g521, g532, g533


def _twice_daily_inclination_terms(sinim, cosim):
    """Fit the tesseral amplitudes to the inclination, for a twelve hour orbit."""
    cosisq = cosim * cosim
    sini2 = sinim * sinim

    f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq)
    f221 = 1.5 * sini2
    f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq)
    f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq)
    f441 = 35.0 * sini2 * f220
    f442 = 39.3750 * sini2 * sini2
    f522 = 9.84375 * sinim * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq)
                              + 0.33333333 * (-2.0 + 4.0 * cosim + 6.0 * cosisq))
    f523 = sinim * (4.92187512 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq)
                    + 6.56250012 * (1.0 + 2.0 * cosim - 3.0 * cosisq))
    f542 = 29.53125 * sinim * (2.0 - 8.0 * cosim
                               + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq))
    f543 = 29.53125 * sinim * (-2.0 - 8.0 * cosim
                               + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq))
    return f220, f221, f321, f322, f441, f442, f522, f523, f542, f543


def _twice_daily_resonance(geometry, eccentricity_at_epoch, mean_motion, aonv):
    """Give the amplitudes with which a twelve hour orbit beats against the Earth."""
    g201, g211, g310, g322, g410, g422, g520, g521, g532, g533 = (
        _twice_daily_eccentricity_terms(eccentricity_at_epoch))
    f220, f221, f321, f322, f441, f442, f522, f523, f542, f543 = (
        _twice_daily_inclination_terms(geometry["sin_inclination"],
                                       geometry["cos_inclination"]))

    temp1 = 3.0 * mean_motion * mean_motion * aonv * aonv
    temp = temp1 * ROOT22
    amplitudes = {"d2201": temp * f220 * g201, "d2211": temp * f221 * g211}
    temp1 = temp1 * aonv
    temp = temp1 * ROOT32
    amplitudes.update(d3210=temp * f321 * g310, d3222=temp * f322 * g322)
    temp1 = temp1 * aonv
    temp = 2.0 * temp1 * ROOT44
    amplitudes.update(d4410=temp * f441 * g410, d4422=temp * f442 * g422)
    temp1 = temp1 * aonv
    temp = temp1 * ROOT52
    amplitudes.update(d5220=temp * f522 * g520, d5232=temp * f523 * g532)
    temp = 2.0 * temp1 * ROOT54
    amplitudes.update(d5421=temp * f542 * g521, d5433=temp * f543 * g533)
    return amplitudes


def _twice_daily_rates(amplitudes, longitude, arg_perigee, mean_longitude_rate):
    """Give how fast a twelve hour orbit's mean longitude and mean motion change."""
    a = amplitudes
    perigee2 = arg_perigee + arg_perigee
    longitude2 = longitude + longitude

    rate = (a["d2201"] * np.sin(perigee2 + longitude - G22)
            + a["d2211"] * np.sin(longitude - G22)
            + a["d3210"] * np.sin(arg_perigee + longitude - G32)
            + a["d3222"] * np.sin(-arg_perigee + longitude - G32)
            + a["d4410"] * np.sin(perigee2 + longitude2 - G44)
            + a["d4422"] * np.sin(longitude2 - G44)
            + a["d5220"] * np.sin(arg_perigee + longitude - G52)
            + a["d5232"] * np.sin(-arg_perigee + longitude - G52)
            + a["d5421"] * np.sin(arg_perigee + longitude2 - G54)
            + a["d5433"] * np.sin(-arg_perigee + longitude2 - G54))
    acceleration = (a["d2201"] * np.cos(perigee2 + longitude - G22)
                    + a["d2211"] * np.cos(longitude - G22)
                    + a["d3210"] * np.cos(arg_perigee + longitude - G32)
                    + a["d3222"] * np.cos(-arg_perigee + longitude - G32)
                    + a["d5220"] * np.cos(arg_perigee + longitude - G52)
                    + a["d5232"] * np.cos(-arg_perigee + longitude - G52)
                    + 2.0 * (a["d4410"] * np.cos(perigee2 + longitude2 - G44)
                             + a["d4422"] * np.cos(longitude2 - G44)
                             + a["d5421"] * np.cos(arg_perigee + longitude2 - G54)
                             + a["d5433"] * np.cos(-arg_perigee + longitude2 - G54)))
    return rate, acceleration * mean_longitude_rate


def _integrate_resonance(minutes_since_epoch, longitude_at_epoch, mean_motion_at_epoch,
                         longitude_rate_offset, arg_perigee_at_epoch, arg_perigee_rate,
                         resonance, amplitudes):
    """March the resonance from the epoch to the wanted time in fixed steps.

    The pull depends on where the orbit has drifted to, so it has to be followed
    step by step rather than evaluated in one go.
    """
    step = RESONANCE_STEP if minutes_since_epoch > 0.0 else -RESONANCE_STEP
    half_step_squared = 0.5 * RESONANCE_STEP * RESONANCE_STEP

    longitude = longitude_at_epoch
    mean_motion = mean_motion_at_epoch
    elapsed = 0.0

    while True:
        rate, acceleration = _resonance_rates(resonance, amplitudes, longitude,
                                              arg_perigee_at_epoch, arg_perigee_rate,
                                              elapsed, mean_motion,
                                              longitude_rate_offset)
        longitude_rate = mean_motion + longitude_rate_offset
        if np.abs(minutes_since_epoch - elapsed) < RESONANCE_STEP:
            break
        longitude += longitude_rate * step + rate * half_step_squared
        mean_motion += rate * step + acceleration * half_step_squared
        elapsed += step

    remaining = minutes_since_epoch - elapsed
    return (longitude + longitude_rate * remaining + rate * remaining * remaining * 0.5,
            mean_motion + rate * remaining + acceleration * remaining * remaining * 0.5)


def _resonance_rates(resonance, amplitudes, longitude, arg_perigee_at_epoch,
                     arg_perigee_rate, elapsed, mean_motion, longitude_rate_offset):
    """Give the rate and the acceleration of the beat at one point of the march."""
    longitude_rate = mean_motion + longitude_rate_offset
    arg_perigee = arg_perigee_at_epoch + arg_perigee_rate * elapsed
    return _twice_daily_rates(amplitudes, longitude, arg_perigee, longitude_rate)


def _lunar_geometry(day):
    """Locate the Moon's orbital plane and perigee on the given day.

    The Moon's node regresses once around in about 18.6 years and its orbit is
    inclined to the ecliptic, so unlike the Sun's its geometry cannot be taken
    as fixed.
    """
    node = np.fmod(4.5236020 - 9.2422029e-4 * day, TWO_PI)
    sin_node, cos_node = np.sin(node), np.cos(node)

    cos_inclination = 0.91375164 - 0.03568096 * cos_node
    sin_inclination = np.sqrt(1.0 - cos_inclination * cos_inclination)
    sin_node_ecliptic = 0.089683511 * sin_node / sin_inclination
    cos_node_ecliptic = np.sqrt(1.0 - sin_node_ecliptic * sin_node_ecliptic)

    mean_longitude = 5.8351514 + 0.0019443680 * day
    argument = np.arctan2(
        SOLAR_SIN_INCLINATION * sin_node / sin_inclination,
        cos_node_ecliptic * cos_node + SOLAR_COS_INCLINATION * sin_node_ecliptic * sin_node,
    )
    argument = mean_longitude + argument - node

    return dict(sin_inclination=sin_inclination, cos_inclination=cos_inclination,
                sin_node=sin_node_ecliptic, cos_node=cos_node_ecliptic,
                sin_arg_perigee=np.sin(argument), cos_arg_perigee=np.cos(argument),
                mean_longitude=mean_longitude)


def _third_body_terms(perturbing, orbit):
    """Project one perturbing body onto the satellite's orbit.

    Returns the ``s`` and ``z`` coefficients of the article, which describe how
    that body's pull varies around the orbit. The same projection serves the
    Sun and the Moon; only the direction of the perturbing body differs.
    """
    zcosg, zsing = perturbing["cos_arg_perigee"], perturbing["sin_arg_perigee"]
    zcosi, zsini = perturbing["cos_inclination"], perturbing["sin_inclination"]
    zcosh, zsinh = perturbing["cos_node"], perturbing["sin_node"]

    cosim, sinim = orbit["cos_inclination"], orbit["sin_inclination"]
    cosomm, sinomm = orbit["cos_arg_perigee"], orbit["sin_arg_perigee"]
    emsq, betasq, rtemsq = orbit["emsq"], orbit["betasq"], orbit["rtemsq"]

    a1 = zcosg * zcosh + zsing * zcosi * zsinh
    a3 = -zsing * zcosh + zcosg * zcosi * zsinh
    a7 = -zcosg * zsinh + zsing * zcosi * zcosh
    a8 = zsing * zsini
    a9 = zsing * zsinh + zcosg * zcosi * zcosh
    a10 = zcosg * zsini
    a2 = cosim * a7 + sinim * a8
    a4 = cosim * a9 + sinim * a10
    a5 = -sinim * a7 + cosim * a8
    a6 = -sinim * a9 + cosim * a10

    x1 = a1 * cosomm + a2 * sinomm
    x2 = a3 * cosomm + a4 * sinomm
    x3 = -a1 * sinomm + a2 * cosomm
    x4 = -a3 * sinomm + a4 * cosomm
    x5 = a5 * sinomm
    x6 = a6 * sinomm
    x7 = a5 * cosomm
    x8 = a6 * cosomm

    z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3
    z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4
    z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4

    z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq
    z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq
    z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq
    z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5)
    z12 = (-6.0 * (a1 * a6 + a3 * a5)
           + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5)))
    z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6)
    z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7)
    z22 = (6.0 * (a4 * a5 + a2 * a6)
           + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8)))
    z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8)
    z1 = z1 + z1 + betasq * z31
    z2 = z2 + z2 + betasq * z32
    z3 = z3 + z3 + betasq * z33

    s3 = perturbing["amplitude"] / orbit["mean_motion"]
    s2 = -0.5 * s3 / rtemsq
    s4 = s3 * rtemsq
    s1 = -15.0 * orbit["eccentricity"] * s4
    s5 = x1 * x3 + x2 * x4
    s6 = x2 * x3 + x1 * x4
    s7 = x2 * x4 - x1 * x3

    return dict(s1=s1, s2=s2, s3=s3, s4=s4, s5=s5, s6=s6, s7=s7,
                z1=z1, z2=z2, z3=z3, z11=z11, z12=z12, z13=z13,
                z21=z21, z22=z22, z23=z23, z31=z31, z32=z32, z33=z33)


def _solar_geometry(sin_node, cos_node):
    """Locate the Sun's orbital plane and perigee, which barely move."""
    return dict(sin_inclination=SOLAR_SIN_INCLINATION, cos_inclination=SOLAR_COS_INCLINATION,
                sin_arg_perigee=SOLAR_SIN_ARG_PERIGEE, cos_arg_perigee=SOLAR_COS_ARG_PERIGEE,
                sin_node=sin_node, cos_node=cos_node, amplitude=SOLAR_PERTURBATION)


def _moon_relative_to_orbit(moon, sin_node, cos_node):
    """Express the Moon's node relative to the satellite's node."""
    return dict(moon,
                sin_node=sin_node * moon["cos_node"] - cos_node * moon["sin_node"],
                cos_node=moon["cos_node"] * cos_node + moon["sin_node"] * sin_node,
                amplitude=LUNAR_PERTURBATION)


def compute_lunisolar_geometry(days_since_1949, eccentricity, arg_perigee,
                               inclination, right_ascension, mean_motion):
    """Describe how the Sun and the Moon pull on the orbit (``dscom``)."""
    day = days_since_1949 + DAYS_FROM_1900

    sin_node, cos_node = np.sin(right_ascension), np.cos(right_ascension)
    emsq = eccentricity * eccentricity
    betasq = 1.0 - emsq
    orbit = dict(sin_inclination=np.sin(inclination), cos_inclination=np.cos(inclination),
                 sin_arg_perigee=np.sin(arg_perigee), cos_arg_perigee=np.cos(arg_perigee),
                 eccentricity=eccentricity, emsq=emsq, betasq=betasq,
                 rtemsq=np.sqrt(betasq), mean_motion=mean_motion)

    moon = _lunar_geometry(day)
    solar = _third_body_terms(_solar_geometry(sin_node, cos_node), orbit)
    lunar = _third_body_terms(_moon_relative_to_orbit(moon, sin_node, cos_node), orbit)

    geometry = dict(orbit)
    geometry["day"] = day
    geometry["solar"] = solar
    geometry["lunar"] = lunar
    geometry["zmol"] = np.fmod(4.7199672 + 0.22997150 * day - moon["mean_longitude"], TWO_PI)
    geometry["zmos"] = np.fmod(6.2565837 + 0.017201977 * day, TWO_PI)
    geometry.update(_periodic_coefficients(solar, lunar, emsq))
    return geometry


# Mean motion of the Sun and of the Moon about the Earth, radians per minute.
SOLAR_MEAN_MOTION = 1.19459e-5
LUNAR_MEAN_MOTION = 1.5835218e-4

# Within three degrees of the equator, prograde or retrograde, the node is too
# poorly defined to perturb.
POLAR_INCLINATION_LIMIT = 5.2359877e-2

# Below this inclination the node and the argument of perigee are replaced by
# the Lyddane variables, which stay well behaved as the orbit flattens out.
LYDDANE_INCLINATION_LIMIT = 0.2


def compute_secular_rates(geometry):
    """Give the drift the Sun and the Moon impose on the orbit (part of ``dsinit``).

    Returned as rates per minute for eccentricity, inclination, mean anomaly,
    argument of perigee and right ascension.
    """
    solar, lunar = geometry["solar"], geometry["lunar"]
    emsq = geometry["emsq"]
    sinim, cosim = geometry["sin_inclination"], geometry["cos_inclination"]

    ses = solar["s1"] * SOLAR_MEAN_MOTION * solar["s5"]
    sis = solar["s2"] * SOLAR_MEAN_MOTION * (solar["z11"] + solar["z13"])
    sls = -SOLAR_MEAN_MOTION * solar["s3"] * (solar["z1"] + solar["z3"] - 14.0 - 6.0 * emsq)
    sghs = solar["s4"] * SOLAR_MEAN_MOTION * (solar["z31"] + solar["z33"] - 6.0)
    shs = -SOLAR_MEAN_MOTION * solar["s2"] * (solar["z21"] + solar["z23"])

    inclination = np.arctan2(sinim, cosim)
    if _node_is_ill_defined(inclination):
        shs = 0.0
    if sinim != 0.0:
        shs = shs / sinim
    sgs = sghs - cosim * shs

    dedt = ses + lunar["s1"] * LUNAR_MEAN_MOTION * lunar["s5"]
    didt = sis + lunar["s2"] * LUNAR_MEAN_MOTION * (lunar["z11"] + lunar["z13"])
    dmdt = sls - LUNAR_MEAN_MOTION * lunar["s3"] * (lunar["z1"] + lunar["z3"] - 14.0 - 6.0 * emsq)
    sghl = lunar["s4"] * LUNAR_MEAN_MOTION * (lunar["z31"] + lunar["z33"] - 6.0)
    shll = -LUNAR_MEAN_MOTION * lunar["s2"] * (lunar["z21"] + lunar["z23"])
    if _node_is_ill_defined(inclination):
        shll = 0.0

    domdt = sgs + sghl
    dnodt = shs
    if sinim != 0.0:
        domdt = domdt - cosim / sinim * shll
        dnodt = dnodt + shll / sinim

    return dict(dedt=dedt, didt=didt, dmdt=dmdt, domdt=domdt, dnodt=dnodt)


def _node_is_ill_defined(inclination):
    """Tell whether the orbit lies too close to the equatorial plane."""
    return (inclination < POLAR_INCLINATION_LIMIT
            or inclination > np.pi - POLAR_INCLINATION_LIMIT)


def _wobble(coefficients, phase, eccentricity_of_orbit):
    """Evaluate one body's periodic contribution at the given phase."""
    angle = phase + 2.0 * eccentricity_of_orbit * np.sin(phase)
    sin_angle = np.sin(angle)
    f2 = 0.5 * sin_angle * sin_angle - 0.25
    f3 = -0.5 * sin_angle * np.cos(angle)
    two, three, four = coefficients
    return two * f2 + three * f3 + (0.0 if four is None else four * sin_angle)


def compute_periodics(geometry, minutes_since_epoch):
    """Give the periodic shift the Sun and the Moon impose (first half of ``dpper``)."""
    g = geometry
    solar_phase = g["zmos"] + SOLAR_MEAN_MOTION * minutes_since_epoch
    lunar_phase = g["zmol"] + LUNAR_MEAN_MOTION * minutes_since_epoch

    def solar(two, three, four=None):
        return _wobble((g[two], g[three], None if four is None else g[four]),
                       solar_phase, SOLAR_ECCENTRICITY)

    def lunar(two, three, four=None):
        return _wobble((g[two], g[three], None if four is None else g[four]),
                       lunar_phase, LUNAR_ECCENTRICITY)

    return dict(
        eccentricity=solar("se2", "se3") + lunar("ee2", "e3"),
        inclination=solar("si2", "si3") + lunar("xi2", "xi3"),
        mean_anomaly=solar("sl2", "sl3", "sl4") + lunar("xl2", "xl3", "xl4"),
        arg_perigee=solar("sgh2", "sgh3", "sgh4") + lunar("xgh2", "xgh3", "xgh4"),
        right_ascension=solar("sh2", "sh3") + lunar("xh2", "xh3"),
    )


def apply_periodics(shift, eccentricity, inclination, right_ascension, arg_perigee,
                    mean_anomaly):
    """Add the periodic shift to the orbital elements (second half of ``dpper``)."""
    inclination = inclination + shift["inclination"]
    eccentricity = eccentricity + shift["eccentricity"]
    sin_inc, cos_inc = np.sin(inclination), np.cos(inclination)

    if inclination >= LYDDANE_INCLINATION_LIMIT:
        node_shift = shift["right_ascension"] / sin_inc
        arg_perigee = arg_perigee + shift["arg_perigee"] - cos_inc * node_shift
        right_ascension = right_ascension + node_shift
        mean_anomaly = mean_anomaly + shift["mean_anomaly"]
        return eccentricity, inclination, right_ascension, arg_perigee, mean_anomaly

    return (eccentricity, inclination,
            *_apply_periodics_near_equator(shift, sin_inc, cos_inc, right_ascension,
                                           arg_perigee, mean_anomaly))


def _turn_to_positive(angle):
    """Bring an angle into a whole turn measured from zero.

    Here the right ascension is used as a number rather than as the argument of
    a sine or a cosine, so whether it is measured from zero or from minus half a
    turn changes the result. The operational software measures it from zero, and
    the element sets this model is fed come from that software.
    """
    return angle + TWO_PI if angle < 0.0 else angle


def _apply_periodics_near_equator(shift, sin_inc, cos_inc, right_ascension, arg_perigee,
                                  mean_anomaly):
    """Add the shift through the Lyddane variables, for a nearly equatorial orbit.

    Close to the equator the node swings wildly for a small change of the
    orbital plane, so the plane itself is shifted instead and the node read back
    from it.
    """
    sin_node, cos_node = np.sin(right_ascension), np.cos(right_ascension)
    plane_x = sin_inc * sin_node + (shift["right_ascension"] * cos_node
                                    + shift["inclination"] * cos_inc * sin_node)
    plane_y = sin_inc * cos_node + (-shift["right_ascension"] * sin_node
                                    + shift["inclination"] * cos_inc * cos_node)

    right_ascension = _turn_to_positive(np.fmod(right_ascension, TWO_PI))
    mean_longitude = (mean_anomaly + arg_perigee + shift["mean_anomaly"] + shift["arg_perigee"]
                      + (cos_inc - shift["inclination"] * sin_inc) * right_ascension)

    previous_node = right_ascension
    right_ascension = _turn_to_positive(np.arctan2(plane_x, plane_y))
    if np.abs(previous_node - right_ascension) > np.pi:
        right_ascension += TWO_PI if right_ascension < previous_node else -TWO_PI

    mean_anomaly = mean_anomaly + shift["mean_anomaly"]
    arg_perigee = mean_longitude - mean_anomaly - cos_inc * right_ascension
    return right_ascension, arg_perigee, mean_anomaly


def _periodic_coefficients(solar, lunar, emsq):
    """Combine the projections into the coefficients ``dpper`` applies."""
    return dict(
        se2=2.0 * solar["s1"] * solar["s6"],
        se3=2.0 * solar["s1"] * solar["s7"],
        si2=2.0 * solar["s2"] * solar["z12"],
        si3=2.0 * solar["s2"] * (solar["z13"] - solar["z11"]),
        sl2=-2.0 * solar["s3"] * solar["z2"],
        sl3=-2.0 * solar["s3"] * (solar["z3"] - solar["z1"]),
        sl4=-2.0 * solar["s3"] * (-21.0 - 9.0 * emsq) * SOLAR_ECCENTRICITY,
        sgh2=2.0 * solar["s4"] * solar["z32"],
        sgh3=2.0 * solar["s4"] * (solar["z33"] - solar["z31"]),
        sgh4=-18.0 * solar["s4"] * SOLAR_ECCENTRICITY,
        sh2=-2.0 * solar["s2"] * solar["z22"],
        sh3=-2.0 * solar["s2"] * (solar["z23"] - solar["z21"]),
        ee2=2.0 * lunar["s1"] * lunar["s6"],
        e3=2.0 * lunar["s1"] * lunar["s7"],
        xi2=2.0 * lunar["s2"] * lunar["z12"],
        xi3=2.0 * lunar["s2"] * (lunar["z13"] - lunar["z11"]),
        xl2=-2.0 * lunar["s3"] * lunar["z2"],
        xl3=-2.0 * lunar["s3"] * (lunar["z3"] - lunar["z1"]),
        xl4=-2.0 * lunar["s3"] * (-21.0 - 9.0 * emsq) * LUNAR_ECCENTRICITY,
        xgh2=2.0 * lunar["s4"] * lunar["z32"],
        xgh3=2.0 * lunar["s4"] * (lunar["z33"] - lunar["z31"]),
        xgh4=-18.0 * lunar["s4"] * LUNAR_ECCENTRICITY,
        xh2=-2.0 * lunar["s2"] * lunar["z22"],
        xh3=-2.0 * lunar["s2"] * (lunar["z23"] - lunar["z21"]),
    )
