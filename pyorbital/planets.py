#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012

# Author(s):

#   Adam Dybbroe <adam.dybbroe@smhi.se>

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Astronomy module.
Parts taken from http://www.stjarnhimlen.se/comp/ppcomp.html
Positions of the sun, the moon, the earth and the other planets
"""

import datetime
import numpy as np


class Moon(object):
    def __init__(self, dtobj):
        # Orbital elements:
        self.orbelem = {"N": 0.0, "i": 0.0, "w": 0.0,
                        "a": 0.0, "e": 0.0, "M": 0.0,
                        "E": 0.0, 'v': 0.0}
        self.datetime = dtobj
        self.day_number = day_number(dtobj)
        #self.lonsun = 0.0
        self.hour_angle = 0.0
        self.right_ascension = 0.0
        self.declination = 0.0

        # Pertubations:
        self.dlon = 0.0
        self.dlat = 0.0
        self.drad = 0.0
        self.pertubations()

        self.get_lonsun()
        self.rdist = 0.0
        self.ecl = obliquity_of_ecliptic(self.day_number)

        # Get the orbital elements of the Sun
        orbelem = orbital_elem(self.day_number, body='moon')

        orbelem = eccentric_anomaly(orbelem)
        self.rdist, orbelem = distance_and_true_anomaly(orbelem)

        self.orbelem = orbelem
        pos = position(orbelem, self.rdist)
        self.geocentric_position = pos
        # Apply the pertubations:
        self.rdist = self.rdist + self.drad

        xg_ = self.geocentric_position['xh']
        yg_ = self.geocentric_position['yh']
        zg_ = self.geocentric_position['zh']

        # Cartesian to polar coordinates:
        lon = np.rad2deg(np.arctan2(yg_, xg_)) + self.dlon
        lat = np.rad2deg(np.arctan2(zg_, np.sqrt(
            xg_ * xg_ + yg_ * yg_))) + self.dlat

        # Polar to cartesian coordinates:
        xg_ = self.rdist * np.cos(np.deg2rad(lon))
        yg_ = self.rdist * np.sin(np.deg2rad(lon))
        zg_ = self.rdist * np.sin(np.deg2rad(lat))
        self.geocentric_position['xh'] = xg_
        self.geocentric_position['yh'] = yg_
        self.geocentric_position['zh'] = zg_

        # The Moon's parallax, i.e. the apparent size of the (equatorial)
        # radius of the Earth, as seen from the Moon:
        self.parallax = np.rad2deg(np.arcsin(1. / self.rdist))

        right_ascension, declination, gdist = ecliptic_to_equatorial(xg_,
                                                                     yg_,
                                                                     zg_,
                                                                     self.ecl)
        self.right_ascension = right_ascension
        self.declination = declination

    def topocentric_position(self, longitude, latitude):
        """The Moon's position, as computed earlier, is geocentric, i.e. as
        seen by an imaginary observer at the center of the Earth. Real
        observers dwell on the surface of the Earth, though, and they will see
        a different position - the topocentric position. This position can
        differ by more than one degree from the geocentric position. To compute
        the topocentric positions, we must add a correction to the geocentric
        position.
        """

        azi, alt_geoc, ha_ = azimuthal_coord(self.datetime,
                                             self.lonsun,
                                             longitude,
                                             latitude,
                                             self.right_ascension,
                                             self.declination)
        #print "azi, alt_geoc, ha_: ", azi, alt_geoc, ha_

        self.hour_angle = ha_
        alt_topoc = alt_geoc - self.parallax * np.cos(np.deg2rad(alt_geoc))

        # Account for the flattening of the Earth:
        gclat = latitude - 0.1924 * np.sin(np.deg2rad(2*latitude))
        rho = 0.99833 + 0.00167 * np.cos(np.deg2rad(2*latitude))

        g__ = np.arctan(np.tan(np.deg2rad(gclat)) /
                        np.cos(np.deg2rad(self.hour_angle)))
        g__ = np.rad2deg(g__)

        # Now we're ready to convert the geocentric Right Ascention and
        # Declination (RA, Decl) to their topocentric values (topRA, topDecl):

        topRA = (self.right_ascension -
                 self.parallax * rho * np.cos(np.deg2rad(gclat)) *
                 np.sin(np.deg2rad(ha_)) /
                 np.cos(np.deg2rad(self.declination))
                 )
        topDecl = (self.declination -
                   self.parallax * rho * np.sin(np.deg2rad(gclat)) *
                   np.sin(np.deg2rad(g__ - self.declination)) /
                   np.sin(np.deg2rad(g__))
                   )

        return topRA, topDecl, alt_topoc, azi

    def get_lonsun(self):
        """Calculate the lonsun"""
        # Ls = Ms + ws       Mean Longitude of the Sun  (Ns=0)
        sun = Sun(self.datetime)
        self.lonsun = sun.lonsun

    def pertubations(self):
        """Pertubations of the moon, to be added to 
        its lon,lat and distance"""

        # Ms, Mm             Mean Anomaly of the Sun and the Moon
        # Nm                 Longitude of the Moon's node
        # ws, wm             Argument of perihelion for the Sun and the Moon
        # Ls = Ms + ws       Mean Longitude of the Sun  (Ns=0)
        # Lm = Mm + wm + Nm  Mean longitude of the Moon
        # D = Lm - Ls        Mean elongation of the Moon
        # F = Lm - Nm        Argument of latitude for the Moon

        sun = Sun(self.datetime)
        Ms_ = sun.orbelem['M']
        Mm_ = self.orbelem['M']
        Nm_ = self.orbelem['N']
        wm_ = self.orbelem['w']
        Lm_ = Mm_ + wm_ + Nm_
        D__ = Lm_ - sun.lonsun
        F__ = Lm_ - Nm_

        # Add these terms to the Moon's longitude (degrees):
        dlon = (-1.274 * np.sin(np.deg2rad(Mm_ - 2*D__))    # (the Evection)
                + 0.658 * np.sin(np.deg2rad(2*D__))         # (the Variation)
                # (the Yearly Equation)
                - 0.186 * np.sin(np.deg2rad(Ms_))
                - 0.059 * np.sin(np.deg2rad(2*Mm_ - 2*D__))
                - 0.057 * np.sin(np.deg2rad(Mm_ - 2*D__ + Ms_))
                + 0.053 * np.sin(np.deg2rad(Mm_ + 2*D__))
                + 0.046 * np.sin(np.deg2rad(2*D__ - Ms_))
                + 0.041 * np.sin(np.deg2rad(Mm_ - Ms_))
                # (the Parallactic Equation)
                - 0.035 * np.sin(np.deg2rad(D__))
                - 0.031 * np.sin(np.deg2rad(Mm_ + Ms_))
                - 0.015 * np.sin(np.deg2rad(2*F__ - 2*D__))
                + 0.011 * np.sin(np.deg2rad(Mm_ - 4*D__))
                )
        self.dlon = dlon

        # Add these terms to the Moon's latitude (degrees):
        dlat = (-0.173 * np.sin(np.deg2rad(F__ - 2*D__))
                - 0.055 * np.sin(np.deg2rad(Mm_ - F__ - 2*D__))
                - 0.046 * np.sin(np.deg2rad(Mm_ + F__ - 2*D__))
                + 0.033 * np.sin(np.deg2rad(F__ + 2*D__))
                + 0.017 * np.sin(np.deg2rad(2*Mm_ + F__))
                )
        self.dlat = dlat

        # Add these terms to the Moon's distance (Earth radii):
        drad = (-0.58 * np.cos(np.deg2rad(Mm_ - 2*D__))
                - 0.46 * np.cos(np.deg2rad(2*D__))
                )
        self.drad = drad


class Sun(object):
    def __init__(self, dtobj):
        # Orbital elements:
        self.orbelem = {"N": 0.0, "i": 0.0, "w": 0.0,
                        "a": 1.0, "e": 0.0, "M": 0.0,
                        "E": 0.0, 'v': 0.0}
        self.datetime = dtobj
        self.day_number = day_number(dtobj)
        self.lonsun = 0.0
        self.rdist = 0.0
        self.ecl = obliquity_of_ecliptic(self.day_number)

        # Get the orbital elements of the Sun
        orbelem = orbital_elem(self.day_number, body='sun')

        orbelem = sun_eccentric_anomaly(orbelem)
        self.rdist, orbelem = distance_and_true_anomaly(orbelem)

        self.orbelem = orbelem

        # True longitude:
        self.lonsun = self.orbelem['v'] + self.orbelem['w']

        self.position = {}

    def get_position(self):
        """Calculate the position of the Sun. 
        Equatorial, rectangular, geocentric coordinates"""

        # Convert lonsun,r to ecliptic rectangular geocentric coordinates xs,ys:
        xs_ = self.rdist * np.cos(np.deg2rad(self.lonsun))
        ys_ = self.rdist * np.sin(np.deg2rad(self.lonsun))
        # To convert this to equatorial, rectangular, geocentric coordinates, compute:
        xe_ = xs_
        ye_ = ys_ * np.cos(np.deg2rad(self.ecl))
        ze_ = ys_ * np.sin(np.deg2rad(self.ecl))

        # Compute the planet's Right Ascension (RA) and Declination (Dec):
        ra_ = np.rad2deg(np.arctan2(ye_, xe_))
        dec = np.rad2deg(np.arctan2(ze_, np.sqrt(xe_ * xe_ + ye_ * ye_)))

        self.position = {'xe': xe_, 'ye': ye_, 'ze': ze_,
                         'right_ascension': ra_, 'declination': dec}


def day_number(dtobj):
    """Get the 'day number' from the date, given by a datetime object
    """
    # d = 367*y - 7 * ( y + (m+9)/12 ) / 4 + 275*m/9 + D - 730530
    # dnumber = (367 * dtobj.year -
    #           7 * ( dtobj.year + (dtobj.month + 9) / 12 ) / 4 +
    #           275 * dtobj.month / 9 + dtobj.day - 730530)
    # return dnumber + dtobj.hour/24.0

    origo = datetime.datetime(1999, 12, 31, 0)
    delta_t = dtobj - origo
    return delta_t.days + (delta_t.seconds +
                           delta_t.microseconds/1000000.)/(3600. * 24.)


def _gmst(lonsun, utc):
    """The Greenwich Mean Sideral Time (GMST) is the LST at Greenwich"""
    # GMST0 = 15 * (Ls + 180_degrees)
    # GMST = GMST0 + UT
    # LST  = GMST + local_longitude/15
    gmst0 = ((lonsun + 180.0) % 360) / 15.0
    ut_ = (utc.hour + utc.minute / 60.0 +
           (utc.second + utc.microsecond / 1000000.) / 3600.)
    return gmst0 + ut_


def _lst(gmst, local_lon):
    """The Local Sideral Time (LST).
    This is simply the RA of your local meridian"""
    return gmst + local_lon / 15.0


def _hour_angle(utc_time, local_lon, lonsun, right_ascension):
    """Hour angle at *utc_time* for the given *local_lon* and
    *right_ascension*. gmst and lst in hours. the hour angle HA is given in
    degrees:
    """
    gmst = _gmst(lonsun, utc_time)
    return (15 * _lst(gmst, local_lon) - right_ascension) % 360


def azimuthal_coord(utc_time, lonsun, local_lon, local_lat,
                    right_ascension, declination):
    """Calculate the azimuthal coordinates (azimuth and altitude)"""

    ha_ = _hour_angle(utc_time, local_lon, lonsun, right_ascension)

    x__ = np.cos(np.deg2rad(ha_)) * np.cos(np.deg2rad(declination))
    y__ = np.sin(np.deg2rad(ha_)) * np.cos(np.deg2rad(declination))
    z__ = np.sin(np.deg2rad(declination))

    xhor = (x__ * np.sin(np.deg2rad(local_lat)) -
            z__ * np.cos(np.deg2rad(local_lat)))
    yhor = y__
    zhor = (x__ * np.cos(np.deg2rad(local_lat)) +
            z__ * np.sin(np.deg2rad(local_lat)))

    azi = np.rad2deg(np.arctan2(yhor, xhor)) + 180.0
    # = atan2( zhor, sqrt(xhor*xhor+yhor*yhor) )
    alt = np.rad2deg(np.arcsin(zhor))

    return azi, alt, ha_


def obliquity_of_ecliptic(dnum):
    """compute the obliquity of the ecliptic"""
    return 23.4393 - 3.563E-7 * dnum


def orbital_elem(dnum, **options):
    """Orbital elements of the Sun, the Moon, and the planets"""

    if "body" in options and options["body"] == 'moon':
        N__ = 125.1228 - 0.0529538083 * dnum
        i__ = 5.1454
        w__ = 318.0634 + 0.1643573223 * dnum
        a__ = 60.2666  # (Earth radii)
        e__ = 0.054900
        M__ = 115.3654 + 13.0649929509 * dnum
    elif "body" in options and options["body"] == 'sun':
        N__ = 0.0
        i__ = 0.0
        w__ = 282.9404 + 4.70935E-5 * dnum
        a__ = 1.000000  # (AU)
        e__ = 0.016709 - 1.151E-9 * dnum
        M__ = 356.0470 + 0.9856002585 * dnum
    else:
        raise IOError("Celestial body not specified or not supported yet!")

    if M__ > 360:
        M__ = M__ - 360 * (int(M__)/360)
    if N__ > 360:
        N__ = N__ - 360 * (int(N__)/360)

    return {"N": N__, "i": i__, "w": w__,
            "a": a__, "e": e__, "M": M__}


def sun_eccentric_anomaly(orbelem):
    """Computation of the eccentric anomaly from the mean anomaly
    for ths Sun"""
    # Note that the formulae for computing E are not exact;
    # however they're accurate enough here.
    E__ = (orbelem["M"] +
           np.rad2deg(orbelem["e"]) *
           np.sin(np.deg2rad(orbelem["M"])) *
           (1.0 + orbelem["e"] * np.cos(np.deg2rad(orbelem["M"]))))

    orbelem["E"] = E__
    return orbelem


def eccentric_anomaly(orbelem):
    """Computation of the eccentric anomaly from the mean anomaly"""

    E__ = (orbelem["M"] +
           np.rad2deg(orbelem["e"]) *
           np.sin(np.deg2rad(orbelem["M"])) *
           (1.0 + orbelem["e"] * np.cos(np.deg2rad(orbelem["M"]))))

    E0_ = E__
    epsilon = 0.001
    delta_e = 2 * epsilon
    niter = 0
    while delta_e > epsilon and niter < 10:
        #print "E0_ = ",E0_
        E1_ = (E0_ -
               (E0_ - np.rad2deg(orbelem["e"]) * np.sin(np.deg2rad(E0_))
                - orbelem["M"]) /
               (1 - orbelem["e"] * np.cos(np.deg2rad(E0_))))
        delta_e = abs(E1_ - E0_)
        #print "Iteration = %d, dev = %f" % (niter, delta_e)
        niter = niter + 1
        E0_ = E1_

    orbelem["E"] = E1_
    return orbelem


def distance_and_true_anomaly(orbelem):
    """Computation of the planets distance 'r' and its true anomaly 'v'"""
    # xv = r * cos(v)
    # yv = r * sin(v)
    xv_ = orbelem["a"] * (np.cos(np.deg2rad(orbelem["E"])) - orbelem["e"])
    yv_ = orbelem["a"] * (np.sqrt(1.0 - orbelem["e"] * orbelem["e"]) *
                          np.sin(np.deg2rad(orbelem["E"])))

    v__ = np.rad2deg(np.arctan2(yv_, xv_))
    r__ = np.sqrt(xv_ * xv_ + yv_ * yv_)
    orbelem["v"] = v__

    return r__, orbelem


def position(orbelem, radius):
    """Computation of the planet's position in 3-dimensional space. For the
    Moon, this is the geocentric (Earth-centered) position in the ecliptic
    coordinate system. For the planets, this is the heliocentric (Sun-centered)
    position, also in the ecliptic coordinate system.
    """

    xh_ = radius * (np.cos(np.deg2rad(orbelem["N"])) *
                    np.cos(np.deg2rad(orbelem["v"] + orbelem["w"])) -
                    np.sin(np.deg2rad(orbelem["N"])) *
                    np.sin(np.deg2rad(orbelem["v"] + orbelem["w"])) *
                    np.cos(np.deg2rad(orbelem["i"])))
    yh_ = radius * (np.sin(np.deg2rad(orbelem["N"])) *
                    np.cos(np.deg2rad(orbelem["v"] + orbelem["w"])) +
                    np.cos(np.deg2rad(orbelem["N"])) *
                    np.sin(np.deg2rad(orbelem["v"] + orbelem["w"])) *
                    np.cos(np.deg2rad(orbelem["i"])))
    zh_ = radius * (np.sin(np.deg2rad(orbelem["v"] + orbelem["w"])) *
                    np.sin(np.deg2rad(orbelem["i"])))

    # Compute the ecliptic longitude and latitude
    lonecl = np.rad2deg(np.arctan2(yh_, xh_))
    latecl = np.rad2deg(np.arctan2(zh_, np.sqrt(xh_ * xh_ + yh_ * yh_)))
    check = np.sqrt(xh_ * xh_ + yh_ * yh_ + zh_ * zh_)

    #print "sqrt(xh*xh+yh*yh+zh*zh) = r: ", check, radius
    return {'xh': xh_, 'yh': yh_, 'zh': zh_,
            'lonecl': lonecl, 'latecl': latecl}


def heliocentric_to_geocentric(lonecl, latecl, radius, lonsun, rsun):
    """Conversion from heliocentric to geocentric coordinates. (Not needed for
    the moon). 
    Now we have computed the heliocentric (Sun-centered) coordinate
    of the planet, and we have included the most important perturbations. We
    want to compute the geocentric (Earth-centerd) position. We should convert
    the perturbed lonecl, latecl, r to (perturbed) xh, yh, zh
    """
    xh_ = radius * np.cos(lonecl) * np.cos(np.deg2rad(latecl))
    yh_ = radius * np.sin(lonecl) * np.cos(np.deg2rad(latecl))
    zh_ = radius * np.sin(np.deg2rad(latecl))

    # If we are computing the Moon's position, this is already the geocentric
    # position, and thus we simply set xg=xh, yg=yh, zg=zh. Otherwise we must
    # also compute the Sun's position: convert lonsun, rs (where rs is the r
    # computed here) to xs, ys:

    xs_ = rsun * np.cos(np.deg2rad(lonsun))
    ys_ = rsun * np.sin(np.deg2rad(lonsun))

    xg_ = xh_ + xs_
    yg_ = yh_ + ys_
    zg_ = zh_

    return {'xg': xg_, 'yg': yg_, 'zg': zg_}


def ecliptic_to_equatorial(xg_, yg_, zg_, ecl):
    """Convert the rectangular, ecliptic coordinates to rectangular, 
    equatorial coordinates: Simply rotate the y-z-plane by ecl, the 
    angle of the obliquity of the ecliptic"""

    xe_ = xg_
    ye_ = (yg_ * np.cos(np.deg2rad(ecl)) -
           zg_ * np.sin(np.deg2rad(ecl)))
    ze_ = (yg_ * np.sin(np.deg2rad(ecl)) +
           zg_ * np.cos(np.deg2rad(ecl)))

    # Compute the planet's Right Ascension (RA) and Declination (Dec):
    ra_ = np.rad2deg(np.arctan2(ye_, xe_))
    dec = np.rad2deg(np.arctan2(ze_, np.sqrt(xe_ * xe_ + ye_ * ye_)))

    # Compute the geocentric distance:
    # = sqrt(xe*xe+ye*ye+ze*ze)
    rg_ = np.sqrt(xg_ * xg_ + yg_ * yg_ + zg_ * zg_)

    return ra_, dec, rg_


if __name__ == "__main__":
    now = datetime.datetime.utcnow()

    #dNum = day_number(now)
    #orb = orbital_elem(dNum, body='moon')
    #orb = eccentric_anomaly(orb)
    #rdist, orb = distance_and_true_anomaly(orb)

    #pos = position(orb, rdist)
    #xgeo, ygeo, zgeo = pos['xh'], pos['yh'], pos['zh']

    #ecliptic = obliquity_of_ecliptic(dNum)

    #tup = ecliptic_to_equatorial(xgeo, ygeo, zgeo, ecliptic)

    # Norrköping: Lat N 58° 35′ 15″ Lon E 16° 11′ 15″
    llon = 16 + 11.0/60. + 15./3600.
    llat = 58. + 35./60. + 15./3600.
    #moon = Moon(now)
    #rasc, decl, alt, azi = moon.topocentric_position(17., 58.0)
    currtime = datetime.datetime(2012, 1, 7, 12, 0)
    endtime = datetime.datetime(2012, 1, 7, 13, 0)
    delta_t = datetime.timedelta(seconds=120)
    heights = []
    dtimes = []
    while currtime < endtime:
        moon = Moon(currtime)
        rasc, decl, alt, azi = moon.topocentric_position(llon, llat)
        heights.append(alt)
        dtimes.append(currtime)
        currtime = currtime + delta_t
        break

    llon = np.arange(100)/10.0 + 10.0
    llat = np.arange(100)/10.0 + 50.0
    rasc, decl, alt, azi = moon.topocentric_position(llon, llat)
