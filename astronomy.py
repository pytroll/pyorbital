#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011 SMHI

# Author(s):

#   Martin Raspaud <martin.raspaud@smhi.se>

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
"""

import datetime
import numpy as np

def _jdays2000(current_time):
    """Get the days since 2000.
    """
    return _days(current_time - datetime.datetime(2000, 1, 1, 12, 0))
    

def _jdays(current_time):
    """Get the julian day of *current_time*.
    """
    return _jdays2000(current_time) + 2451545

def _days(d_t):
    """Get the days (floating point) from *d_t*.
    """
    return (d_t.days +
            (d_t.seconds +
             d_t.microseconds / (1000000.0)) / (24 * 3600.0))

def gmst(current_time):
    """Greenwich mean sidereal current_time, in radians.
    http://celestrak.com/columns/v02n02/
    """
    now = current_time
    #now = datetime.datetime(1995, 10, 1, 9, 0)
    now0 = datetime.datetime(now.year, now.month, now.day)
    epoch = datetime.datetime(2000, 1, 1, 12, 0)
    du2 = _days(now - epoch)
    d_u = _days(now0 - epoch)

    dus = (du2 - d_u) * 86400
    t_u = d_u / 36525.0
    theta_g_0 = (24110.54841 + t_u * (8640184.812866 +
                                      t_u * (0.093104 - t_u * 6.2 * 10e-6)))
    theta_g = (theta_g_0 + dus * 1.00273790934) % 86400
    return (theta_g / 86400.0) * 2 * np.pi

def lmst(current_time, longitude):
    return gmst(current_time) + longitude


def hour_angle(current_time, longitude, right_ascension):
    return ((lmst(current_time, longitude) * np.pi / 12.0 - right_ascension)
            % (2 * np.pi))
    

def sun_mean_longitude(current_time):
    jd = _jdays(current_time)
    t=jd/36525.0 + 1

    g2r = np.pi / 180

    al=279.696680+36000.768920*t+0.00030250*t*t
    am=358.47583+35999.04975*t-0.000150*t*t+0.0000033*t*t*t
    print al, am
    al = (al % 360.0) * g2r
    am = (am % 360.0) * g2r
    print np.rad2deg((al, am))


    c1=(1.919460-0.004789*t-0.000014*t*t)
    c2=(0.020094-0.000100*t)
    ec=(c1*np.sin(am)+c2*np.sin(2.0*am)+0.000293*np.sin(3.0*am))*g2r

    alon=al+ec

    # Perturbations

    a=(153.23+22518.7541*t)
    b=(216.57+45037.5082*t)
    c=(312.69+32964.3577*t)
    d=(350.74+445267.1142*t-0.00144*t*t)
    e=(231.19+20.20*t)

    a=(a%360.0)*g2r
    b=(b%360.0)*g2r
    c=(c%360.0)*g2r
    d=(d%360.0)*g2r
    e=(e%360.0)*g2r

    alon=alon+(  0.00134*np.cos(a)
                 + 0.00154*np.cos(b)
                 + 0.00200*np.cos(c)
                 + 0.00179*np.sin(d)
                 + 0.00178*np.sin(e)) * g2r

    om=(259.18-1934.142*t)
    om=(om % 360.0)*g2r
    if(om<0):
        om += 2.0*np.pi
    alon=alon-(0.00569+0.00479*np.sin(om))*g2r

    alon = (alon % (2.0*np.pi))
    
    return(alon)

def ecliptic_obliquity(jd):
    u = (jd)/3652500.0
    
    a1 = 2.18-3375.70*u+0.36*u*u
    a2 = 3.51+125666.39*u+0.10*u*u
    e  = (0.4090928
          -     0.0226938*u
          -    75.0E-7*u*u
          + 96926.0E-7*u*u*u
          -  2491.0E-7*u*u*u*u
          - 12104.0E-7*u*u*u*u*u
          + (446.0*np.cos(a1)+28.0*np.cos(a2))*1.0e-7)
    return(e)

def sml(current_time):
    jd = _jdays2000(current_time)
    w = 282.9404 + 4.70935e-5 * jd
    print "longitude of perihelion", w
    a = 1
    e = 0.016709 - 1.151e-9 * jd
    print "eccentricity", e
    M = (356.047 + 0.9856002585 * jd) % 360
    print "mean anomaly", M
    oblecl = 23.4393 - 3.563e-7 * jd
    L = (w + M)%360
    print L

def ecliptic2equatorial(current_time, lon, lat):
    eps = ecliptic_obliquity(_jdays2000(current_time))
    ra = np.arctan2(-np.sin(lat) * np.sin(eps) +
                    np.cos(lat) * np.cos(eps) * np.sin(lon),
                    np.cos(lon) * np.cos(lat)) % (2 * np.pi)
    dec = np.arcsin(np.sin(lat) * np.cos(eps) +
                    np.cos(lat) * np.sin(eps) * np.sin(lon))
    return ra, dec

def _cos_sun_zenith_angles(current_time, lons, lats):
    ra, dec = ecliptic2equatorial(current_time,
                                  sun_mean_longitude(current_time),
                                  0)
    #Alt=asin(sinlat*sindelta+coslat*cosdelta*cosh);
    #Az = atan2(cosdelta*sinh,-sindelta*coslat+cosdelta*sinlat*cosh);
    return (np.sin(lats) * np.sin(dec) +
            np.cos(lats) * np.cos(dec) *
            np.cos(hour_angle(current_time, lons, ra)))

def sun_zenith_angles(current_time, lons, lats):
    return np.arccos(_cos_sun_zenith_angles(current_time, lons, lats))


# FROM http://www.geoastro.de/elevaz/basics/index.htm


def sun_ecliptic_longitude(current_time):
    T = _jdays2000(current_time) / 36525.0
    # mean anomaly, rad
    M = np.deg2rad(357.52910 + 35999.05030*T - 0.0001559*T*T - 0.00000048*T*T*T)
    print np.rad2deg(M)
    # mean longitude, deg
    L0 = 280.46645 + 36000.76983*T + 0.0003032*T*T
    DL = (1.914600 - 0.004817*T - 0.000014*T*T)*np.sin(M) + (0.019993 - 0.000101*T)*np.sin(2*M) + 0.000290*np.sin(3*M)
    # true longitude, deg
    L = L0 + DL
    return np.deg2rad(L)

def cmp_right_ascension(current_time):
    T = _jdays2000(current_time) / 36525.0
    eps = np.deg2rad(23.0 + 26.0/60.0 + 21.448/3600.0 - (46.8150*T + 0.00059*T*T - 0.001813*T*T*T)/3600)
    L = sun_ecliptic_longitude(current_time)
    X = np.cos(L)
    Y = np.cos(eps)*np.sin(L)
    Z = np.sin(eps)*np.sin(L)
    R = np.sqrt(1.0-Z*Z)
    # sun declination
    delta = np.arctan2(Z, R)
    # right ascension
    RA = 2 * np.arctan2(Y, (X+R))
    return RA, delta

def local_hour_angle(current_time, lon, ra):
    return lmst(current_time, lon) - ra

def get_alt_az(current_time, lon, lat):
    ra, dec = cmp_right_ascension(current_time)
    h = local_hour_angle(current_time, lon, ra)
    return np.arcsin(np.sin(lat)*np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h)), np.arctan2(-np.sin(h), (np.cos(lat)*np.tan(dec) - np.sin(lat)*np.cos(h)))

def cos_zen(current_time, lon, lat):
    ra, dec = cmp_right_ascension(current_time)
    h = local_hour_angle(current_time, lon, ra)
    return (np.sin(lat)*np.sin(dec) + np.cos(lat) * np.cos(dec) * np.cos(h))
    
if __name__ == '__main__':
    t = datetime.datetime(2009, 8, 10, 14, 30)
#    t = datetime.datetime(1990, 4, 19, 0, 0)
    lat = np.deg2rad(58)
    lon = np.deg2rad(16)
    #sml(t)
    
    #print "solar zenith angle", np.rad2deg(sun_zenith_angles(t, lon, lat))
    
    print "ra, declination", ecliptic2equatorial(t, lon, lat)

    import ephem
    o = ephem.Observer()
    o.lat, o.long, o.date = lat, lon, t
    sun = ephem.Sun(o)
    #print "jd", _jdays(t)

    #t = datetime.datetime(2011, 1, 1, 14, 30)
    #ra, dec = np.rad2deg(ecliptic2equatorial(t, lon, lat))
    #sun_mean_longitude(t)

    d = _jdays2000(t)

    i = 0.0
    w = np.deg2rad(282.9404 + 4.70935E-5 * d)
    a = 1.000000 
    e = np.deg2rad(0.016709 - 1.151E-9 * d)
    M = np.deg2rad(356.0470 + 0.9856002585 * d)
