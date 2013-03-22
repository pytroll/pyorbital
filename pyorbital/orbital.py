#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2011, 2012, 2013.

# Author(s):

#   Esben S. Nielsen <esn@dmi.dk>
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

"""Module for computing the orbital parameters of satellites.
"""

from datetime import datetime, timedelta
import numpy as np
from pyorbital import tlefile
from pyorbital import astronomy


ECC_EPS = 1.0e-6	# Too low for computing further drops.
ECC_LIMIT_LOW = -1.0e-3
ECC_LIMIT_HIGH = 1.0 - ECC_EPS	# Too close to 1 
ECC_ALL = 1.0e-4

EPS_COS = 1.5e-12
EPS_SIN = 1.0e-12
EPS_COSIO = 1.5e-12

NR_EPS = 1.0e-12

CK2 = 5.413080e-4
CK4 = 0.62098875e-6
E6A = 1.0e-6
QOMS2T = 1.88027916e-9
S = 1.01222928
S0 = 78.0
XJ3 = -0.253881e-5
#XKE = 0.743669161e-1
XKE = 7.43669161331734132e-2
XKMPER = 6378.135
XMNPDA = 1440.0
#MFACTOR = 7.292115E-5
AE = 1.0
SECDAY = 8.6400E4

F = 1 / 298.257223563 # Earth flattening WGS-84
A = 6378.137 # WGS84 Equatorial radius    


SGDP4_ZERO_ECC = 0
SGDP4_NEAR_SIMP = 1 
SGDP4_NEAR_NORM = 2
SGDP4_DEEP_NORM = 3
SGDP4_DEEP_RESN = 4
SGDP4_DEEP_SYNC = 5

KS = AE * (1.0 + S0 / XKMPER)
A3OVK2 = (-XJ3 / CK2) * AE ** 3

ZNS = 1.19459e-5
C1SS = 2.9864797e-6
ZES = 0.01675

ZNL = 1.5835218e-4
C1L = 4.7968065e-7
ZEL = 0.0549

ZCOSIS = 0.91744867
ZSINIS = 0.39785416
ZCOSGS = 0.1945905
ZSINGS = -0.98088458

ROOT22 = 1.7891679e-6
ROOT32 = 3.7393792e-7
ROOT44 = 7.3636953e-9
ROOT52 = 1.1428639e-7
ROOT54 = 2.1765803e-9

G22 = 5.7686396
G32 = 0.95240898
G44 = 1.8014998
G52 = 1.0508330
G54 = 4.4108898

THDT = 4.37526908801129966e-3

STEP = 720.0


class OrbitalError(Exception):
    pass


class Orbital(object):
    """Class for orbital computations.

    The *satellite* parameter is the name of the satellite to work on and is
    used to retreive the right TLE data for internet or from *tle_file* in case
    it is provided.
    """

    def __init__(self, satellite, tle_file=None, line1=None, line2=None):
        satellite = satellite.upper()
        self.satellite_name = satellite
        self.tle = tlefile.read(satellite, tle_file=tle_file, 
                                line1=line1, line2=line2)
        self.orbit_elements = OrbitElements(self.tle)
        self._sgdp4 = _SGDP4(self.orbit_elements)
        
        
        pos_epoch, vel_epoch = self.get_position(self.tle.epoch, 
                                                 normalize=False)
        
        """
        if np.abs(pos_epoch[2]) > 1 or not vel_epoch[2] > 0:
            # Epoch not at ascending node
            self.orbit_elements.an_time = self.get_last_an_time(self.tle.epoch)
        else:
            # Epoch at ascending node (z < 1 km) and positive v_z
            self.orbit_elements.an_time = self.tle.epoch
            
        self.orbit_elements.an_period = self.orbit_elements.an_time - \
                        self.get_last_an_time(self.orbit_elements.an_time)"""
            

    def __str__(self):
        return self.satellite_name + " " + str(self.tle)

    def get_last_an_time(self, utc_time):
        """Calculate time of last ascending node relative to the
        specified time
        """        
        
        # Propagate backwards to ascending node
        dt = timedelta(minutes=10)
        t_old = utc_time
        t_new = t_old - dt
        pos0, vel0 = self.get_position(t_old, normalize=False)
        pos1, vel1 = self.get_position(t_new, normalize=False)
        while not (pos0[2] > 0 and pos1[2] < 0):
            pos0, vel0 = pos1, vel1
            t_old = t_new
            t_new = t_old - dt
            pos1, vel1 = self.get_position(t_new, normalize=False)
        
        # Return if z within 1 km of an
        if np.abs(pos0[2]) < 1:
            return t_old
        elif np.abs(pos1[2]) < 1:
            return t_new
    
        # Bisect to z within 1 km
        while np.abs(pos1[2]) > 1:
            pos0, vel0 = pos1, vel1
            dt = (t_old - t_new) / 2
            t_mid = t_old - dt
            pos1, vel1 = self.get_position(t_mid, normalize=False)
            if pos1[2] > 0:
                t_old = t_mid
            else:
                t_new = t_mid
        
        return t_mid                                

    def get_position(self, utc_time, normalize=True):
        """Get the cartesian position and velocity from the satellite.
        """
        
        kep = self._sgdp4.propagate(utc_time)
        pos, vel = kep2xyz(kep)
        
        if normalize:
            pos /= XKMPER
            vel /= XKMPER * XMNPDA / SECDAY

        return pos, vel


    def get_lonlatalt(self, utc_time):
        """Calculate sublon, sublat and altitude of satellite.
        http://celestrak.com/columns/v02n03/
        """
        (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = self.get_position(utc_time, normalize=True)

        lon = ((np.arctan2(pos_y * XKMPER, pos_x * XKMPER) - astronomy.gmst(utc_time))
               % (2 * np.pi))

        lon = np.where(lon > np.pi, lon - np.pi * 2, lon)
        lon = np.where(lon <= -np.pi, lon + np.pi *2, lon)

        r = np.sqrt(pos_x ** 2 + pos_y ** 2)
        lat = np.arctan2(pos_z, r)
        e2 = F * (2 - F)
        while True:
            lat2 = lat
            c = 1/(np.sqrt(1 - e2 * (np.sin(lat2) ** 2)))
            lat = np.arctan2(pos_z + c * e2 *np.sin(lat2), r)
            if np.all(abs(lat - lat2) < 1e-10):
                break
        alt = r / np.cos(lat)- c;
        alt *= A
        return np.rad2deg(lon), np.rad2deg(lat), alt

    def find_aos(self, time, lon, lat):
        pass

    def find_aol(self, time, lon, lat):
        pass

    def get_observer_look(self, time, lon, lat, alt):
        """Calculate observers look angle to a satellite.
        http://celestrak.com/columns/v02n02/

        time: Observation time (datetime object)
        lon: Longitude of observer position on ground
        lat: Latitude of observer position on ground
        alt: Altitude above sea-level (geoid) of observer position on ground

        Return: (Azimuth, Elevation)
        """
        (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = self.get_position(time, normalize=False)
        (opos_x, opos_y, opos_z), (ovel_x, ovel_y, ovel_z) = \
                                    astronomy.observer_position(time, lon, lat, alt)
        
        lon = np.deg2rad(lon)
        lat = np.deg2rad(lat)
        
        theta = (astronomy.gmst(time) + lon) % (2 * np.pi)

        rx = pos_x - opos_x
        ry = pos_y - opos_y
        rz = pos_z - opos_z
        rvx = vel_x - ovel_x
        rvy = vel_y - ovel_y
        rvz = vel_z - ovel_z

        sin_lat = np.sin(lat)
        cos_lat = np.cos(lat)
        sin_theta = np.sin(theta)
        cos_theta = np.cos(theta)

        top_s = sin_lat * cos_theta * rx + sin_lat * sin_theta * ry - cos_lat * rz
        top_e = -sin_theta * rx + cos_theta * ry
        top_z = cos_lat * cos_theta * rx + cos_lat * sin_theta * ry + sin_lat * rz

        az = np.arctan(-top_e / top_s)
        if top_s > 0:
            az = az + np.pi
        if az < 0:
            az = az + 2 * np.pi

        rg = np.sqrt(rx * rx + ry * ry + rz * rz)
        el = np.arcsin(top_z / rg)
        w = (rx * rvx + ry * rvy + rz * rvz) / rg

        return np.rad2deg(az), np.rad2deg(el)

    def get_orbit_number(self, utc_time, tbus_style=False):
        """Calculate orbit number at specified time.
        Optionally use TBUS-style orbit numbering (TLE orbit number + 1)
        """
        
        dt = astronomy._days(utc_time - self.orbit_elements.an_time)
        orbit_period = astronomy._days(self.orbit_elements.an_period)
                 
        orbit = int(self.tle.orbit + dt / orbit_period + 
                 self.tle.mean_motion_derivative * dt**2 + 
                 self.tle.mean_motion_sec_derivative * dt**3)
                 
        if tbus_style:
            orbit += 1
        return orbit
        
    def get_zenith_overpass(self, utc_time, obslon, obslat):
        """Get the time when the satellite is highest on the horizon relative
        to the observer position on ground given by *obslon*,*obslat* closest
        in time (in the future) to the time given by *utc_time*. Using the
        method of gradient ascent with variable step parameter (decreasing in
        size as we apprach 90 degrees zenith angle.
        """

        # First check if the elevation is above zero. If not continue until it
        # is, using larger steps when the ange is far from zero and shorter
        # ones when the elevation gets closer to zero:
        el_start = self.get_observer_look(utc_time, obslon, obslat, 0.0)[1]
        if el_start < 0:
            start_time = utc_time
            elev = el_start
            idx = 0
            NIDX = 100
            while elev < 0 and idx < NIDX:
                var_scale = np.exp(np.square(np.sin(elev * np.pi/180.)))
                t_step = timedelta(seconds = (300 * var_scale))
                start_time = start_time + t_step
                elev = self.get_observer_look(start_time, obslon, obslat, 0.0)[1]
                idx = idx + 1
                #print idx, start_time, var_scale, elev

            utc_time = start_time - t_step

        precision = timedelta(seconds=0.01)
        sec_step = 0.5
        t_step = timedelta(seconds=sec_step/2.0)

        # Local derivative:
        def fprime(timex):
            el0 = self.get_observer_look(timex - t_step, 
                                         obslon, obslat, 0.0)[1]
            el1 = self.get_observer_look(timex + t_step, 
                                         obslon, obslat, 0.0)[1]
            return el0, (el1 - el0) / sec_step 

        tx0 = utc_time - timedelta(seconds=1.0)
        tx1 = utc_time
        idx = 0
        NIDX = 100
        eps = 1000. # Step size scaling
        while abs(tx1 - tx0) > precision and idx < NIDX:
            tx0 = tx1
            fpr = fprime(tx0)
            #var_scale = abs(90.0-fpr[0])
            var_scale = np.abs(np.cos(fpr[0] * np.pi/180.))
            tx1 = tx0 + timedelta(seconds = (eps * var_scale * fpr[1]))
            #print idx, tx0, tx1, fpr
            idx = idx + 1
    
        if abs(tx1 - tx0) <= precision and idx < NIDX:
            return tx1, idx
        else:
            return None
        
    def get_risetime(self, utc_time, obslon, obslat, **kwargs):
        """Get the risetime of the satellite closest in time in the future to
        the time *utc_time* for a reception station at *obslon*, *obslat*
        position on the ground
        """
        retv = self.get_zenith_overpass(utc_time, obslon, obslat)
        if retv:
            zenith_time = retv[0]
        else:
            zenith_time = utc_time
        one_minute = timedelta(seconds = 60)
        return self._get_time_at_horizon(zenith_time - one_minute, 
                                         obslon, obslat, **kwargs)
        
    def get_falltime(self, utc_time, obslon, obslat, **kwargs):
        """Get the falltime of the satellite closest in time in the future to
        the time *utc_time* for a reception station at *obslon*, *obslat*
        position on the ground
        """
        retv = self.get_zenith_overpass(utc_time, obslon, obslat)
        if retv:
            zenith_time = retv[0]
        else:
            zenith_time = utc_time
        one_minute = timedelta(seconds = 60)
        return self._get_time_at_horizon(zenith_time + one_minute, 
                                         obslon, obslat, **kwargs)
        
    def _get_time_at_horizon(self, utc_time, obslon, obslat, **kwargs):
        """Get the time closest in time to *utc_time* when the
        satellite is at the horizon relative to the position of an observer on
        ground (altitude = 0)
        """
        if "precision" in kwargs:
            precision = kwargs['precision']
        else:
            precision = timedelta(seconds=0.001)
        if "max_iterations" in kwargs:
            nmax_iter = kwargs["max_iterations"]
        else:
            nmax_iter = 100

        sec_step = 0.5
        t_step = timedelta(seconds=sec_step/2.0)

        # Local derivative:
        def fprime(timex):
            el0 = self.get_observer_look(timex - t_step, 
                                         obslon, obslat, 0.0)[1]
            el1 = self.get_observer_look(timex + t_step, 
                                         obslon, obslat, 0.0)[1]
            return el0, (abs(el1) - abs(el0)) / sec_step 

        tx0 = utc_time - timedelta(seconds=1.0)
        tx1 = utc_time
        idx = 0
        #eps = 500.
        eps = 100.
        while abs(tx1 - tx0) > precision and idx < nmax_iter:
            tx0 = tx1
            fpr = fprime(tx0)
            # When the elevation is high the scale is high, and when
            # the elevation is low the scale is low
            #var_scale = np.abs(np.sin(fpr[0] * np.pi/180.))
            #var_scale = np.sqrt(var_scale)
            var_scale = np.abs(fpr[0])
            tx1 = tx0 - timedelta(seconds = (eps * var_scale * fpr[1]))
            idx = idx + 1
            #print idx, tx0, tx1, var_scale, fpr
            if abs(tx1 - utc_time) < precision and idx < 2:
                tx1 = tx1 + timedelta(seconds=1.0)
                
        if abs(tx1 - tx0) <= precision and idx < nmax_iter:
            return tx1
        else:
            return None


class OrbitElements(object):
    """Class holding the orbital elements.
    """
    
    def __init__(self, tle):
        self.epoch = tle.epoch
        self.excentricity = tle.excentricity
        self.inclination = np.deg2rad(tle.inclination)
        self.right_ascension = np.deg2rad(tle.right_ascension)
        self.arg_perigee = np.deg2rad(tle.arg_perigee)
        self.mean_anomaly = np.deg2rad(tle.mean_anomaly)

        self.mean_motion = tle.mean_motion * (np.pi * 2 / XMNPDA)
        self.mean_motion_derivative = tle.mean_motion_derivative * np.pi * 2 / XMNPDA ** 2
        self.mean_motion_sec_derivative = tle.mean_motion_sec_derivative * np.pi * 2 / XMNPDA ** 3
        self.bstar = tle.bstar * AE

        n_0 = self.mean_motion
        k_e = XKE
        k_2 = CK2
        i_0 = self.inclination
        e_0 = self.excentricity

        a_1 = (k_e / n_0) ** (2.0/3)
        delta_1 = ((3/2.0) * (k_2 / a_1**2) * ((3 * np.cos(i_0)**2 - 1) /
                                              (1 - e_0**2)**(2.0/3)))

        a_0 = a_1 * (1 - delta_1/3 - delta_1**2 - (134.0/81) * delta_1**3)

        delta_0 = ((3/2.0) * (k_2 / a_0**2) * ((3 * np.cos(i_0)**2 - 1) /
                                              (1 - e_0**2)**(2.0/3)))

        # original mean motion
        n_0pp = n_0 / (1 + delta_0)
        self.original_mean_motion = n_0pp

        # semi major axis
        a_0pp = a_0 / (1 - delta_0)
        self.semi_major_axis = a_0pp

        self.period = np.pi * 2 / n_0pp

        self.perigee = (a_0pp * (1 - e_0) / AE - AE) * XKMPER

        self.right_ascension_lon = (self.right_ascension
                                           - astronomy.gmst(self.epoch))

        if self.right_ascension_lon > np.pi:
            self.right_ascension_lon -= 2 * np.pi    


class _SGDP4(object):
    """Class for the SGDP4 computations.
    """

    def __init__(self, orbit_elements):
        self.mode = None

        perigee = orbit_elements.perigee
        self.eo = orbit_elements.excentricity
        print 'eo', self.eo
        self.xincl = orbit_elements.inclination
        self.xno = orbit_elements.original_mean_motion
        k_2 = CK2
        k_4 = CK4
        k_e = XKE
        self.bstar = orbit_elements.bstar
        self.omegao = orbit_elements.arg_perigee
        self.xmo = orbit_elements.mean_anomaly
        self.xnodeo = orbit_elements.right_ascension
        self.t_0 = orbit_elements.epoch
        self.xn_0 = orbit_elements.mean_motion
        A30 = -XJ3 * AE**3

        if not(0 <= self.eo < ECC_LIMIT_HIGH):
            raise OrbitalError('Eccentricity out of range: %e' % self.eo)
        elif not((0.0035 * 2 * np.pi / XMNPDA) < self.xn_0 < (18 * 2 * np.pi / XMNPDA)):
            raise OrbitalError('Mean motion out of range: %e' % self.xn_0)
        elif not(0 < self.xincl < np.pi):
            raise OrbitalError('Inclination out of range: %e' % self.xincl)

        self.cosIO = np.cos(self.xincl)
        self.sinIO = np.sin(self.xincl)
        theta2 = self.cosIO**2
        theta4 = theta2 ** 2 
        self.x3thm1 = 3.0 * theta2 - 1.0
        self.x1mth2 = 1.0 - theta2
        self.x7thm1 = 7.0 * theta2 - 1.0

        a1 = (XKE / self.xn_0) ** (2. / 3)
        betao2 = 1.0 - self.eo**2
        betao = np.sqrt(betao2)
        temp0 = 1.5 * CK2 * self.x3thm1 / (betao * betao2)
        del1 = temp0 / (a1**2)
        a0 = a1 * (1.0 - del1 * (1.0 / 3.0 + del1 * (1.0 + del1 * 134.0 / 81.0)))
        del0 = temp0 / (a0**2)
        self.xnodp = self.xn_0 / (1.0 + del0)
        self.aodp = (a0 / (1.0 - del0))
        self.perigee = (self.aodp * (1.0 - self.eo) - AE) * XKMPER
        self.apogee = (self.aodp * (1.0 + self.eo) - AE) * XKMPER
        self.period = (2 * np.pi * 1440.0 / XMNPDA) / self.xnodp 
     
        if self.eo == 0:
            self.mode = self.SGDP4_ZERO_ECC
            return     
     
        if self.period >= 225:
            # Deep-Space model
            self.mode = SGDP4_DEEP_NORM
        elif self.perigee < 220:
            # Near-space, simplified equations
            self.mode = SGDP4_NEAR_SIMP
        else:
            # Near-space, normal equations
            self.mode = SGDP4_NEAR_NORM

        if self.perigee < 156:
            s4 = self.perigee - 78
            if s4 < 20:
                s4 = 20
            
            qoms24 = ((120 - s4) * (AE / XKMPER))**4
            s4 = (s4 / XKMPER + AE)
        else:
            s4 = KS
            qoms24 = QOMS2T

        pinvsq = 1.0 / (self.aodp**2 * betao2**2)
        tsi = 1.0 / (self.aodp - s4)
        self.eta = self.aodp * self.eo * tsi
        etasq = self.eta**2
        eeta = self.eo * self.eta
        psisq = np.abs(1.0 - etasq)
        coef = qoms24 * tsi**4
        coef_1 = coef / psisq**3.5

        self.c2 = (coef_1 * self.xnodp * (self.aodp *
             (1.0 + 1.5 * etasq + eeta * (4.0 + etasq)) +
             (0.75 * CK2) * tsi / psisq * self.x3thm1 *
             (8.0 + 3.0 * etasq * (8.0 + etasq))))

        self.c1 = self.bstar * self.c2

        self.c4 = (2.0 * self.xnodp * coef_1 * self.aodp * betao2 * (self.eta *
             (2.0 + 0.5 * etasq) + self.eo * (0.5 + 2.0 *
             etasq) - (2.0 * CK2) * tsi / (self.aodp * psisq) * (-3.0 *
             self.x3thm1 * (1.0 - 2.0 * eeta + etasq *
             (1.5 - 0.5 * eeta)) + 0.75 * self.x1mth2 * (2.0 *
             etasq - eeta * (1.0 + etasq)) * np.cos(2.0 * self.omegao))))

        self.c5, self.c3, self.omgcof = 0.0, 0.0, 0.0


        if self.mode == SGDP4_NEAR_NORM:
            self.c5 = (2.0 * coef_1 * self.aodp * betao2 *
             (1.0 + 2.75 * (etasq + eeta) + eeta * etasq))
            if self.eo > ECC_ALL:
                self.c3 = coef * tsi * A3OVK2 * self.xnodp * AE * self.sinIO / self.eo
            self.omgcof = self.bstar * self.c3 * np.cos(self.omegao)

        temp1 = 3.0 * CK2 * pinvsq * self.xnodp
        temp2 = temp1 * CK2 * pinvsq
        temp3 = 1.25 * CK4 * pinvsq**2 * self.xnodp

        self.xmdot = (self.xnodp + (0.5 * temp1 * betao * self.x3thm1 + 0.0625 *
                temp2 * betao * (13.0 - 78.0 * theta2 +
                137.0 * theta4)))

        x1m5th = 1.0 - 5.0 * theta2

        self.omgdot = (-0.5 * temp1 * x1m5th + 0.0625 * temp2 *
                 (7.0 - 114.0 * theta2 + 395.0 * theta4) +
                 temp3 * (3.0 - 36.0 * theta2 + 49.0 * theta4))

        xhdot1 = -temp1 * self.cosIO
        self.xnodot = (xhdot1 + (0.5 * temp2 * (4.0 - 19.0 * theta2) +
                 2.0 * temp3 * (3.0 - 7.0 * theta2)) * self.cosIO)

        if self.eo > ECC_ALL:
            self.xmcof = (-(2. / 3) * AE) * coef * self.bstar / eeta
        else:    
            self.xmcof = 0.0

        self.xnodcf = 3.5 * betao2 * xhdot1 * self.c1
        self.t2cof = 1.5 * self.c1
        
        # Check for possible divide-by-zero for X/(1+cos(xincl)) when calculating xlcof */
        temp0 = 1.0 + self.cosIO
        if np.abs(temp0) < EPS_COS:
            temp0 = np.sign(temp0) * EPS_COS
    	    
        self.xlcof = 0.125 * A3OVK2 * self.sinIO * (3.0 + 5.0 * self.cosIO) / temp0

        self.aycof = 0.25 * A3OVK2 * self.sinIO
        
        self.cosXMO = np.cos(self.xmo)
        self.sinXMO = np.sin(self.xmo)
        self.delmo = (1.0 + self.eta * self.cosXMO)**3
        
        self._deep_space = None
        
        if self.mode == SGDP4_NEAR_NORM:        
            c1sq = self.c1**2
            self.d2 = 4.0 * self.aodp * tsi * c1sq
            temp0 = self.d2 * tsi * self.c1 / 3.0
            self.d3 = (17.0 * self.aodp + s4) * temp0
            self.d4 = 0.5 * temp0 * self.aodp * tsi * (221.0 * self.aodp + 31.0 * s4) * self.c1
            self.t3cof = self.d2 + 2.0 * c1sq
            self.t4cof = 0.25 * (3.0 * self.d3 + self.c1 * (12.0 * self.d2 + 10.0 * c1sq))
            self.t5cof = (0.2 * (3.0 * self.d4 + 12.0 * self.c1 * self.d3 + 6.0 * self.d2**2 + 
                    15.0 * c1sq * (2.0 * self.d2 + c1sq)))

        elif self.mode == SGDP4_DEEP_NORM:
            self._deep_space = _DeepSpace(self.t_0, self.omegao, self.xnodeo, 
                                          self.xmo, self.eo, self.xincl, 
                                          self.aodp, self.xmdot, 
                                          self.omgdot, self.xnodot, 
                                          self.xnodp)
            #raise NotImplementedError('Deep space calculations not supported')
        
    def propagate(self, utc_time):
        kep = {}

        ts = astronomy._days(utc_time - self.t_0) * XMNPDA
        print 'ts', ts
        em = self.eo
        xinc = self.xincl
        
        xmp   = self.xmo + self.xmdot * ts
        xnode = self.xnodeo + ts * (self.xnodot + ts * self.xnodcf)
        omega = self.omegao + self.omgdot * ts
        
        if self.mode == SGDP4_ZERO_ECC:
            raise NotImplementedError('Mode SGDP4_ZERO_ECC not implemented')
        elif self.mode == SGDP4_NEAR_SIMP:
            raise NotImplementedError('Mode "Near-space, simplified equations"'
                                      ' not implemented')
        elif self.mode == SGDP4_NEAR_NORM:
            delm  = self.xmcof * ((1.0 + self.eta * np.cos(xmp))**3 - self.delmo)
            temp0 = ts * self.omgcof + delm
            xmp += temp0
            omega -= temp0
            tempa = 1.0 - (ts * (self.c1 + ts * (self.d2 + ts * (self.d3 + ts * self.d4))))
            tempe = self.bstar * (self.c4 * ts + self.c5 * (np.sin(xmp) - self.sinXMO))
            templ = ts * ts * (self.t2cof + ts * (self.t3cof + ts * (self.t4cof + ts * self.t5cof)))
            a = self.aodp * tempa**2
            e = em - tempe
            xl = xmp + omega + xnode + self.xnodp * templ

            x3thm1 = self.x3thm1
            x1mth2 = self.x1mth2
            x7thm1 = self.x7thm1
            xlcof = self.xlcof
            aycof = self.aycof
            sinIO = self.sinIO
            cosIO = self.cosIO
            
        else:
            #raise  NotImplementedError('Deep space calculations not supported')
            tempa = 1.0 - ts * self.c1
            tempe = self.bstar * ts * self.c4
            templ = ts * ts * self.t2cof
            xn = self.xnodp
            xmp, omega, xnode, em, xinc, xn = \
                self._deep_space.dpsec(xmp, omega, xnode, em, xinc, xn, ts)
            print 'xmp', xmp, 'omega', omega, 'xnode', xnode, 'em', em, 'xinc', xinc, 'xn', xn
            a = (XKE / xn) ** (2. / 3.) * tempa ** 2
            e = em - tempe
            xmam = xmp + self.xnodp * templ
            print 'xnodp', self.xnodp, 'templ', templ
            print 'a', a, 'e', e, 'xmam', xmam
            e, xinc, omega, xnode, xmam = \
                self._deep_space.dpper(e, xinc, omega, xnode, xmam, ts)
            print 'dpper'
            print 'e', e, 'xinc', xinc, 'omega', omega, 'xnode', xnode, 'xmam', xmam
            
            if xinc < 0:
                xinc = -xinc
                xnode += np.pi
                omega -= np.pi
                
            xl = xmam + omega + xnode

            # Re-compute the perturbed values            
            sinIO = np.sin(xinc)
            cosIO = np.cos(xinc)
            theta2 = cosIO * cosIO
            x3thm1 = 3.0 * theta2 - 1.0
            x1mth2 = 1.0 - theta2
            x7thm1 = 7.0 * theta2 - 1.0

            # Check for possible divide-by-zero for X/(1+cosIO) when calculating xlcof            
            temp0 = 1.0 + cosIO
            temp0 = np.where(np.abs(temp0) < EPS_COSIO, np.sign(temp0) * EPS_COSIO, temp0)

            xlcof = 0.125 * A3OVK2 * sinIO * (3.0 + 5.0 * cosIO) / temp0
            aycof = 0.25 * A3OVK2 * sinIO
            
        if np.any(a < 1):
            raise Exception('Satellite crased at time %s', utc_time)
        elif np.any(e < ECC_LIMIT_LOW):
            raise ValueError('Satellite modified eccentricity to low: %e < %e'
                             % (e, ECC_LIMIT_LOW))

        e = np.where(e < ECC_EPS, ECC_EPS, e)
        e = np.where(e > ECC_LIMIT_HIGH, ECC_LIMIT_HIGH, e)
            
        beta2 = 1.0 - e ** 2
        
        # Long period periodics
        sinOMG = np.sin(omega)
        cosOMG = np.cos(omega) 

        temp0 = 1.0 / (a * beta2)
        axn = e * cosOMG
        ayn = e * sinOMG + temp0 * aycof
        print 'axn', axn, 'ayn', ayn
        xlt = xl + temp0 * xlcof * axn

        elsq = axn ** 2 + ayn ** 2
        print 'elsq', elsq
        
        if np.any(elsq >= 1):
            raise Exception('e**2 >= 1 at %s', utc_time)
            
        kep['ecc'] = np.sqrt(elsq)
        
        epw = np.fmod(xlt - xnode, 2 * np.pi)
        print 'pre epw', epw
        # needs a copy in case of an array
        capu = np.array(epw)
        maxnr = kep['ecc']
        for i in range(10):
            sinEPW = np.sin(epw)
            cosEPW = np.cos(epw)

            ecosE = axn * cosEPW + ayn * sinEPW
            esinE = axn * sinEPW - ayn * cosEPW
            f = capu - epw + esinE
            if np.all(np.abs(f) < NR_EPS):
                break
                
            df = 1.0 - ecosE

            # 1st order Newton-Raphson correction. 
            nr = f / df
            
            # 2nd order Newton-Raphson correction.
            nr = np.where(np.logical_and(i == 0, np.abs(nr) > 1.25 * maxnr),
                          np.sign(nr) * maxnr,
                          f / (df + 0.5*esinE*nr))
            epw += nr
            
        print 'post epw', epw    
        print 'sinEPW', sinEPW
        print 'cosEPW', cosEPW
        print 'ecosE', ecosE
        print 'axn * cosEPW', axn * cosEPW
        print 'ayn * sinEPW', ayn * sinEPW
        
        # Short period preliminary quantities 
        temp0 = 1.0 - elsq
        betal = np.sqrt(temp0)
        pl = a * temp0
        r = a * (1.0 - ecosE)
        print 'r', r
        invR = 1.0 / r
        temp2 = a * invR
        print 'temp2', temp2
        temp3 = 1.0 / (1.0 + betal)
        print 'temp3', temp3
        cosu = temp2 * (cosEPW - axn + ayn * esinE * temp3)
        sinu = temp2 * (sinEPW - ayn - axn * esinE * temp3)
        print 'cosu', cosu, 'sinu', sinu  
        u = np.arctan2(sinu, cosu)
        sin2u = 2.0 * sinu * cosu
        cos2u = 2.0 * cosu**2 - 1.0
        temp0 = 1.0 / pl
        temp1 = CK2 * temp0
        temp2 = temp1 * temp0
        
        # Update for short term periodics to position terms. 

        rk = r * (1.0 - 1.5 * temp2 * betal * x3thm1) + 0.5 * temp1 * x1mth2 * cos2u
        uk = u - 0.25 * temp2 * x7thm1 * sin2u
        print 'u', u, 'temp2', temp2, 'x7thm1', x7thm1, 'sin2u', sin2u
        xnodek = xnode + 1.5 * temp2 * cosIO * sin2u
        xinck = xinc + 1.5 * temp2 * cosIO * sinIO * cos2u
        
        if np.any(rk < 1):
            raise Exception('Satellite crashed at time %s', utc_time)
        
        temp0 = np.sqrt(a)
        temp2 = XKE / (a * temp0)
        rdotk = ((XKE * temp0 * esinE * invR -temp2 * temp1 * x1mth2 * sin2u) * 
                 (XKMPER / AE * XMNPDA / 86400.0))
        rfdotk = ((XKE * np.sqrt(pl) * invR + temp2 * temp1 * 
                  (x1mth2 * cos2u + 1.5 * x3thm1)) * 
                  (XKMPER / AE * XMNPDA / 86400.0))
            
        kep['radius'] = rk * XKMPER / AE
        kep['theta'] = uk
        kep['eqinc'] = xinck
        kep['ascn'] = xnodek 
        kep['argp'] = omega 
        kep['smjaxs'] = a * XKMPER / AE  
        kep['rdotk'] = rdotk 
        kep['rfdotk'] = rfdotk

        print 'kep:', kep
        return kep   


class _DeepSpace(object):
    
    def __init__(self, epoch, omegao, xnodeo, xmo, orb_eo, orb_xincl, aodp, 
                 xlldot, omgdot, xnodot, xnodp):

        def mod2pi(a):
            b = np.fmod(a, 2 * np.pi)
            if b < 0:
                return b + 2 * np.pi
            else:
                return b

        eo = eq = orb_eo
        xincl = orb_xincl
        
        # Decide on direct or Lyddane Lunar-Solar perturbations. 
        self.ilsd = False
        if xincl >= 0.2:
            self.ilsd = True
            
        # Drop some terms below 3 deg inclination.
        self.ishq = False
        if xincl >= 0.052359877:
            self.ishq = True
            
        sinomo = np.sin(omegao)
        cosomo = np.cos(omegao)
        sinq = np.sin(xnodeo)
        cosq = np.cos(xnodeo)
        siniq = np.sin(xincl)
        cosiq = np.cos(xincl)
        
        if np.abs(siniq) <= EPS_SIN:
            siniq = np.sign(siniq) * EPS_SIN
            
        cosiq2 = cosiq ** 2
        siniq2 = siniq ** 2
        print cosiq2, siniq2
        
        ao = aodp
        self.omgdt = omgdot
        eqsq = eo ** 2
        bsq = 1 - eqsq
        rteqsq = np.sqrt(bsq)
        print 'rteqsq', rteqsq
        self.thgr = astronomy.gmst(epoch)
        print 'thgr', self.thgr
        
        xnq = xnodp
        aqnv = 1. / ao
        xmao = xmo
        xpidot = self.omgdt + xnodot
        print 'xpidot', xpidot
        self.omegaq = omegao
        
        # Initialize lunar terms.
        # Note: the Dundee reference code d50 is off by 1 day
        jd1900 = astronomy.jdays1900(epoch)
        jd1900 += 1 #TODO: compability?
        print 'epoch', epoch
        print 'd1900', jd1900
        xnodce = 4.523602 - jd1900 * 9.2422029e-4
        temp0 = np.fmod(xnodce, 2 * np.pi)
        stem = np.sin(temp0)
        ctem = np.cos(temp0)
        
        zcosil = 0.91375164 - ctem * 0.03568096
        zsinil = np.sqrt(1.0 - zcosil * zcosil)
        zsinhl = stem * 0.089683511 / zsinil
        zcoshl = np.sqrt(1.0 - zsinhl * zsinhl)
        print 'zcoshl', zcoshl
        
        c = jd1900 * 0.2299715 + 4.7199672
        gam = jd1900 * 0.001944368 + 5.8351514
        print 'c', c
        print 'gam', gam
        self.zmol = mod2pi(c - gam)
        print 'zmol',  self.zmol
        zx = stem * 0.39785416 / zsinil
        zy = zcoshl * ctem + zsinhl * 0.91744867 * stem
        zx = np.arctan2(zx, zy)
        zx = np.fmod(gam + zx - xnodce, 2 * np.pi)
        print 'zx', zx
        zsingl = np.sin(zx)
        zcosgl = np.cos(zx)
        self.zmos = mod2pi(jd1900 * 0.017201977 + 6.2565837)
        print 'zsingl', zsingl 
        print 'zcosgl', zcosgl
        print 'zmos', self.zmos
        
        # Do solar terms
        zcosg = ZCOSGS
        zsing = ZSINGS
        zcosi = 0.91744867
        zsini = 0.39785416
        zcosh = cosq
        zsinh = sinq
        cc = C1SS
        zn = ZNS
        ze = ZES
        zmo = self.zmos
        xnoi = 1. / xnq
        
        for ls in range(2):
            a1 = zcosg * zcosh + zsing * zcosi * zsinh
            a3 = -zsing * zcosh + zcosg * zcosi * zsinh
            a7 = -zcosg * zsinh + zsing * zcosi * zcosh
            a8 = zsing * zsini
            a9 = zsing * zsinh + zcosg * zcosi * zcosh
            a10 = zcosg * zsini
            a2 = cosiq * a7 + siniq * a8
            a4 = cosiq * a9 + siniq * a10
            a5 = -siniq * a7 + cosiq * a8
            a6 = -siniq * a9 + cosiq * a10
    
            x1 = a1 * cosomo + a2 * sinomo
            x2 = a3 * cosomo + a4 * sinomo
            x3 = -a1 * sinomo + a2 * cosomo
            x4 = -a3 * sinomo + a4 * cosomo
            x5 = a5 * sinomo
            x6 = a6 * sinomo
            x7 = a5 * cosomo
            x8 = a6 * cosomo
            
            z31 = x1 * 12.0 * x1 - x3 * 3.0 * x3
            z32 = x1 * 24.0 * x2 - x3 * 6.0 * x4
            z33 = x2 * 12.0 * x2 - x4 * 3.0 * x4
            z1 = (a1 * a1 + a2 * a2) * 3.0 + z31 * eqsq
            z2 = (a1 * a3 + a2 * a4) * 6.0 + z32 * eqsq
            z3 = (a3 * a3 + a4 * a4) * 3.0 + z33 * eqsq
            z11 = a1 * -6.0 * a5 + eqsq * (x1 * -24.0 * x7 - x3 * 6.0 * x5)
            z12 = ((a1 * a6 + a3 * a5) * -6.0 + eqsq * ((x2 * x7 
                   + x1 * x8) * -24.0 - (x3 * x6 + x4 * x5) * 6.0))
            z13 = a3 * -6.0 * a6 + eqsq * (x2 * -24.0 * x8 - x4 * 6.0 * x6)
            z21 = a2 * 6.0 * a5 + eqsq * (x1 * 24.0 * x5 - x3 * 6.0 * x7)
            z22 = ((a4 * a5 + a2 * a6) * 6.0 + eqsq * ((x2 * x5 + x1 * x6) 
                    * 24.0 - (x4 * x7 + x3 * x8) * 6.0))
            z23 = a4 * 6.0 * a6 + eqsq * (x2 * 24.0 * x6 - x4 * 6.0 * x8)
            z1 = z1 + z1 + bsq * z31
            z2 = z2 + z2 + bsq * z32
            z3 = z3 + z3 + bsq * z33
            s3 = cc * xnoi
            s2 = s3 * -0.5 / rteqsq
            s4 = s3 * rteqsq
            s1 = eq * -15.0 * s4
            s5 = x1 * x3 + x2 * x4
            s6 = x2 * x3 + x1 * x4
            s7 = x2 * x4 - x1 * x3
            se = s1 * zn * s5
            si = s2 * zn * (z11 + z13)
            sl = -zn * s3 * (z1 + z3 - 14.0 - eqsq * 6.0)
            sgh = s4 * zn * (z31 + z33 - 6.0)
            print 'sgh', sgh
            
            shdq = 0
            if self.ishq:
                sh = -zn * s2 * (z21 + z23)
                shdq = sh / siniq
                
            self.ee2 = s1 * 2.0 * s6
            self.e3 = s1 * 2.0 * s7
            self.xi2 = s2 * 2.0 * z12
            self.xi3 = s2 * 2.0 * (z13 - z11)
            self.xl2 = s3 * -2.0 * z2
            self.xl3 = s3 * -2.0 * (z3 - z1)
            self.xl4 = s3 * -2.0 * (-21.0 - eqsq * 9.0) * ze
            self.xgh2 = s4 * 2.0 * z32
            self.xgh3 = s4 * 2.0 * (z33 - z31)
            self.xgh4 = s4 * -18.0 * ze
            self.xh2 = s2 * -2.0 * z22
            self.xh3 = s2 * -2.0 * (z23 - z21)
            print 'xh3', self.xh3
            
            # Break on second iteration
            if ls == 1:
                break
            
            # Do lunar terms
            self.sse = se
            self.ssi = si
            self.ssl = sl
            self.ssh = shdq
            self.ssg = sgh - cosiq * self.ssh
            self.se2 = self.ee2
            self.si2 = self.xi2
            self.sl2 = self.xl2
            self.sgh2 = self.xgh2
            self.sh2 = self.xh2
            self.se3 = self.e3
            self.si3 = self.xi3
            self.sl3 = self.xl3
            self.sgh3 = self.xgh3
            self.sh3 = self.xh3
            self.sl4 = self.xl4
            self.sgh4 = self.xgh4
            zcosg = zcosgl
            zsing = zsingl
            zcosi = zcosil
            zsini = zsinil;
            zcosh = zcoshl * cosq + zsinhl * sinq
            zsinh = sinq * zcoshl - cosq * zsinhl
            zn = ZNL
            cc = C1L
            ze = ZEL
            zmo = self.zmol
        
        self.sse += se
        self.ssi += si
        self.ssl += sl
        self.ssg += sgh - cosiq * shdq
        self.ssh += shdq
        print 'ssg', self.ssg
        
        if 0.0034906585 < xnq < 0.0052359877:
            raise NotImplementedError('24h resonance not implemented')
        elif 0.00826 <= xnq <= 0.00924 and eq >= 0.5:
            #raise NotImplementedError('12h resonance not implemented')
            self.iresfl = 1
            self.isynfl = 0
            eoc = eq * eqsq
            g201 = -0.306 - (eq - 0.64) * 0.44

            if eq <= 0.65:
                g211 = 3.616 - eq * 13.247 + eqsq * 16.29
                g310 = eq * 117.39 - 19.302 - eqsq * 228.419 + eoc * 156.591
                g322 = eq * 109.7927 - 18.9068 - eqsq * 214.6334 + eoc * 146.5816
                g410 = eq * 242.694 - 41.122 - eqsq * 471.094 + eoc * 313.953
                g422 = eq * 841.88 - 146.407 - eqsq * 1629.014 + eoc * 1083.435
                g520 = eq * 3017.977 - 532.114 - eqsq * 5740.032 + eoc * 3708.276
            else:
                g211 = eq * 331.819 - 72.099 - eqsq * 508.738 + eoc * 266.724
                g310 = eq * 1582.851 - 346.844 - eqsq * 2415.925 + eoc * 1246.113
                g322 = eq * 1554.908 - 342.585 - eqsq * 2366.899 + eoc * 1215.972
                g410 = eq * 4758.686 - 1052.797 - eqsq * 7193.992 + eoc * 3651.957
                g422 = eq * 16178.11 - 3581.69 - eqsq * 24462.77 + eoc * 12422.52
    
                if eq <= 0.715:
                    g520 = 1464.74 - eq * 4664.75 + eqsq * 3763.64
                else:
                    g520 = eq * 29936.92 - 5149.66 - eqsq * 54087.36 + eoc * 31324.56
    
            if eq < 0.7:
                g533 = eq * 4988.61 - 919.2277 - eqsq * 9064.77 + eoc * 5542.21
                g521 = eq * 4568.6173 - 822.71072 - eqsq * 8491.4146 + eoc * 5337.524
                g532 = eq * 4690.25 - 853.666 - eqsq * 8624.77 + eoc * 5341.4
            else:
                g533 = eq * 161616.52 - 37995.78 - eqsq * 229838.2 + eoc * 109377.94
                g521 = eq * 218913.95 - 51752.104 - eqsq * 309468.16 + eoc * 146349.42
                g532 = eq * 170470.89 - 40023.88 - eqsq * 242699.48 + eoc * 115605.82

            f220 = (cosiq * 2.0 + 1.0 + cosiq2) * 0.75
            f221 = siniq2 * 1.5
            f321 = siniq * 1.875 * (1.0 - cosiq * 2.0 - cosiq2 * 3.0)
            f322 = siniq * -1.875 * (cosiq * 2.0 + 1.0 - cosiq2 * 3.0)
            f441 = siniq2 * 35.0 * f220
            f442 = siniq2 * 39.375 * siniq2
            f522 = (siniq * 9.84375 * (siniq2 * (1.0 - cosiq *
                    2.0 - cosiq2 * 5.0) + (cosiq * 4.0 - 2.0 + cosiq2 * 6.0) * 0.33333333))
            f523 = (siniq * (siniq2 * 4.92187512 * (-2.0 - cosiq *
                    4.0 + cosiq2 * 10.0) + (cosiq * 2.0 +
                    1.0 - cosiq2 * 3.0) * 6.56250012))
            f542 = (siniq * 29.53125 * (2.0 - cosiq * 8.0 +
                    cosiq2 * (cosiq * 8.0 - 12.0 + cosiq2 * 10.0)))
            f543 = (siniq * 29.53125 * (-2.0 - cosiq * 8.0 +
                    cosiq2 * (cosiq * 8.0 + 12.0 - cosiq2 * 10.0)))
            xno2 = xnq * xnq
            ainv2 = aqnv * aqnv
            temp1 = xno2 * 3.0 * ainv2
            temp0 = temp1 * ROOT22
            self.d2201 = temp0 * f220 * g201
            self.d2211 = temp0 * f221 * g211
            temp1 *= aqnv
            temp0 = temp1 * ROOT32
            self.d3210 = temp0 * f321 * g310
            self.d3222 = temp0 * f322 * g322
            temp1 *= aqnv
            temp0 = temp1 * 2.0 * ROOT44
            self.d4410 = temp0 * f441 * g410
            self.d4422 = temp0 * f442 * g422
            temp1 *= aqnv
            temp0 = temp1 * ROOT52
            self.d5220 = temp0 * f522 * g520
            self.d5232 = temp0 * f523 * g532
            temp0 = temp1 * 2.0 * ROOT54
            self.d5421 = temp0 * f542 * g521
            self.d5433 = temp0 * f543 * g533
            xlamo = xmao + xnodeo + xnodeo - self.thgr - self.thgr
            bfact = xlldot + xnodot + xnodot - THDT - THDT
            bfact += self.ssl + self.ssh + self.ssh
            
        else:
            self.iresfl = False
            self.isynfl = False
            
        if not self.iresfl:
            self.mode = SGDP4_DEEP_NORM
        else:
            #raise NotImplementedError('Only normal deep space init')
            # Integrator initialization
            self.xfact = bfact - xnq
            self.xli = xlamo
            self.xni = xnq
            
            self.xnddt0, self.xndot0, self.xldot0, = self.dot_terms_calculated(0.0, self.xli, self.xni)
            
            if self.isynfl:
                self.mode = SGDP4_DEEP_SYNC
            else:
                self.mode = SGDP4_DEEP_RESN

    def dpsec(self, xmp, omega, xnode, ema, xinca, xna, tsince):
        xll = xmp + self.ssl * tsince
        omgasm = omega + self.ssg * tsince
        xnodes = xnode + self.ssh * tsince
        em = ema + self.sse * tsince
        xinc = xinca + self.ssi * tsince
        
        #if not self.iresfl:
        #    return xll, omgasm, xnodes, em, xinc, xna
        #else:
        if self.iresfl:
            #raise NotImplementedError('Only normal deep space')
            xni = self.xni
            xli = self.xli
            
            xnddt = self.xnddt0
            xndot = self.xndot0
            xldot = self.xldot0
            
            print 'xni', xni, 'xli', xli
            
            delt = np.sign(tsince) * STEP
            ft = 0
            while np.abs(tsince - ft) >= STEP:
                xli += delt * (xldot + delt * 0.5 * xndot)
                xni += delt * (xndot + delt * 0.5 * xnddt)
                ft += delt
                xnddt, xndot, xldot = self.dot_terms_calculated(ft, xli, xni)
                
            
            print 'integrated xni', xni, 'xli', xli, 'ft', ft
            
            dt = tsince - ft 
            xl = xli + dt * (xldot + dt * 0.5 * xndot)
            xna = xni + dt * (xndot + dt * 0.5 * xnddt)
            
            temp0 = -xnodes + self.thgr + tsince * THDT
            
            if not self.isynfl:
                xll = xl + 2 * temp0
            else:
                xll = xl - omgasm + temp0
                
        return xll, omgasm, xnodes, em, xinc, xna

    def dpper(self, e, xinca, omega, xnode, xmam, ts):
        pgh, ph, pe, pinc, pl = self.compute_lunar_solar(ts)
        print 'pgh', pgh, 'ph', ph, 'pe', pe, 'pinc', pinc, 'pl', pl 
        xinc = xinca + pinc
        em = e + pe
        sinis = np.sin(xinc)
        cosis = np.cos(xinc)
        
        if self.ilsd:
            tmp_ph = ph / sinis
            omgasm = omega + pgh - cosis * tmp_ph
            xnodes = xnode + tmp_ph
            xll = xmam + pl            
        else:
            raise NotImplementedError('Lyddane modifications not implemented')    
        
        return em, xinc, omgasm, xnodes, xll
    
    def compute_lunar_solar(self, tsince):
        
        # Update Solar terms
        zm = self.zmos + ZNS * tsince
        zf = zm + ZES * 2. * np.sin(zm)
        sinzf = np.sin(zf)
        coszf = np.cos(zf)
        f2 = sinzf * 0.5 * sinzf - 0.25
        f3 = sinzf * -0.5 * coszf
        ses  = self.se2 * f2 + self.se3 * f3
        sis  = self.si2 * f2 + self.si3 * f3
        sls  = self.sl2 * f2 + self.sl3 * f3 + self.sl4 * sinzf
        
        sghs = self.sgh2 * f2 + self.sgh3 * f3 + self.sgh4 * sinzf
        shs  = self.sh2  * f2 + self.sh3  * f3
        
        # Update Lunar terms
        zm = self.zmol + ZNL * tsince
        zf = zm + ZEL * 2.0 * np.sin(zm)
        sinzf = np.sin(zf)
        coszf = np.cos(zf)
        f2 = sinzf * 0.5 * sinzf - 0.25
        f3 = sinzf * -0.5 * coszf
        
        sel = self.ee2 * f2 + self.e3 * f3
        sil = self.xi2 * f2 + self.xi3 * f3
        sll = self.xl2 * f2 + self.xl3 * f3 + self.xl4 * sinzf
        
        sghl = self.xgh2 * f2 + self.xgh3 * f3 + self.xgh4 * sinzf
        shl  = self.xh2  * f2 + self.xh3  * f3
        
        pgh = sghs + sghl
        ph = shs + shl
        pe = ses + sel
        pinc = sis + sil
        pl = sls + sll
        
        return pgh, ph, pe, pinc, pl

    def dot_terms_calculated(self, atime, xli, xni):
        if self.isynfl:
            xndot = (self.del1 * np.sin(xli - self.fasx2)
                     + self.del2 * np.sin((xli - self.fasx4) * 2.0)
                     + self.del3 * np.sin((xli - self.fasx6) * 3.0))

            xnddt = (self.del1 * np.cos(xli - self.fasx2)
                     + self.del2 * np.cos((xli - self.fasx4) * 2.0) * 2.0
                     + self.del3 * np.cos((xli - self.fasx6) * 3.0) * 3.0)
        else:
            xomi = self.omegaq + self.omgdt * atime
            x2omi = 2 * xomi
            x2li = 2 * xli

            xndot = (self.d2201 * np.sin(x2omi + xli - G22)
                     + self.d2211 * np.sin(xli - G22)
                     + self.d3210 * np.sin(xomi + xli - G32)
                     + self.d3222 * np.sin(-xomi + xli - G32)
                     + self.d5220 * np.sin(xomi + xli - G52)
                     + self.d5232 * np.sin(-xomi + xli - G52)
                     + self.d4410 * np.sin(x2omi + x2li - G44)
                     + self.d4422 * np.sin(x2li - G44)
                     + self.d5421 * np.sin(xomi + x2li - G54)
                     + self.d5433 * np.sin(-xomi + x2li - G54))

            xnddt = (self.d2201 * np.cos(x2omi + xli - G22)
                     + self.d2211 * np.cos(xli - G22)
                     + self.d3210 * np.cos(xomi + xli - G32)
                     + self.d3222 * np.cos(-xomi + xli - G32)
                     + self.d5220 * np.cos(xomi + xli - G52)
                     + self.d5232 * np.cos(-xomi + xli - G52)
                     + (self.d4410 * np.cos(x2omi + x2li - G44)
                     +  self.d4422 * np.cos(x2li - G44)
                     +  self.d5421 * np.cos(xomi + x2li - G54)
                     +  self.d5433 * np.cos(-xomi + x2li - G54)) * 2.0)

        xldot = xni + self.xfact
        xnddt *= xldot
    
        return xnddt, xndot, xldot         
        
def kep2xyz(kep):
    sinT = np.sin(kep['theta'])
    cosT = np.cos(kep['theta'])
    sinI = np.sin(kep['eqinc'])
    cosI = np.cos(kep['eqinc'])    
    sinS = np.sin(kep['ascn'])
    cosS = np.cos(kep['ascn'])
    
    xmx = -sinS * cosI
    xmy = cosS * cosI

    ux = xmx * sinT + cosS * cosT
    uy = xmy * sinT + sinS * cosT
    uz = sinI * sinT
    
    x = kep['radius'] * ux
    y = kep['radius'] * uy
    z = kep['radius'] * uz
    
    vx = xmx * cosT - cosS * sinT
    vy = xmy * cosT - sinS * sinT
    vz = sinI * cosT
    
    v_x = kep['rdotk'] * ux + kep['rfdotk'] * vx
    v_y = kep['rdotk'] * uy + kep['rfdotk'] * vy
    v_z = kep['rdotk'] * uz + kep['rfdotk'] * vz
    
    return np.array((x, y, z)), np.array((v_x, v_y, v_z))
        
if __name__ == "__main__":
    obs_lon, obs_lat = np.deg2rad((12.4143, 55.9065))
    obs_alt = 0.02
    o = Orbital(satellite="METOP-B")

    t_start = datetime.now()
    t_stop = t_start + timedelta(minutes=20)
    t = t_start
    while t < t_stop:
        t += timedelta(seconds=15)
        lon, lat, alt = o.get_lonlatalt(t)
        lon, lat = np.rad2deg((lon, lat))
        az, el = o.get_observer_look(t, obs_lon, obs_lat, obs_alt)
        ob = o.get_orbit_number(t, tbus_style=True)
        print az, el, ob
