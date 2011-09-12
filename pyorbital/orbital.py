"""
Current Day Number: 11364.541666667 (10 February 2011, 13:00:00 UTC)

// SGP4 test

NORAD Catalog No: 25544
Satellite Name: ISS (ZARYA)
TLE Line 1: 1 25544U 98067A   11036.41885377  .00015023  00000-0  11493-3 0  7736
TLE Line 2: 2 25544  51.6463 116.5688 0003364 277.8443 128.8911 15.72389581700166

SDP4: 0
Period: 91.580357527185
Age: 5.1228128965013
Lat: -51.62128386358
Lon: 18.52912756499
Alt: 367.39281472872
Range: 10318.551827092
Range Rate: 1.2960949888636
Azi: 168.25676645006
Ele: -50.891962381116
Velocity: 7.6898969280099
Visible: 0
Eclipsed: 0
Eclipse Depth: -1.1299921383153
"""
TLE = """1 25544U 98067A   11036.41885377  .00015023  00000-0  11493-3 0  7736
2 25544  51.6463 116.5688 0003364 277.8443 128.8911 15.72389581700166"""

TLE = """1 27424U 02022A   11045.03153664  .00000197  00000-0  53804-4 0  8714
2 27424  98.2146 347.6229 0001550  78.4076 281.7310 14.57100380467139"""

import datetime
import numpy as np
import tlefile
import astronomy


CK2 = 5.413080e-4
CK4 = 0.62098875e-6
E6A = 1.0e-6
QOMS2T = 1.88027916e-9
S = 1.01222928
XJ3 = -0.253881e-5
XKE = 0.743669161e-1
XKMPER = 6378.137
XMNPDA = 1440.0
MFACTOR = 7.292115E-5
AE = 1.0
SECDAY = 8.6400E4
# earth flattening
F = 1/298.257223563

class OrbitalError(Exception):
    pass


class Orbital(object):

    def __init__(self, satellite, tle_file=None):
        satellite = satellite.upper()
        self.satellite_name = satellite
        self.tle = tlefile.read(satellite, tle_file)
        self.orbit_elements = OrbitElements(self.tle)

    def __str__(self):
        return self.satellite_name + " " + str(self.tle)

    def get_position(self, time, normalize=True):
        # for near earth orbits, period must be < 255 minutes

        if self.orbit_elements.period < 255:
            pos, vel = _sgp4(self.orbit_elements, time)
        else:
            raise NotImplementedError, "Currently only handles near earth orbits."

        if not normalize:
            pos = [v*XKMPER for v in pos]
            vel = [v*(XKMPER*XMNPDA/SECDAY) for v in vel]

        return pos, vel

    def get_lonlatalt(self, time):
        (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = self.get_position(time, normalize=True)
        del vel_x, vel_y, vel_z
        lon = ((np.arctan2(pos_y * XKMPER, pos_x * XKMPER) - astronomy.gmst(time))
               % (2 * np.pi))

        if lon > np.pi:
            lon -= np.pi * 2
        if lon <= -np.pi:
            lon += np.pi * 2

        r = np.sqrt(pos_x ** 2 + pos_y ** 2)
        lat = np.arctan2(pos_z, r)
        e2 = F * (2 - F)
        while True:
            lat2 = lat
            c = 1/(np.sqrt(1 - e2 * (np.sin(lat2) ** 2)))
            lat = np.arctan2(pos_z + c * e2 *np.sin(lat2), r)
            if abs(lat - lat2) < 1e-10:
                break
        alt = r / np.cos(lat)- c;
        alt *= XKMPER
        return lon, lat, alt

    def find_aos(self, time, lon, lat):
        pass

    def find_aol(self, time, lon, lat):
        pass

    def get_observer_look(self, time, lon, lat, alt):
        """Calculate observers look angle to a satellite.
        http://celestrak.com/columns/v02n02/
        """
        (pos_x, pos_y, pos_z), (vel_x, vel_y, vel_z) = self.get_position(time, normalize=False)
        (opos_x, opos_y, opos_z), (ovel_x, ovel_y, ovel_z) = \
                                    astronomy.observer_position(time, lon, lat, alt)
        theta = (astronomy.gmst(time) + lon) % (2*np.pi)

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

        top_s = sin_lat*cos_theta*rx + sin_lat*sin_theta*ry - cos_lat*rz
        top_e = -sin_theta*rx + cos_theta*ry
        top_z = cos_lat*cos_theta*rx + cos_lat*sin_theta*ry + sin_lat*rz

        az = np.arctan(-top_e/top_s)
        if top_s > 0:
            az = az + np.pi
        if az < 0:
            az = az + 2*np.pi

        rg = np.sqrt(rx*rx + ry*ry + rz*rz)
        el = np.arcsin(top_z/rg)
        w = (rx*rvx + ry*rvy + rz*rvz)/rg

        return az, el, rg, w

    def get_orbit_number(self, time, tbus_style=False):
        #TODO: Handled corner cases of non AN epoch and propagation to near AN
        # and use node periode instead of revolutions 
        lon, lat, alt = self.get_lonlatalt(self.tle.epoch)
        dt = astronomy._days(time - self.tle.epoch)
        orbit = int(self.tle.orbit + self.tle.mean_motion * dt + 
                 self.tle.mean_motion_derivative * dt**2 + 
                 self.tle.mean_motion_sec_derivative * dt**3)
        if tbus_style:
            orbit += 1
        return orbit

class OrbitElements(object):

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


def _sgp4(orbit_elements, time):
    # for near earth orbits
    # see fx http://celestrak.com/

    perigee = orbit_elements.perigee
    a_0pp = orbit_elements.semi_major_axis
    e_0 = orbit_elements.excentricity
    i_0 = orbit_elements.inclination
    n_0pp = orbit_elements.original_mean_motion
    k_2 = CK2
    k_4 = CK4
    k_e = XKE
    bstar = orbit_elements.bstar
    w_0 = orbit_elements.arg_perigee
    M_0 = orbit_elements.mean_anomaly
    W_0 = orbit_elements.right_ascension
    t_0 = orbit_elements.epoch
    A30 = -XJ3 * AE**3

    if perigee < 98:
        s = 20/XKMPER + AE
        qoms2t = (QOMS2T ** 0.25 + S - s) ** 4
    elif perigee < 156:
        s = a_0pp * (1 - e_0) - S + AE 
        qoms2t = (QOMS2T ** 0.25 + S - s) ** 4
    else:
        qoms2t = QOMS2T
        s = S

    theta = np.cos(i_0)
    xi = 1 / (a_0pp - s)
    beta_0 = np.sqrt(1 - e_0 ** 2)
    eta = a_0pp * e_0 * xi

    C_2 = (qoms2t * xi**4 * n_0pp * (1 - eta**2)**(-3.5) *
           (a_0pp * (1 + 1.5 * eta**2 + 4 * e_0 * eta + e_0 * eta**3) +
            1.5 * (k_2 * xi) / (1 - eta**2) * (-0.5 + 1.5 * theta**2)*
            (8 + 24 * eta**2 + 3 * eta**4)))

    C_1 = bstar * C_2

    C_3 = (qoms2t * xi ** 5 * A30 * n_0pp * AE * np.sin(i_0) / (k_2 * e_0))

    coef = 2 * qoms2t * xi**4 * a_0pp * beta_0**2*(1-eta**2)**(-7/2.0)

    C_4 = (coef * n_0pp *
           ((2 * eta * (1 + e_0 * eta) + e_0/2.0 + (eta**3)/2.0) -
            2 * k_2 * xi / (a_0pp * (1 - eta**2)) *
            (3*(1-3*theta**2) * (1 + (3*eta**2)/2.0 - 2*e_0*eta - e_0*eta**3/2.0) +
             3/4.0*(1-theta**2)*(2*eta**2 - e_0*eta - e_0*eta**3)*np.cos(2*w_0))))
    
    C_5 = coef * (1 + 11/4.0 * eta * (eta + e_0) + e_0 * eta**3)
    D_2 = 4 * a_0pp * xi * C_1**2
    D_3 = 4/3.0 * a_0pp * xi**2 * (17*a_0pp + s) * C_1**3
    D_4 = 2/3.0 * a_0pp * xi**3 * (221*a_0pp + 31*s) * C_1**4

    # Secular effects of atmospheric drag and gravitation
    dt = astronomy._days(time - t_0) * XMNPDA

    M_df = (M_0 + (1 +
                   3*k_2*(-1 + 3*theta**2)/(2*a_0pp**2 * beta_0**3) +
                   3*k_2**2*(13 - 78*theta**2 + 137*theta**4)/
                   (16*a_0pp**4*beta_0**7))*
            n_0pp*dt)
    w_df = (w_0 + (-3*k_2*(1 - 5*theta**2)/(2*a_0pp**2*beta_0**4) +
                    3 * k_2**2 * (7 - 114*theta**2 + 395*theta**4)/
                    (16*a_0pp*beta_0**8) +
                    5*k_4*(3-36*theta**2+49*theta**4)/(4*a_0pp**4*beta_0**8))*
            n_0pp*dt)
    W_df = (W_0 + (-3*k_2*theta/(a_0pp**2*beta_0**4) +
                    3*k_2**2*(4*theta- 19*theta**3)/(2*a_0pp**4*beta_0**8) +
                    5*k_4*theta*(3-7*theta**2)/(2*a_0pp**4*beta_0**8))*
            n_0pp*dt)
    deltaw = bstar * C_3 * np.cos(w_0)*dt
    deltaM = (-2/3.0 * qoms2t * bstar * xi**4 * AE / (e_0*eta) *
               ((1 + eta * np.cos(M_df))**3 - (1 + eta * np.cos(M_0))**3))
    M_p = M_df + deltaw + deltaM
    w = w_df - deltaw - deltaM
    W = (W_df - 21/2.0 * (n_0pp * k_2 * theta)/(a_0pp**2 * beta_0**2) *
         C_1 * dt**2)

    e = (e_0 -
         bstar * C_4 * dt -
         bstar * C_5 * (np.sin(M_p) - np.sin(M_0)))

    a = a_0pp * (1 - C_1 * dt - D_2 * dt**2 - D_3 * dt**3 - D_4 * dt**4)**2
    L = M_p + w + W + n_0pp * (3/2.0 * C_1 * dt**2 +
                               (D_2 + 2 * C_1 ** 2) * dt**3 +
                               1/4.0 * (3*D_3 + 12*C_1*D_2 + 10*C_1**3)*dt**4 +
                               1.0/5 * (3*D_4 + 12*C_1*D_3 + 6*D_2**2 +
                                        30*C_1**2*D_2 + 15*C_1**4)*dt**5)
    beta = np.sqrt(1 - e**2)
    n = k_e / (a ** (3/2.0))
    
    # Long-period periodic terms
    a_xN = e * np.cos(w)
    a_yNL = A30 * np.sin(i_0) / (4.0 * k_2 * a * beta**2)
    L_L = a_yNL/2 * a_xN * ((3 + 5 * theta) / (1 + theta))
    L_T = L + L_L
    a_yN = e * np.sin(w) + a_yNL
    
    U = (L_T - W) % (np.pi * 2)
    
    Epw = U
    for i in range(10):
        DeltaEpw = ((U - a_yN * np.cos(Epw) + a_xN  * np.sin(Epw) - Epw) /
                    (-a_yN * np.sin(Epw) - a_xN * np.cos(Epw) + 1))
        Epw = Epw + DeltaEpw
        if DeltaEpw < 10e-12:
            break

    # preliminary quantities for short-period periodics
    
    ecosE = a_xN * np.cos(Epw) + a_yN * np.sin(Epw)
    esinE = a_xN * np.sin(Epw) - a_yN * np.cos(Epw)

    e_L = (a_xN**2 + a_yN**2)**(0.5)
    p_L = a * (1 - e_L**2)
    r = a * (1 - ecosE)
    rdot = k_e * np.sqrt(a)/r * esinE
    rfdot = k_e * np.sqrt(p_L) / r
    cosu = a / r * (np.cos(Epw) - a_xN +
                    (a_yN * (esinE) / (1 + np.sqrt(1 - e_L**2))))
    sinu = a / r * (np.sin(Epw) - a_yN +
                    (a_xN * (esinE) / (1 + np.sqrt(1 - e_L**2))))
    u = np.arctan2(sinu, cosu)


    cos2u = np.cos(2*u)
    sin2u = np.sin(2*u)

    Deltar = k_2/(2*p_L) * (1 - theta**2) * cos2u
    Deltau = -k_2/(4*p_L**2) * (7*theta**2 - 1) * sin2u
    DeltaW = 3*k_2 * theta / (2 * p_L**2) * sin2u
    Deltai = 3*k_2 * theta / (2 * p_L**2) * cos2u * np.sin(i_0)
    Deltardot = - k_2 * n / p_L * (1 - theta**2) * sin2u
    Deltarfdot = k_2 * n / p_L * ((1 - theta**2) * cos2u -
                                  3/2.0 * (1 - 3*theta**2))

    # osculating quantities

    r_k = r * (1 - 3/2.0 * k_2 * np.sqrt(1 - e_L**2)/p_L**2 *
               (3 * theta**2 - 1)) + Deltar
    u_k = u + Deltau
    W_k = W + DeltaW
    i_k = i_0 + Deltai
    rdot_k = rdot + Deltardot
    rfdot_k = rfdot + Deltarfdot

    M_x = -np.sin(W_k) * np.cos(i_k)
    M_y = np.cos(W_k) * np.cos(i_k)
    M_z = np.sin(i_k)

    N_x = np.cos(W_k)
    N_y = np.sin(W_k)
    N_z = 0
    
    U_x = M_x * np.sin(u_k) + N_x * np.cos(u_k)
    U_y = M_y * np.sin(u_k) + N_y * np.cos(u_k)
    U_z = M_z * np.sin(u_k) + N_z * np.cos(u_k)
    
    V_x = M_x * np.cos(u_k) - N_x * np.sin(u_k)
    V_y = M_y * np.cos(u_k) - N_y * np.sin(u_k)
    V_z = M_z * np.cos(u_k) - N_z * np.sin(u_k)


    r_x = r_k * U_x
    r_y = r_k * U_y
    r_z = r_k * U_z
    
    rdot_x = rdot_k * U_x + rfdot_k * V_x
    rdot_y = rdot_k * U_y + rfdot_k * V_y
    rdot_z = rdot_k * U_z + rfdot_k * V_z

    return (r_x, r_y, r_z), (rdot_x, rdot_y, rdot_z)


if __name__ == "__main__":
    obs_lon, obs_lat = np.deg2rad((12.4143, 55.9065))
    obs_alt = 0.02
    o = Orbital(satellite="noaa 18") #, tle_file="/net/prodsat/datadb/sat/orbit/tle/tle_20110327.txt")
    print o

    t_start = datetime.datetime(2011, 3, 28, 2, 15)
    t_stop = t_start + datetime.timedelta(minutes=17)
    t = t_start
    while t < t_stop:
        t += datetime.timedelta(seconds=15)
        lon, lat, alt = o.get_lonlatalt(t)
        lon, lat = np.rad2deg((lon, lat))
        az, el, rg, w = o.get_observer_look(t, obs_lon, obs_lat, obs_alt)
        az, el = np.rad2deg((az, el))
        print str(t) + ': ', "%6.2f, %6.2f, %6.2f - %6.2f, %5.2f"%(lon, lat, alt, az, el), "%7.2f"%rg, "%7.4f"%w
