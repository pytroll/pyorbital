#
#

import datetime

tle_urls = ('http://celestrak.com/NORAD/elements/weather.txt',
            'http://celestrak.com/NORAD/elements/resource.txt')

def read(satellite, tle_file=None):
    return Tle(satellite, tle_file)

class Tle(object):

    def __init__(self, satellite, tle_file=None):
        satellite = satellite.strip().upper()

        if tle_file:
            urls = (tle_file,)
            open_func = open
        else:
            import urllib2
            urls = tle_urls
            open_func = urllib2.urlopen
        
        tle = ""
        for url in urls:
            fp = open_func(url)
            for l0 in fp:
                l1, l2 = fp.next(), fp.next()
                if l0.strip() == satellite:
                    tle = l1.strip() + "\n" + l2.strip()
                    break
            fp.close()
            if tle:
                break
        
        if not tle:
            raise AttributeError, "Found no TLE entry for '%s'" % satellite

        self._read_tle(tle)

    def _read_tle(self, lines):

        def _read_tle_decimal(rep):
            num = int(rep[:-2]) * 1.0e-5
            exp = int(rep[-2:])
            return num * 10 ** exp

        lines = lines.split()
        self.satnumber = lines[1][:5]
        self.classification = lines[1][5:]
        self.id_launch_year = lines[2][:2]
        self.id_launch_number = lines[2][2:5]
        self.id_launch_piece = lines[2][5:]
        self.epoch_year = int(lines[3][:2])
        self.epoch_day = float(lines[3][2:])
        self.epoch = (datetime.datetime.strptime(lines[3][:2], "%y") +
                             datetime.timedelta(days=float(lines[3][2:]) - 1))
        self.mean_motion_derivative = float(lines[4])
        self.mean_motion_sec_derivative = _read_tle_decimal(lines[5])
        self.bstar = _read_tle_decimal(lines[6])
        self.ephemeris_type = int(lines[7])
        self.element_number = int(lines[8][:-1])

        self.inclination = float(lines[11])
        self.right_ascension = float(lines[12])
        self.excentricity = int(lines[13]) * 10 ** -7
        self.arg_perigee = float(lines[14])
        self.mean_anomaly = float(lines[15])
        self.mean_motion = float(lines[16][:11])
        self.orbit = int(lines[16][11:-1])

    def __str__(self):
        import pprint, StringIO
        s = StringIO.StringIO()
        pprint.pprint(self.__dict__, s)
        return s.getvalue()[:-1]

if __name__ == '__main__':
    tle = read('noaa 19')
    print tle
