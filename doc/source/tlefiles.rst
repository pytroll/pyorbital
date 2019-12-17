TLE files
---------
Pyorbital has a module for parsing NORAD TLE-files

    >>> from pyorbital import tlefile
    >>> tle = tlefile.read('noaa 18', '/path/to/my/tle_file.txt')
    >>> tle.inclination
    99.043499999999995

If no path is given pyorbital tries to read the earth observation TLE-files from celestrak.com

