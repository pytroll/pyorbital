    
Computing satellite postion
---------------------------
The orbital module enables computation of satellite position and velocity at a
specific time:

    >>> from pyorbital.orbital import Orbital
    >>> from datetime import datetime
    >>> orb = Orbital("noaa 18")
    >>> now = datetime.utcnow()
    >>> # Get normalized position and velocity of the satellite:
    >>> orb.get_position(now)
    ([0.57529384846822862, 0.77384005228105424, 0.59301408257897559],
    [0.031846489698768146, 0.021287993461926374, -0.05854106186659274])
    >>> # Get longitude, latitude and altitude of the satellite:
    >>> orb.get_lonlatalt(now)
    (-1.1625895579622014, 0.55402132517640568, 847.89381184656702)
