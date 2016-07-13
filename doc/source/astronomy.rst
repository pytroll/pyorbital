Computing astronomical parameters
---------------------------------
The astronomy module enables computation of certain parameters of interest for
satellite remote sensing for instance the Sun-zenith angle:

    >>> from pyorbital import astronomy
    >>> from datetime import datetime
    >>> utc_time = datetime(2012, 5, 15, 15, 45)
    >>> lon, lat = 12, 56
    >>> astronomy.sun_zenith_angle(utc_time, lon, lat)
    62.685986438071602

