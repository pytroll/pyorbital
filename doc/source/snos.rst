Finding Simultaneous Nadir Overpasses for two different platforms
------------------------------------------------------------------

Pyorbital facilitates the identification and calculation of simultaneous nadir
overpasses (SNOs) for two different platforms. This is for example useful when
you want to colocate observations in time and space and where the viewing
conditions are the same (nadir viewing) for both platforms. One example is when
validating clous parameters derived from Imager sensors like VIIRS, AVHRR or
MODIS against the active nadir viewing sensors of the Calipso and CloudSat in
the A-train constellation.

**Example usage:**

Below is an example searching for SNO points between the EOS-Aqua and the
Metop-B platforms. We allow a maximum of 10 minutes between the two
observations where the two nadir viewsare crossing. Thus here *simultaneous*
means within 10 minutes. In the example we limit the time period to January
2015 and we provide a directory with daily tle-files for the entire month.

.. code-block:: bash
                
  %> python ./bin/get_snos.py -s 20150101 -e 20150131 -t 10 -p Metop-B -c ./examples/snos.yaml

A few basic parameters need to be defined in the configuration yaml file **snos.yaml**:

.. code-block:: bash

  station:
    longitude: 16.1465
    latitude: 58.5780
    altitude: 0.03

  tle-dirs:
    - /home/a000680/data/tles/2015/201501

  platform_names:
    - NOAA-18
    - NOAA-19
    - NOAA-20
    - Suomi-NPP
    - Metop-B
    - Metop-C

  tle-file-format: 'tle-%Y%m%d.txt'

  reference_platform:
    name: EOS-Aqua

  
