Changelog
=========

v1.1.0 (2016-10-27)
-------------------

- Update changelog. [Martin Raspaud]

- Bump version: 1.0.1 → 1.1.0. [Martin Raspaud]

- Merge branch 'master' into develop. [Martin Raspaud]

- Enable travis testing for py3. [Antonio Valentino]

- Fix regression in TLE reading. [Antonio Valentino]

- Python 3 compatibility. [Antonio Valentino]

v1.0.1 (2016-02-17)
-------------------

- Update changelog. [Martin Raspaud]

- Bump version: 1.0.0 → 1.0.1. [Martin Raspaud]

- Change sun_angle test to AlmostEqual. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


v1.0.0 (2015-08-25)
-------------------

- Update changelog. [Martin Raspaud]

- Bump version: 0.3.2 → 1.0.0. [Martin Raspaud]

- Cleanup. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Fix version number. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Cosmetics. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge pull request #2 from pnuu/feature_tle_lookup. [Martin Raspaud]

  Use NORAD catalog numbers for TLE reading

- Example file for mapping OSCAR platform names and NORAD catalog
  numbers. [Panu Lahtinen]

- Add setup.cfg for easy rpm generation. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge branch 'develop' of github.com:mraspaud/pyorbital into develop.
  [Martin Raspaud]

- Merge pull request #1 from spareeth/develop. [Martin Raspaud]

  changes to avhrr_gacfunction and read_tle_decimal

- Added '+' as a condition in the read_tle function. [Sajid Pareeth]

- Renaming the variable scans_nb to scan_times in offset in avhrr_gac
  function. [Sajid Pareeth]

- Bugfix: eccentricity too low message formatting. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Allow reading TLE from the most recent file described by the TLES env.
  [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Change decimate to frequency in avhrr instruments. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Add the avhrr instrument, gac version. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Accept missing zeros in TLE (old noaa compatibility). [Martin Raspaud]

- Add the horizon parameter to get_next_passes to get the
  risetime/falltime at given angle. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge branch 'master' into develop. [Martin Raspaud]

- Fix backwards numpy compatibility. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


v0.3.2 (2014-04-10)
-------------------

- Merge branch 'develop' [Martin Raspaud]

- Bump up version number. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge branch 'feature-no-scipy' into develop. [Martin Raspaud]

- Remove scipy dependencies. [Martin Raspaud]

  Was depending on scipy.optimize, brent and brentq function.
  Replaced by secant method root finding and successive parabolic
  interpolation local minimum finding.

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Correcting the travis file. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


v0.3.1 (2014-02-24)
-------------------

- Bugfix in travis file. [Martin Raspaud]

- Bump up version number. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Fixed documentation. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Cleanup. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- New nadir computations for geoloc. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- More unit tests. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


v0.3.0 (2014-01-07)
-------------------

- Auto update version number in documentation. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Change to version file and bump up to v0.3.0. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Cleanup the testfiles. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Add a test to read tle from file. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Fix doc path in MANIFEST.in. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


v0.2.4 (2014-01-07)
-------------------

- Merge branch 'feature-travis' into pre-master. [Martin Raspaud]

- Add test for tle reading, cleanup and make ready for travis. [Martin
  Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Cleanup. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Add function to fetch the tle files from internet manually. [Martin
  Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Adding the viirs instrument. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Change sphinx theme. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Fix doc for readthedocs. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Remove unused old file. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge branch 'geoloc' into pre-master. [Martin Raspaud]

- Work on geolocation. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Numpyze the orbital computation. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Add some logging in tle file fetching. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Fix syntax error in doc/conf.py. [Martin Raspaud]

- Make the scan angle of avhrr an argument. [Martin Raspaud]

- Factorize avhrr code (geoloc definition) [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Add Mikhail's definition of AMSU-A. [Martin Raspaud]

- Add instrument examples for geoloc. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge branch 'geoloc' of github.com:mraspaud/pyorbital into geoloc.
  [Martin Raspaud]

- Try fixing nadir. [Martin Raspaud]

- Fix attitude. [Martin Raspaud]

- Updated doc and copyright. [Martin Raspaud]

- Add geoloc example. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Merge branch 'feature-vectorize' into geoloc. [Martin Raspaud]

- Vectorize the days function. [Martin Raspaud]

- Merge branch 'master' into geoloc. [Martin Raspaud]

- Merge branch 'pre-master' into geoloc. [Martin Raspaud]

- Cosmetics. [Martin Raspaud]

- Computations for true nadir. [Martin Raspaud]

- Bugfix in the example and added attitude correction (roll and pitch
  for now). [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Cosmetic, be consistent in name og time argument as 'utc_time' [Lars
  Orum Rasmussen]

- Get_zenith_overpass replaced by Martin's get_next_passes. [Lars Orum
  Rasmussen]

- Add sun_earth_distance_correction function. [Martin Raspaud]

v0.2.3 (2013-03-07)
-------------------

- Merge branch 'release-0.2.3' [Martin Raspaud]

- Merge branch 'pre-master' into release-0.2.3. [Martin Raspaud]

- Bumped up version number. [Martin Raspaud]

- Corrected search for previous an_time with a substracted 10 min. dt.
  [Esben S. Nielsen]

- Merge branch 'release-0.2.2' [Martin Raspaud]

- Import with_statement in test_aiaa.py for python 2.5 compliance.
  [Esben S. Nielsen]

- Made unit tests python 2.5 and 2.6 compliant. [Esben S. Nielsen]

- Removed download URL from setup.py. [Esben S. Nielsen]

- Bumped version number and marked as stable. [Esben S. Nielsen]

- Better handling of time deltas in test_aiaa.py. [Esben S. Nielsen]

- Updated equator test with position check. [Esben S. Nielsen]

- Now uses nodal period for orbit number calculation instead of revs/day
  for mean motion. [Esben S. Nielsen]

- Orbit number now handles epoch AN mis-match. Made AIAA unit test path
  agnostic. [Esben S. Nielsen]

- Better __main__ [Lars Orum Rasmussen]

- Adding risetime and falltime functions, and improving the
  get_zenith_overpass function. [Adam Dybbroe]

- Editorial. [Adam Dybbroe]

- Cleanup. [Martin Raspaud]

  Signed-off-by: Martin Raspaud <martin.raspaud@smhi.se>


- Feature: Correcting/adding test cases from the aiaa. [Martin Raspaud]

- Style: raises NotImplementedErrors instead of just Exceptions. [Martin
  Raspaud]

- Merge branch 'pre-master' of github.com:mraspaud/pyorbital into pre-
  master. [Martin Raspaud]

- Adding new function get_zenith_overpass to get the time when the
  satellite passes over zenith relative to an observer on ground. [Adam
  Dybbroe]

- Feature: Added checksum for tle lines. [Martin Raspaud]

v0.2.1 (2012-06-01)
-------------------

- Updated version number. [Martin Raspaud]

- Added pyorbital path to doc/source/conf.py. [Esben S. Nielsen]

- Updated docs and added license and manifest. [Esben S. Nielsen]

- Merge branch 'pre-master' of https://github.com/mraspaud/pyorbital
  into pre-master. [Adam Dybbroe]

- Merge branch 'pre-master' of https://github.com/mraspaud/pyorbital
  into pre-master. [Lars Orum Rasmussen]

- Added access to line1 and line2 in a Tle instance. [Lars Orum
  Rasmussen]

  Change satellite to platform


- Spelling error. [Adam Dybbroe]

v0.2.0 (2012-05-14)
-------------------

- Prepared for pypi. [Martin Raspaud]

- Merge branch 'geoloc' into pre-master. [Martin Raspaud]

- Added now compute pixels on the ellipsoid, not on the sphere anymore.
  [Martin Raspaud]

- Merge branch 'master' into geoloc. [Martin Raspaud]

- Updated the geoloc todo list. [Martin Raspaud]

- Added the geoloc module. [Martin Raspaud]

- Merge branch 'master' into pre-master. [Martin Raspaud]

  Conflicts:
  	pyorbital/tlefile.py


- Corrected handling of mean motion and orbitnumber fields in
  tlefiles.py. [Esben S. Nielsen]

- Testing getting the orbit number from the TLEs. [Adam.Dybbroe]

- Fixing bug in tle file reading, so that also NPP and other satellites
  with orbit numbers less than 9999 can be handled. [Adam.Dybbroe]

- Typo. [Adam.Dybbroe]

- Merge branch 'master' into pre-master. [Martin Raspaud]

- Removed html submodule. [Martin Raspaud]

- Fixing bug in function sun_zenith_angle. Changing interfaces so that
  all public functions expects lon,lat in degrees. All internal
  functions us radians. Made the lsmt and local_hour_angle functions
  private. [Adam.Dybbroe]

- Adding main. [Adam.Dybbroe]

- Gathering unit tests to the tests-directory. [Adam.Dybbroe]

- Added separate test-script for astronomy.py. [Adam.Dybbroe]

- Collected all unit test scripts under the tests directory.
  [Adam.Dybbroe]

- Merge branch 'release-0.2.0' [Martin Raspaud]

  Conflicts:
  	doc/build
  	setup.py


- Bumped version number to 0.2.0. [Martin Raspaud]

- Added html documentation. [Martin Raspaud]

- Corrected sgp4's propagate in the case of array as input, and cleaned
  up. [Martin Raspaud]

- Fixed calling test_aiaa from another directory. [Martin Raspaud]

- Vectorize merge. [Martin Raspaud]

- Merging master branch. [Martin Raspaud]

- Remove html submodule. [Martin Raspaud]

- Remove html submodule. [Martin Raspaud]

- Added Esben in the author field. [Martin Raspaud]

- Removed unneded .pyc file. [Martin Raspaud]

- Added unittests. [Esben S. Nielsen]

- Corrected observer_look function and added first unittest. [Esben S.
  Nielsen]

- Corrected observer_pos in astronomy. [Esben S. Nielsen]

- Setting up documentation. [Martin Raspaud]

v0.1.0 (2011-10-03)
-------------------

- Merge branch 'release-0.1.0' [Martin Raspaud]

- Bumped version number to 0.1.0. [Martin Raspaud]

- Merge branch 'dundee_port' into pre-master. [Martin Raspaud]

- Cleanup and documentation. [Martin Raspaud]

- Now using unittest module for aiaa test cases. [Martin Raspaud]

- Added licences, and removed prints. [Martin Raspaud]

- Added basic tests to pyorbital. [Martin Raspaud]

- Ported SGP4 code to Dundee implementation. [Esben S. Nielsen]

- Ported sgp4 init. [Esben S. Nielsen]

- Added the first unit test :) [Martin Raspaud]

- New gmst function (from AIAA paper). Cleaning. [Martin Raspaud]

- Merged DMI and SMHI versions. [Esben S. Nielsen]

- Made the package more package-like. [Martin Raspaud]

- Cleanup of astronomy file. [Martin Raspaud]

- Added a readme file. [Martin Raspaud]

- Added astronomy.py file. [Martin Raspaud]


