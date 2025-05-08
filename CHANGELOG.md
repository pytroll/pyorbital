## Version 1.10.1 (2025/05/08)


### Pull Requests Merged

#### Features added

* [PR 189](https://github.com/pytroll/pyorbital/pull/189) - Work with gcps for avhrr

In this release 1 pull request was closed.


## Version 1.10.0 (2025/04/10)

### Issues Closed

* [Issue 184](https://github.com/pytroll/pyorbital/issues/184) - avhrr_gac instrument definition is broken ([PR 185](https://github.com/pytroll/pyorbital/pull/185) by [@oembury](https://github.com/oembury))

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 185](https://github.com/pytroll/pyorbital/pull/185) - Fix avhrr gac geometry ([184](https://github.com/pytroll/pyorbital/issues/184))

#### Features added

* [PR 188](https://github.com/pytroll/pyorbital/pull/188) - Cleanup for #185
* [PR 186](https://github.com/pytroll/pyorbital/pull/186) - Support yaw arrays for ScanGeometry with same broadcasting as roll / pitch
* [PR 177](https://github.com/pytroll/pyorbital/pull/177) - Fix datetime imports

In this release 4 pull requests were closed.


## Version 1.9.2 (2024/12/12)

### Pull Requests Merged

#### Bugs fixed

* [PR 175](https://github.com/pytroll/pyorbital/pull/175) - Remove pytz dependency

In this release 1 pull request was closed.


## Version 1.9.1 (2024/11/29)

### Pull Requests Merged

#### Bugs fixed

* [PR 172](https://github.com/pytroll/pyorbital/pull/172) - Update CI build action to work without setup.py

#### Features added

* [PR 170](https://github.com/pytroll/pyorbital/pull/170) - Deprecate PPP_CONFIG_DIR for specifying config path

In this release 2 pull requests were closed.

## Version 1.9.0 (2024/11/29)

### Issues Closed

* [Issue 160](https://github.com/pytroll/pyorbital/issues/160) - SENTINEL-2 TLE files ([PR 161](https://github.com/pytroll/pyorbital/pull/161) by [@simonrp84](https://github.com/simonrp84))
* [Issue 119](https://github.com/pytroll/pyorbital/issues/119) - tlefile.py seem to use satellite name rather than the international designator number ([PR 118](https://github.com/pytroll/pyorbital/pull/118) by [@adybbroe](https://github.com/adybbroe))
* [Issue 80](https://github.com/pytroll/pyorbital/issues/80) - NotImplementedError: Mode "Near-space, simplified equations" not implemented ([PR 124](https://github.com/pytroll/pyorbital/pull/124) by [@JonathanMaes](https://github.com/JonathanMaes))

In this release 3 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 159](https://github.com/pytroll/pyorbital/pull/159) - Bugfix use of datetime.datetime objects in call to `get_last_an_time`

#### Features added

* [PR 171](https://github.com/pytroll/pyorbital/pull/171) - Refactor _SGDB4 class
* [PR 169](https://github.com/pytroll/pyorbital/pull/169) - Switch to use pyproject.toml instead of setup.py, and skip versioneer
* [PR 161](https://github.com/pytroll/pyorbital/pull/161) - Add new platforms and clarify help message for `check_platform` ([160](https://github.com/pytroll/pyorbital/issues/160))
* [PR 146](https://github.com/pytroll/pyorbital/pull/146) - Update CI to use Python 3.10 - 3.12 and plain Miniforge
* [PR 124](https://github.com/pytroll/pyorbital/pull/124) - Implement SGDP4_NEAR_SIMP propagation ([80](https://github.com/pytroll/pyorbital/issues/80))
* [PR 118](https://github.com/pytroll/pyorbital/pull/118) - Refactor and improve tests concerning the TLE file handling ([119](https://github.com/pytroll/pyorbital/issues/119))
* [PR 98](https://github.com/pytroll/pyorbital/pull/98) - Add slstr ([81](https://github.com/pytroll/pyorbital/issues/81))

In this release 8 pull requests were closed.

## Version 1.8.3 (2024/06/25)

### Issues Closed

* [Issue 151](https://github.com/pytroll/pyorbital/issues/151) - Issue Calculating Accurate View Zenith Angles on Terra Satellite Overpasses

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 156](https://github.com/pytroll/pyorbital/pull/156) - Fix dtype preservation in astronomy functions

In this release 1 pull request was closed.


## Version 1.8.2 (2024/02/05)

### Issues Closed

* [Issue 140](https://github.com/pytroll/pyorbital/issues/140) - pyorbital cannot read TLE for MTG-I1 / Meteosat-12 ([PR 141](https://github.com/pytroll/pyorbital/pull/141) by [@gerritholl](https://github.com/gerritholl))
* [Issue 139](https://github.com/pytroll/pyorbital/issues/139) - `Orbital` cannot get TLEs from the internet

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 143](https://github.com/pytroll/pyorbital/pull/143) - Fix a bug in using TLES env variable

#### Features added

* [PR 141](https://github.com/pytroll/pyorbital/pull/141) - Add Meteosat-12 to platforms.txt ([140](https://github.com/pytroll/pyorbital/issues/140))

In this release 2 pull requests were closed.


## Version 1.8.1 (2024/01/05)

### Pull Requests Merged

#### Bugs fixed

* [PR 138](https://github.com/pytroll/pyorbital/pull/138) - Update celestrak urls ([139](https://github.com/pytroll/pyorbital/issues/139))

#### Features added

* [PR 137](https://github.com/pytroll/pyorbital/pull/137) - Prettify the RTD pages

#### Documentation changes

* [PR 137](https://github.com/pytroll/pyorbital/pull/137) - Prettify the RTD pages
* [PR 132](https://github.com/pytroll/pyorbital/pull/132) - Add .readthedocs.yaml

In this release 4 pull requests were closed.


## Version 1.8.0 (2023/07/12)

### Issues Closed

* [Issue 112](https://github.com/pytroll/pyorbital/issues/112) - Is the TLES environment variable described? ([PR 113](https://github.com/pytroll/pyorbital/pull/113) by [@adybbroe](https://github.com/adybbroe))

In this release 1 issue was closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 129](https://github.com/pytroll/pyorbital/pull/129) - Fix bug getting local tlefiles
* [PR 128](https://github.com/pytroll/pyorbital/pull/128) - Fix typo in VIIRS geoloc definition
* [PR 121](https://github.com/pytroll/pyorbital/pull/121) - fixed geoloc_example and added variable descriptions

#### Features added

* [PR 120](https://github.com/pytroll/pyorbital/pull/120) - Update versioneer to stop using deprecated distutils module.
* [PR 113](https://github.com/pytroll/pyorbital/pull/113) - Make use of env variables free from satpy ([112](https://github.com/pytroll/pyorbital/issues/112))

#### Documentation changes

* [PR 113](https://github.com/pytroll/pyorbital/pull/113) - Make use of env variables free from satpy ([112](https://github.com/pytroll/pyorbital/issues/112))

In this release 6 pull requests were closed.


## Version 1.7.3 (2022/07/11)

### Pull Requests Merged

#### Bugs fixed

* [PR 103](https://github.com/pytroll/pyorbital/pull/103) - www.celestrak.org → celestrak.org ([101](https://github.com/pytroll/pyorbital/issues/101))

#### Documentation changes

* [PR 104](https://github.com/pytroll/pyorbital/pull/104) - Fixed comment in sun_earth_distance_correction

In this release 2 pull requests was closed.


## Version 1.7.2 (2022/07/07)

### Issues Closed

* [Issue 100](https://github.com/pytroll/pyorbital/issues/100) - TLE URIs land at redirect page ([PR 101](https://github.com/pytroll/pyorbital/pull/101) by [@gerritholl](https://github.com/gerritholl))
* [Issue 97](https://github.com/pytroll/pyorbital/issues/97) - sun_earth_distance_correction returns the square of the sun-earth distance relativ to 1 AU ([PR 99](https://github.com/pytroll/pyorbital/pull/99) by [@depion](https://github.com/depion))

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 101](https://github.com/pytroll/pyorbital/pull/101) - Update celestrak URIs ([100](https://github.com/pytroll/pyorbital/issues/100))
* [PR 99](https://github.com/pytroll/pyorbital/pull/99) - Fixed earth sun distance ([97](https://github.com/pytroll/pyorbital/issues/97))

In this release 2 pull requests were closed.


## Version 1.7.1 (2021/12/22)

### Pull Requests Merged

#### Bugs fixed

* [PR 92](https://github.com/pytroll/pyorbital/pull/92) - Fix bogus designator assignment

In this release 1 pull request was closed.


## Version 1.7.0 (2021/12/20)

### Issues Closed

* [Issue 90](https://github.com/pytroll/pyorbital/issues/90) - get_observer_look raises IndexError on numpy 1.21.4
* [Issue 85](https://github.com/pytroll/pyorbital/issues/85) - Azimuth/elevation output not changing with 1-second increments
* [Issue 79](https://github.com/pytroll/pyorbital/issues/79) - ModuleNotFoundError: 'pyorbital' is not a package
* [Issue 72](https://github.com/pytroll/pyorbital/issues/72) - Unexpected Nans in get_observer_look_no_tle
* [Issue 38](https://github.com/pytroll/pyorbital/issues/38) - Issue with pyorbital.planets

In this release 5 issues were closed.

### Pull Requests Merged

#### Features added

* [PR 91](https://github.com/pytroll/pyorbital/pull/91) - Add get_observer_look test for scalar case and update stickler config ([91](https://github.com/pytroll/pyorbital/issues/91))
* [PR 89](https://github.com/pytroll/pyorbital/pull/89) - Change tested Python versions to 3.8, 3.9 and 3.10
* [PR 83](https://github.com/pytroll/pyorbital/pull/83) - Add Sentinel-5P to platform_names
* [PR 78](https://github.com/pytroll/pyorbital/pull/78) - Add parser to read TLEs from Multi Mission Administrative Messages

#### Documentation changes

* [PR 82](https://github.com/pytroll/pyorbital/pull/82) - Add historical TLE files link

In this release 5 pull requests were closed.


## Version 1.6.1 (2021/04/12)

### Issues Closed

* [Issue 63](https://github.com/pytroll/pyorbital/issues/63) - Runtime error in get_next_passes ([PR 64](https://github.com/pytroll/pyorbital/pull/64))
* [Issue 62](https://github.com/pytroll/pyorbital/issues/62) -  can  this tool run
* [Issue 22](https://github.com/pytroll/pyorbital/issues/22) - get_next_passes returns max-elevation-time time not between rise & fall time ([PR 76](https://github.com/pytroll/pyorbital/pull/76))

In this release 3 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 76](https://github.com/pytroll/pyorbital/pull/76) - Fix parabolic estimation ([22](https://github.com/pytroll/pyorbital/issues/22))
* [PR 65](https://github.com/pytroll/pyorbital/pull/65) - Add requests to the requirements
* [PR 64](https://github.com/pytroll/pyorbital/pull/64) - Fix inappropriate runtime warning ([63](https://github.com/pytroll/pyorbital/issues/63))
* [PR 60](https://github.com/pytroll/pyorbital/pull/60) - Skip tests if data are not available
* [PR 59](https://github.com/pytroll/pyorbital/pull/59) - Fix tests on i386

#### Features added

* [PR 75](https://github.com/pytroll/pyorbital/pull/75) - Fix numpy deprecation warnings
* [PR 67](https://github.com/pytroll/pyorbital/pull/67) - Added CALIPSO among platforms

In this release 7 pull requests were closed.


## Version 1.6.0 (2020/06/24)

### Issues Closed

* [Issue 52](https://github.com/pytroll/pyorbital/issues/52) - Pyorbital fails to find TLE For some satellites ([PR 53](https://github.com/pytroll/pyorbital/pull/53))
* [Issue 28](https://github.com/pytroll/pyorbital/issues/28) - Unknown units in return value for get_alt_az ([PR 46](https://github.com/pytroll/pyorbital/pull/46))

In this release 2 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 46](https://github.com/pytroll/pyorbital/pull/46) - Fix doc about get_alt_az() return units ([28](https://github.com/pytroll/pyorbital/issues/28))

#### Features added

* [PR 56](https://github.com/pytroll/pyorbital/pull/56) - Add a script to download TLEs and store them to a database
* [PR 53](https://github.com/pytroll/pyorbital/pull/53) - Added active.txt tle path to TLE_URLS ([52](https://github.com/pytroll/pyorbital/issues/52))
* [PR 50](https://github.com/pytroll/pyorbital/pull/50) - docstring fixes
* [PR 49](https://github.com/pytroll/pyorbital/pull/49) - Equatorial Crossing Time
* [PR 47](https://github.com/pytroll/pyorbital/pull/47) - Add support for MWHS-2 (FY-3) and skip edge-functions
* [PR 45](https://github.com/pytroll/pyorbital/pull/45) - Adds engineering.txt TLE source ([15](https://github.com/pytroll/pyorbital/issues/15))

#### Documentation changes

* [PR 50](https://github.com/pytroll/pyorbital/pull/50) - docstring fixes

In this release 8 pull requests were closed.


## Version 1.5.0 (2018/11/16)

### Pull Requests Merged

#### Features added

* [PR 40](https://github.com/pytroll/pyorbital/pull/40) - Add platforms.txt to package data

In this release 1 pull request was closed.

## Version 1.4.0 (2018/10/23)

### Issues Closed

* [Issue 36](https://github.com/pytroll/pyorbital/issues/36) - Issue(s) with get_next_passes
* [Issue 34](https://github.com/pytroll/pyorbital/issues/34) - Get root secant converging to wrong solution ([PR 35](https://github.com/pytroll/pyorbital/pull/35))
* [Issue 30](https://github.com/pytroll/pyorbital/issues/30) - get_observer_look turns xarray.DataArray objects into dask.array objects
* [Issue 29](https://github.com/pytroll/pyorbital/issues/29) - URL error
* [Issue 27](https://github.com/pytroll/pyorbital/issues/27) - satellite.get_lonlatalt(now) returns wrong longitude with numpy < 1.11
* [Issue 18](https://github.com/pytroll/pyorbital/issues/18) - Sun-satellite angle ranges are not consistent

In this release 6 issues were closed.

### Pull Requests Merged

#### Bugs fixed

* [PR 39](https://github.com/pytroll/pyorbital/pull/39) - Bugfix python3
* [PR 35](https://github.com/pytroll/pyorbital/pull/35) - Use Scipy brentq method instead of secant method to perform root-finding ([34](https://github.com/pytroll/pyorbital/issues/34))

#### Features added

* [PR 37](https://github.com/pytroll/pyorbital/pull/37) - Switch to versioneer
* [PR 33](https://github.com/pytroll/pyorbital/pull/33) - Remove Develop branch

In this release 4 pull requests were closed.
