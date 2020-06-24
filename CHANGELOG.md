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
