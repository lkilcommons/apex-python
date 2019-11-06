# ApexPython

## What is this?
This code is a python wrapper around Art Richmond's original Modified Magnetic Apex magnetic coordinate system, as implemented in legacy Fortran. This is the same Fortran code that is available from the [CEDAR Database](https://cedarweb.vsp.ucar.edu/wiki/index.php/Tools_and_Models:Empirical_Models).

## Who should use this?
If you don't have a specific reason to use this Apex coordinates wrapper (i.e. it's a dependancy of something else you'd like to install), I'd suggest you **USE (Apexpy)[http://github.com/aburrell/apexpy] INSTEAD**. This is an older implementation of Apex in Fortran code, the ApexPython wrapper is older and less maintained, and Apexpy can do almost everything this library can do.

## Contents of this package

I've written a Python class called apex_converter, which expands upon the basic functionality in the Fortran library in several ways:

1. Vectorization - all of the functions in the apex_coverter class are equipped to handle numpy arrays
2. Result Caching - after a large call to the fortran library has completed, the results for ALL variables, not just the ones the user requested, are stored in a dictionary (called last_run).


## Installation
Tested using the Anaconda Python distribution on Ubuntu 14.04.
Also have tested using OSX system python in a Python virualenv.
If you test in windows and find problems, file an issue in the issue tracker and I'll work on it.

Currently only Python 2.7 is supported. Working on Python 3 support.

Dependanices:
- numpy
- gfortran (for f2py wrapper)

Optional Libraries:
- matplotlib
- [geospacepy-lite](https://bitbucket.org/amienext/geospacepy-lite)

From the directory in which you have cloned the repo, issue:
```{sh}
python setup.py install
```

## Usage

```{python}
# Import
import apexpy
# Instantiate the class
# Epoch is the time frame for which you want to transform data
# altmin and altmax are the minimum and maximum valid altitudes for geographic locations
ac = apexpy.apex_converter(epoch=2010.5,altmin=300.,altmax=900.)
# Transform some data
lats = np.ones((10,))*80.
lons = np.linspace(-180.,180.,10)
alts = np.linspace(400.,600.,10)
alat,alon,qdlat = ac.geo2apex(lats,lons,alts)
```

See apex_converter_demo Jupyter notebook for more usage. Also see documentation.

Some later example in the notebook require [geospacepy-lite](https://github.com/lkilcommons/geospacepy-lite) for special plot types

## References
- Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
- VanZandt, T.E., W.L.Clark, and J.M.Warnock, Magnetic Apex Coordinates: A magnetic coordinate system for the ionospheric F2 layer, J Geophys. Res., 77, 2406-2411, 1972.

## Acknowledgements
- Liam Kilcommons (University of Colorado, Aerospace Engineering Sciences), took over this project and extended Peter Schmitt's port with the apex-converter class (apexpy.apex_converter)
       
The following people have helped review and update this project:
- Roy Barnes
- Astrid Maute
- Peter Schmitt
    
## Copyright and License
This software is part of the NCAR APEX coordinate system
transformation software.  Use is governed by the Open Source
Academic Research License Agreement contained in the file LICENSE.
