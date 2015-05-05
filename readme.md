# Apexpy
This is a set of python convenience tools that wrap Art Richmonds original Apex Fortran codes.

The F2Py extension bindings of the original tools is a apexpy.apex. This portion was
written by Peter Schmitt, formerly of NCAR High Altitude Observatory.

Liam Kilcommons (University of Colorado, Aerospace Engineering Sciences), took over this
project and extended Peter's port with the apex-converter class (apexpy.apex_converter).
This class simplify processing large amounts of spacecraft data by providing simplified interface
to the apex routines, along with convenience routines for converting from Apex Longitude to Magnetic Local Time,
finding the subsolar point, etc.

No modifications other than those in the following original README
have been made.

# Installation
You must have gfortran installed. Apex_converter requires numpy.
From the directory in which you have cloned the repo, issue:
`python setup.py install`


# Original README 
```
NAME
    Apex Coordinate System Transformations

VERSION
    version 0.0.0.  Revision 954 from TGCM repository at HAO. July 3, 2013.

AUTHOR
    Art Richmond

WEBSITE
    https://bugs.hao.ucar.edu/projects/apex

DESCRIPTION
    The Apex package implements transformations between the following
    coordinate systems:

      * Apex
      * Modified Apex
      * Quasi-Dipole Apex
      * Geographic

    The Apex package is written in Fortran 77 and ships with Python
    bindings.

IMPORTANT NOTE
    The Apex code has been under development for more than 20 years.
    The codebase was cleaned up in 2012 to prepare for an official 1.0.0
    release.  Part of this cleanup may have changed the interface
    slightly: The original unmodified implementation is at
    src/apxntrp_legacy.f.  The published API is available in
    src/apxntrp.f90.  This API is implemented in Fortran 90 and
    subroutines may have slightly different parameters from the legacy
    code.

REQUIREMENTS
    The Apex code has been tested on Linux (CentOS-5) and Mac (10.7.x)
    workstations and requires the following software:

    * GNU Make
    * A valid Fortran compiler (IBM XLF, gfortran, Portland Group or Intel) 
    * Optional Python bindings require:
       - Python >= 2.4
       - Numpy >= 1.4
       - Nose >= 0.11.3

INSTALLATION
    Edit Makefile to set the following variables: FCOMPILER, FC, OPT.
    Type `gmake` to get a list of build targets related to the Apex
    package.

    UNSUPPORTED: Edit src/apxntrp_legacy.f to set the Inter Record Length Factor (IRLF).
    The IRLF is a system-dependent value used for unformatted read &
    writes used by the APEX code. IRLF is 1 when RECL units are words,
    otherwise it is a function of the computer word length; e.g.,

          IRLF = 1 on DEC           (32-bit words, RECL units are words)
          IRLF = 1 on SGI using f77 (32-bit words, RECL units are words)
          IRLF = 4 on SGI using f90 (32-bit words, RECL units are bytes)
          IRLF = 4 on PC            (32-bit words, RECL units are bytes) **default**
          IRLF = 4 on Sun           (32-bit words, RECL units are bytes)
          IRLF = 8 on Cray          (64-bit words, RECL units are bytes)

VERIFICATION
    `gmake test` will build & execute example & test programs.  Python
    bindings are required as we use the Python unittest framework to

    manage the tests.  You can tell if the tests ran succesfully if
    you see output like:

       ......
       ----------------------------------------------------------------------
       Ran 6 tests in 0.442s
       
       OK

   The tests had a problem if you see either "E" or "F" and a stack
   trace followed by the text "FAILED".  Here's an example of a failed
   tests:

       ...F.E
       ======================================================================
       ERROR: Quasi-Dipole Apex to Geographic (works with global grid only)
       ----------------------------------------------------------------------
       Traceback (most recent call last):
         File "/Users/schmitt/src/apex/trunk/src/test_apxntrp.py", line 104, in test_APXQ2G
           (gdlat,gdlon, ist) = apex.apxq2gzzz(qdlat=-6.8306035995483398,
       AttributeError: 'module' object has no attribute 'apxq2gzzz'
       
       ======================================================================
       FAIL: Modified Apex to Geographic (works with global grid only)
       ----------------------------------------------------------------------
       Traceback (most recent call last):
         File "/Users/schmitt/src/apex/trunk/src/test_apxntrp.py", line 96, in test_APXM2G
           self.assertAlmostEqual(-4.45752, gdlat, places=4)
       AssertionError: -4.45752 != -4.825495719909668 within 4 places
       
       ----------------------------------------------------------------------
       Ran 6 tests in 0.436s
       
       FAILED (errors=1, failures=1)

USAGE
    Supported routines are:
      * ggrid: generate grids for transformations
      * mka: create interpolation tables
      * ggc: get grid coordinates
      * apxall: Transform geographic to Apex
      * apxa2g: Transform Apex to geographic
      * apxm2g: Transform Modified Apex to geographic
      * apxq2g: Transform Quasi-dipole Apex to geographic 
      * apxmall: Transform Geographic to Modified Apex & Quasi-dipole Apex

    See test programs at
       src/test_apxntrp.f
       src/test_apxntrp.py
       bin/demo_apxntrp.py
    for example usage, required variables & calling sequence.

    Additional documentation may be found as comments in the Fortran
    source code.  See src/apxntrp.f90.  You may find more verbose
    documentation at src/apxntrp_legacy.f.

    To use the Apex routines in your existing Fortran code, compile
    libapex.a (ie. enter `gmake libapex.a`) and link against the Apex
    library: Append to your link flags (ie. LDFLAGS):
        -L /path/to/apex/lib -lapex

    Please note as of version 1.0.0, the routines APXRDA and APXWRA
    are old implementations and have not been tested.  We do not
    recommend use with APXRDA/APXWRA.  A future release will reimplement
    these subroutines with NetCDF or HDF I/O.

REFERENCES
    Richmond, A. D., Ionospheric Electrodynamics Using Magnetic Apex
       Coordinates, J. Geomag. Geoelectr., 47, 191-212, 1995.
    VanZandt, T.E., W.L.Clark, and J.M.Warnock, Magnetic Apex
       Coordinates: A magnetic coordinate system for the ionospheric
       F2 layer, J Geophys. Res., 77, 2406-2411, 1972.

ACKNOWLEDGEMENTS
    The following people have helped review and update this project:
       Roy Barnes
       Astrid Maute
       Peter Schmitt

COPYRIGHT AND LICENSE
    This software is part of the NCAR APEX coordinate system
    transformation software.  Use is governed by the Open Source
    Academic Research License Agreement contained in the file LICENSE.

SEMANTIC VERSIONING
    This is versioned in accordance with semantic versioning  <http://semver.org/>.
```