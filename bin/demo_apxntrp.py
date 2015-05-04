#!/usr/bin/env python

import numpy
import os
import sys
try:
    # Find & import the Apex package:
    dirname = os.path.dirname( os.path.abspath(__file__) )
    prefix = os.path.split(dirname)[0]
    sys.path.insert(0, os.path.join(prefix, 'lib'))
    import apex
except ImportError:
    print "Error importing Apex module.  Did you build the Python bindings to Apex?"
    sys.exit(-1)



# due to a shortcomming in f2py (documented here http://cens.ioc.ee/pipermail/f2py-users/2008-December/001764.html), 
# we must know dimensions apriori.   Hard-coding them here:
nlat=13
nlon=21
nalt=2

workArray = numpy.empty(nlat*nlon*nalt*6, numpy.float32)
# Define grid points suitable for preparing interpolation
# tables for Apex coordinates, Modified Apex Coordinates and
# Quasi-Dipole coordinates
(gplat,gplon,gpalt) = apex.ggrid(nvert=4,
                                 glamn=-90.0, glamx=90.0,
                                 glomn=-180.0,   glomx=180.0,
                                 altmn=300.0,   altmx=300.0,
                                 nlat=nlat, nlon=nlon,  nalt=nalt)

# Create interpolation tables for a single time
ist = apex.apxmka(msgun=6,
                  epoch=1994.5,
                  gplat=gplat, gplon=gplon, gpalt=gpalt,
                  wk=workArray,
                  nlat=nlat,   nlon=nlon,   nalt=nalt)


# Apex to Geographic
(gdlat, gdlon, ist) = apex.apxa2g(alat=-13.840052604675293,
                                  alon=3.7327430248260498,
                                  alt=300.0,
                                  wk=workArray,
                                  lwk=len(workArray))

print 'After a2g, geographic coordinates are', gdlat, gdlon


# Modified Apex to Geographic
(gdlat, gdlon, ist) = apex.apxm2g(xlatm=-13.840052604675293,
                                  alon=3.7327430248260498,
                                  alt=300,
                                  hr=8.5,
                                  wk=workArray,
                                  lwk=len(workArray))

print 'After m2g, geographic coordinates are', gdlat, gdlon

# Quasi-Dipole to Geographic
(gdlat,gdlon, ist) = apex.apxq2g(qdlat=-6.8306035995483398,
                                 qdlon=3.7327430248260498,
                                 alt=300.0,
                                 wk=workArray,
                                 lwk=len(workArray))

print 'After q2g, geographic coordinates are', gdlat, gdlon
