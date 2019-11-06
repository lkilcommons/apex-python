#!/usr/bin/python
# Install script for Python bindings to Apex coordinate system transformation software.
import glob, os
import sys

from setuptools import setup, Extension
from setuptools.command import install as _install
from numpy.distutils.core import setup, Extension # Use numpy.distutils for this, it's much more full featured

os.environ['DISTUTILS_DEBUG'] = "1"

#The PyF wrapper and unit tests were made by Peter Schmitt
#Liam Kilcommons made the apex-converter vectorization and optimization
sources_f = glob.glob('apexpython/*.f')

#For some strange reason the compilation automation
#includes multiple copies of apex-f2pywrappers.o in the linking
#step if you don't remove it from the sources list, it must
#get automatically pulled in by including apex.pyf or something
#if 'apexpy/apex-f2pywrappers.f' in sources_f:
#	sources_f.remove('apexpy/apex-f2pywrappers.f')

apex = Extension(name = 'apexpython.apex',
                 sources = ['apexpython/apex.pyf','apexpython/apxntrp.f90'] + sources_f,
                 extra_f77_compile_args=['-g'])

setup(name = 'apexpython',
      version = '1.1',
      description       = 'Python bindings for Apex',
      author            = 'Peter Schmitt / Liam Kilcommons',
      author_email      = 'schmitt@ucar.edu / liam.kilcommons@colorado.edu',
      install_requires	= ['numpy'],
      packages			= ['apexpython'],
      ext_modules = [apex]
      )

