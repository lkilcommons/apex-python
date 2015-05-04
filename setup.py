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
apex = Extension(name = 'apexpy.apex',
                 sources = ['apexpy/apex.pyf','apexpy/apxntrp.f90'] + glob.glob('apexpy/*.f'))

setup(name = 'apexpy',
      version = '1.1',
      description       = 'Python bindings for Apex',
      author            = 'Peter Schmitt / Liam Kilcommons',
      author_email      = 'schmitt@ucar.edu / liam.kilcommons@colorado.edu',
      install_requires	= ['numpy'],
      packages			= ['apexpy'],
      ext_modules = [apex]
      )
