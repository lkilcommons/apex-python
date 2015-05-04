#!/usr/bin/python
# Install script for Python bindings to Apex coordinate system transformation software.
import glob
import sys

from numpy.distutils.core import setup, Extension

#The PyF wrapper and unit tests were made by Peter Schmitt
#Liam Kilcommons made the apex-converter vectorization and optimization
apex = Extension(name = 'apex-python.apex',
                 sources = ['apex-python/apex.pyf','apex-python/apxntrp.f90'] + glob.glob('apex-python/*.f'))

if __name__ == "__main__":
    setup(name = 'apex-python',
          version = '1.0',
          description       = 'Python bindings for Apex',
          author            = 'Peter Schmitt / Liam Kilcommons',
          author_email      = 'schmitt@ucar.edu / liam.kilcommons@colorado.edu',
          install_requires	= ['numpy','matplotlib'],
          packages			= ['apex-python'],
          ext_modules = [apex]
          )
