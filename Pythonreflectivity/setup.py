#<Pythonreflectivityl: A Python Package for simulation of x-ray reflectivities of Heterostructures>
#    Copyright (C) <2017>  <Martin Zwiebler>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU Lesser General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU Lesser General Public License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

#import numpy
#
#
#from distutils.core import setup
#from distutils.extension import Extension
#
#from Cython.Build import cythonize
#setup(ext_modules = cythonize(Extension(
#           "*",                                # the extension name
#           sources=["*.pyx"],
#           language="c++",
#           include_dirs=[numpy.get_include()]
#      )))
#

import numpy

##try:
##    from setuptools import setup
##    from setuptools import Extension
#except ImportError:
from distutils.core import setup
from distutils.extension import Extension

from Cython.Build import cythonize
setup(ext_modules = cythonize(Extension(
           "*",                                # the extension name
           sources=["*.pyx"],
           language="c++",
           extra_compile_args = ["-O3"],
           include_dirs=[numpy.get_include()]
      )))


