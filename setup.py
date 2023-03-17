import numpy

from distutils.core import setup, Extension
from Cython.Build import cythonize

setup(ext_modules = cythonize(Extension(
           "Pythonreflectivity",                                # the extesion name
           sources=["Pythonreflectivity.pyx"],
           language="c++",
           include_dirs=[numpy.get_include()]
      )))


