from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

setup(
    ext_modules=cythonize(Extension(
           "Pythonreflectivity",                                # the extension name
           sources=["Pythonreflectivity.pyx"],
           language="c++",
           include_dirs=[numpy.get_include()]
      ))

)