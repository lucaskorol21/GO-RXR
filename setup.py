from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy


setup(
    name='GO-RXR',
    version='0.3',
    packages=[''],
    url='https://github.com/lucaskorol21/GO-RXR',
    license='',
    author='Lucas Korol',
    author_email='lsk601@usask.ca',
    maintainer='Robert Green',
    description='Global Optimization of Resonant X-ray Reflectometry Tool',
    python_requires='==3.7',
    install_requires=[
        'Cython==0.29.24',
        'PyQt5==5.15.7',
        'h5py==2.9.0',
        'matplotlib==3.4.3',
        'numba==0.55.2',
        'pyqtgraph==0.12.4',
        'scipy==1.7.1'
    ],
    ext_modules=cythonize(Extension(
           "Pythonreflectivity",                                # the extension name
           sources=["Pythonreflectivity.pyx"],
           language="c++",
           include_dirs=[numpy.get_include()]
      ))

)
