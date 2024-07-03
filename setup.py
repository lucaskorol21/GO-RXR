from setuptools import setup, find_packages

setup(
    name='GO-RXR',
    version='1.0',
    packages=find_packages(include=['Magnetic_Scattering_Factor', 'Scattering_Factor', 'Testing', 'Tutorial']),
    # packages=['Magnetic_Scattering_Factor', 'Scattering_Factor', 'Testing', 'Tutorial'],
    url='https://github.com/lucaskorol21/GO-RXR',
    author='Lucas Korol',
    author_email='lsk601@usask.ca',
    maintainer='Robert J. Green',
    description='Global Optimization of Resonant X-ray Reflectometry Tool',
    python_requires='>=3.7',
    install_requires=[
        'Cython==3.0.8',
        'numpy==1.26.3',
        'h5py==3.10.0',
        'numba==0.59.0',
        'pyqtgraph==0.13.3',
        'scipy==1.12.0',
        'matplotlib==3.8.1',
        'PyQt5==5.15.2'
    ]
)
