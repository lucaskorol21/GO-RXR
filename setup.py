from setuptools import setup, find_packages

setup(
    name='GO-RXR',
    version='0.3',
    packages=find_packages(include=['Magnetic_Scattering_Factor', 'Scattering_Factor', 'Testing', 'Tutorial']),
    url='https://github.com/lucaskorol21/GO-RXR',
    author='Lucas Korol',
    author_email='lsk601@usask.ca',
    maintainer='Robert J. Green',
    description='Global Optimization of Resonant X-ray Reflectometry Tool',
    python_requires='>=3.7',
    install_requires=[
        'Cython==3.0.2',
        'numpy==1.21.4',
        'PyQt5==5.15.7',
        'h5py==2.9.0',
        'matplotlib==3.4.3',
        'numba==0.55.2',
        'pyqtgraph==0.12.4',
        'scipy==1.7.1'
    ]
)
