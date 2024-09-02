import subprocess
import sys

try:
    import numpy
except ImportError:
    print("Numpy is not installed. Installing now...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "numpy==1.26.3"])
    import numpy

try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    print("Setuptools is not installed. Installing now...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "--upgrade", "pip", "setuptools"])
    from setuptools import setup, find_packages, Extension

try:
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext
except ImportError:
    print("Cython is not installed. Installing now...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "Cython==3.0.8"])
    from Cython.Build import cythonize
    from Cython.Distutils import build_ext

# Include directories for the Cython extensions
extra_include_dirs = [numpy.get_include()]

# Extra arguments for the extension modules
ext_args = dict(
    include_dirs=extra_include_dirs,
    extra_compile_args=['-O3'],
    extra_link_args=['-O3'],
    language='c++',
)

# Reading long descriptions from files
with open('README.md', encoding="utf-8") as readme_file:
    readme = readme_file.read()

# Define the Cython extensions
ext_modules = cythonize([
    Extension(
        "Pythonreflectivity",  # Name of the extension
        sources=["UTILS/Pythonreflectivity.pyx"],  # Path to your .pyx file
        **ext_args
    ),
])

# Define requirements
requirements = [
    'Cython==3.0.8',
    'numpy==1.26.3',
    'h5py==3.10.0',
    'numba==0.59.0',
    'pyqtgraph==0.13.3',
    'scipy==1.12.0',
    'matplotlib==3.8.1',
    'PyQt5==5.15.2'
]

# Define setup and test requirements
setup_requirements = ['Cython>=3.0.8', 'pytest-runner']
test_requirements = ['pytest>=3', ]

# Setup configuration
setup(
    author="Lucas Korol, Robert J. Green, Jesus P. Curbelo, Raymond J. Spiteri",
    author_email='lsk601@usask.ca',
    python_requires='>=3.10',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'Topic :: Scientific/Engineering :: Physics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: BSD 3-Clause License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.10',
    ],
    description="GO-RXR: Global Optimization of Resonant X-ray Reflectometry Tool for Analyzing Quantum Materials",
    long_description=readme,
    long_description_content_type='text/markdown',
    install_requires=requirements,
    ext_modules=ext_modules,  # Include the Cython extensions here
    license="BSD 3-Clause License",
    include_package_data=True,
    keywords=[
        'resonant x-ray reflectometry',
        'global optimization',
        'thin films',
        'materials science',
        'quantum materials',
        'python'
    ],
    name='GO-RXR',
    packages=find_packages(include=['Magnetic_Scattering_Factor', 'Scattering_Factor', 'Testing', 'Tutorial']),
    setup_requires=setup_requirements,
    test_suite='TESTS',
    tests_require=test_requirements,
    url='https://github.com/lucaskorol21/GO-RXR',
    version='1.0',
    zip_safe=False,
    cmdclass={'build_ext': build_ext},
)