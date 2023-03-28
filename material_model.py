"""
Library: material_model
Version: 0.1
Author: Lucas Korol
Institution: University of Saskatchewan
Last Updated: March 22nd, 2023
Python: version 3.7

Purpose: This script contains all functions related to retrieving form factors and calculating the structural and
        magnetic optical constants.

Imported Libraries

material_structure (version 0.1) - Part of the 'name' software package and is used to construct the sample model and calculate the
                     reflectivity spectra. This is only used for testing of the library functions

pickle - Pickle is used to load the data base of form factors, which are of a dictionary format.

numba (version 0.55.2) - Not currently used in this script.

scipy (version 1.7.1) - This library is used to interpolate the form factors for an arbitrary energy

os - used to access database of atomic masses
"""

from material_structure import *  # version
import pickle
from numba import *
from scipy import interpolate
import os


# Loads the non-magnetic form factors stored in our database
with open('form_factor.pkl', 'rb') as f:
    global ff
    ff = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
f.close()

# Loads the magnetic form factors stored in our database
with open('form_factor_magnetic.pkl','rb') as f:
    global ffm
    ffm = pickle.load(f)
f.close()

def _use_given_ff(directory):
    """
    Purpose: Scan cwd for form factors files and return their names (with .ff and .ffm stripped)
    :param directory: Current working directory
    :return: names of form factors (magnetic and structural) found in directory
    """
    global ff
    global ffm

    struct_names = []  # stores list of structural form factors
    mag_names = []  # stores list of magnetic form factors

    # loops through project directory for form factor files
    for file in os.listdir(directory):
        if file.endswith(".ff"):
            name = directory +'/' + file
            element = file.strip(".ff")

            struct_names.append(element)
            with open(name,'rb') as f:
                ff[element] = np.loadtxt(name)

            f.close()

        elif file.endswith(".ffm"):
            name = directory + '/' +file
            element = file.strip(".ffm")
            mag_names.append(element)

            with open(name, 'rb') as f:
                ffm[element] = np.loadtxt(name)
            f.close()

    return struct_names, mag_names

def form_factor(f,E):

    """
    Purpose: Determines form factors with energy E using linear interpolation
    :param f: List of form factors of form np.array([E, f_real, f_imag]) where E, f_real, and f_imag are arrays
    :param E: Float or list of desired energy
    :return: Array that contains the real and imaginary values of the form factor at energy E: f=[real, imaginary].
             If user inputs an array of energies [E_1,E_2,..,E_n] the output will be [[fr_1,fi_1],[fr_2,fi_2],...,[fr_n,fi_n]].
             Note - All values are real and are converted to imaginary numbers elsewhere.
    """

    # Linear interpolation
    fr = interpolate.interp1d(f[:,0],f[:,1])  # real component
    fi = interpolate.interp1d(f[:,0],f[:,2])  # imaginary component

    if isinstance(E, list) or isinstance(E, np.ndarray):  # handle multiple energy case (energy scan)
        F = np.array([np.array([fr(x), fi(x)]) if x > f[0, 0] and x < f[-1, 0] else np.array([0, 0]) for x in E])
    else:  # handle single energy case (reflectivity scan)
        F = np.array([fr(E), fi(E)]) if E>f[0,0] and E<f[-1,0] else np.array([0,0])

    return F

def find_form_factor(element, E, mag):
    """
    Purpose: Return the magnetic or non-magnetic form factor of a selected element and energy
    :param element: String containing element symbol
    :param E: Energy in electron volts
    :param mag: Boolean used to identify if non-magnetic or magnetic form factor is requested
    :return:
    """
    #global ffm
    #global ff

    if mag:  # magnetic form factor
        mag_keys = list(ffm.keys())
        if element not in mag_keys:
            raise NameError(element + " not found in magnetic form factors")
        F = form_factor(ffm[element],E)
    else:  # non-magnetic form factor
        struc_keys = list(ff.keys())
        if element not in struc_keys:
            raise NameError(element + " not found in structural form factors")
        F = form_factor(ff[element], E)

    return F

def MOC(rho, sfm, E, n):
    """
    Purpose: Computes the magneto-optical constant for the energy scan
    :param rho: dictionary containing the element symbol as the key and a numpy array as the value
    :param sfm: dictionary that contains the element symbol as the key and the absorptive and dispersive form factor components
    :param E: a numpy array containing energy values in eV
    :return: The absorptive and dispersive magnetic-optical constants
    """
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum



    elements = list(rho.keys())  # retrieves all the magnetic elements in the layer

    delta_m = np.array([np.zeros(n) for x in range(len(E))])  # pre-initialization
    beta_m = np.array([np.zeros(n) for x in range(len(E))])  # pre-initialization
    # Computes the dispersive and absorptive components of the magnetic-optical constant using list comprehensions
    for element in elements:
        delta_m = delta_m + np.array(
            [constant[x] * sfm[element][x, 0] * rho[element] for x in range(len(sfm[element][:, 0]))])
        beta_m = beta_m + np.array(
            [constant[x] * sfm[element][x, 1] * rho[element] for x in range(len(sfm[element][:, 1]))])


    return delta_m, beta_m

def magnetic_optical_constant(rho, sfm, E):
    """
    Purpose: Calculate the magnetic optical constants
    :param rho: Magnetic density in mol/cm^3
    :param sf: dictionary relates elements to scattering factor sf = {'ele1':'ffm1',...,'eleN':'ffmN'}
    :param E: Desired energy in units of eV
    :return: delta_m - magneto-optic dispersive component
             beta_m  - magneto-optic absorptive component
    """

    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * np.pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * np.pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    f1 = 0  # for dispersive component computation
    f2 = 0  # for absorptive component computation

    elements = list(rho.keys())  # retrieves all the magnetic elements in the layer

    # Computes the dispersive and absorptive components of the index of refraction
    for element in elements:
        f1 = f1 + sfm[element][0] * rho[element]
        f2 = f2 + sfm[element][1] * rho[element]



    delta_m = constant * f1  # dispersive component
    beta_m = constant * f2  # absorptive component


    return delta_m, beta_m

def IoR(rho,sf,E):
    """
    Purpose: compute the refractive index for multiple energies
    :param rho: dictionary containing element symbol as key and numpy array as value
    :param sf: dictionary containing element symbol as key and numpy array of dispersive and absorptive form factors
    :param E: numpy array of energies (eV)
    :return: The absorptive and dispersive components of the refractive index
    """
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * np.pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * np.pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    elements = list(rho.keys())  # retrieves all the magnetic elements in the layer
    delta = np.array([np.zeros(len(rho[elements[0]])) for x in range(len(E))])  # initialization
    beta = np.array([np.zeros(len(rho[elements[0]])) for x in range(len(E))])  # initialization

    # Computes the dispersive and absorptive components of the index of refraction using list comprehensions
    for element in elements:
        delta = delta + np.array([constant[x]*sf[element][x, 0] * rho[element] for x in range(len(sf[element][:, 0]))])
        beta = beta + np.array([constant[x] * sf[element][x, 1] * rho[element] for x in range(len(sf[element][:, 1]))])
    return delta, beta

def index_of_refraction(rho, sf, E):
    """
    Purpose: Calculates the dispersive and absorptive components of the index of refraction
    :param rho: Dictionary containing density profile of elements {'ele1':rho1,...,'eleN':rhoN}
    :param sf: Dictionary of scattering factors {'ele1':ff1, ... , 'eleN':ffN}
    :param E: Desired energy in units of eV
    :return: delta - dispersive component of the refractive index
             beta - absorptive component of the refractive index
    """

    mag = False  # statement for retrieval of non=magnetic form factors
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    #h = 4.1357e-15
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    #c = 2.9979e10
    #re = 2.8179e-13
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    #avocado = 6.0221e23
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]
    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    f1 = 0  # dispersive form factor
    f2 = 0  # absorptive form factor

    elements = list(rho.keys())  # retrieves element symbols within layer
    """
    F = dict()  # dictionary used to contain the form factors
    for element in elements:
        F[ element] = find_form_factor(sf[element],E, mag)
    """
    #  Computes the dispersive and absorptive components
    for element in elements:
        f1 = f1 + sf[element][0] * rho[element]
        f2 = f2 + sf[element][1] * rho[element]

    delta = f1*constant  # dispersive component
    beta = f2*constant  # absorptive component

    return delta, beta




if __name__ == "__main__":


    E = 800 # Xray Energy
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e18  # Speed of light in vacuum [A/s]

    Theta = np.linspace(0.1, 89.9, 899)  # Angles
    wavelength = (h*c)/E  # Wavelength (same unit as roughness) (Angstroms or nm)
    test = np.loadtxt('test_example.txt')  # Data from ReMagX

    sample = slab(2)


    sample.addlayer(0, 'SrTiO3', 50, density=1.56366)
    sample.addlayer(1, 'LaMnO3', 38, density=1.56366)
    sample.energy_shift()


    qz = (0.001013546247)*E*np.sin(Theta*np.pi/180)

    sample.reflectivity(E,qz)


