import cython
import Pythonreflectivity as pr
import numpy
import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
import os
from material_structure import *


def form_factor(f,E):
    """
    Purpose: Determines form factors with energy E using linear interpolation
    :param f: List of form factors
    :param E: Desired energy
    :return: Array contains real and imaginary component of form factor at energy E: f=[real, imaginary]
    """

    from scipy import interpolate

    """
    #E_axis = f[:,0]
    # Splice interpolation
    tck_real = interpolate.splrep(E_axis,f[:,1])
    f_real = interpolate.splev(E, tck_real)

    tck_imag = interpolate.splrep(E_axis, f[:,2])
    f_imag = interpolate.splev(E, tck_imag)

    """
    # Linear interpolation
    fr = interpolate.interp1d(f[:,0],f[:,1])
    fi = interpolate.interp1d(f[:,0],f[:,2])

    # Check if energy range is within the range specified by
    if f[0,1] > E or f[-1,1] < E:
        f_real = 0
        f_imag = 0
    else:
        f_real = fr(E)
        f_imag = fi(E)


    return [f_real, f_imag]

def find_form_factor(element, E, mag):
    """
    Purpose: Finds the form factors at energy 'E' from the form factor database for the specified element
    :param element: The name of the element in abbreviated from (string)
    :param E: Energy of scan
    :return: The complex form factor
    """
    F = 0

    if mag:
        my_dir = os.getcwd() + r'\Magnetic_Scattering_Factor'
    else:
        my_dir = os.getcwd() + r'\Scattering_Factor'

    for ifile in os.listdir(my_dir):

        if ifile.endswith(element + '.txt'):
            F = form_factor(np.loadtxt(my_dir +  "\\" + ifile),E)

    return F

def magnetic_optical_constant(rho, sf, E):
    mag = True
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    #value = 0
    f1 = 0
    f2 = 0
    elements = list(rho.keys())  # get's elements in layer
    F = dict()
    for element in elements:
        F[element] = find_form_factor(sf[element], E, mag)
    if len(elements) == 1:
        f1 = F[element[0]][0]*rho[element[0]]
        f2 = F[element[0]][1] * rho[element[0]]
        #value = F[element[0]] * rho[element[0]]
    else:
        for element in elements:
            f1 = f1 + F[element][0]*rho[element]
            f2 = f2 + F[element][1] * rho[element]
            #value = value + F[element] * rho[element]


    delta_m = constant * f1
    beta_m = constant * f2


    return delta_m, beta_m

def index_of_refraction(rho, sf, E):
    """
    Purpose: Compute the dielectric tensor constant
    :param L: Dictionary {Element 1: [density array], Element 2: [density array],..., Element N: {density array}}
    :param E: Energy of incoming photon (eV)
    :return: Dielectric constant (f_real + i*f_imag)
    """
    mag = False
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

    f1 = 0
    f2 = 0
    #value = 0
    elements = list(rho.keys())  # get's elements in layer
    F = dict()
    for element in elements:
        F[element] = find_form_factor(sf[element],E, mag)

    if len(elements) == 1:
        f1 = F[element[0]][0]*rho[element]
        f2 = F[element[0]][1] * rho[element]
        #value = F[element[0]]*rho[element[0]]
    else:
        for element in elements:
            #value = value + F[element]*rho[element]
            f1 = f1 + F[element][0]*rho[element]
            f2 = f2 + F[element][1]*rho[element]

    delta = f1*constant
    beta = f2*constant

    return delta, beta




if __name__ == "__main__":

    """
    E = 800 # Xray Energy
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e18  # Speed of light in vacuum [A/s]

    Theta = np.linspace(0.1, 89.9, 899)  # Angles
    wavelength = (h*c)/E  # Wavelength (same unit as roughness) (Angstroms or nm)
    test = np.loadtxt('test_example.txt')  # Data from ReMagX

    sample = slab(2)

    sample.addlayer(0, 'Fe', 50, density=1.56366)
    sample.addlayer(1, 'Fe', 38, density=1.56366)

    thickness, density, mag_density = sample.density_profile()
    eps = dielectric_constant(density, E, mag)

    A = pr.Generate_structure(2)  # initializes slab structure
    A[0].seteps(eps[0])  # creates the substrate layer
    A[1].seteps(eps[0])  # creates film layer
    A[1].setd(38)  # sets thickness


    R1 = pr.Reflectivity(A, Theta, wavelength, MultipleScattering=True)  # Computes the reflectivity

    plt.figure()
    qz = (0.001013546247)*E*sin(Theta*pi/180)
    Sigma, = plt.plot(qz, R1[0], 'k-',label='Python')
    plt.yscale("log")
    plt.xlabel('qz')
    plt.ylabel('Reflectivity')
    plt.title('ReMagX vs. Python Script (800 eV)')
    plt.show()
    """
    my_dir = os.getcwd() + r'\Magnetic_Scattering_Factor'
    file = my_dir + "\\" + "Ni.txt"
    print(file)
    np.loadtxt(file)
