import cython
import Pythonreflectivity as pr
import numpy
import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
import os
from material_structure import *


def interpolate(fi,ff, Ei, Ef, E):
    """
    Purpose: Finds the form factor for a specific energy using linear interpolation
    :param fi: Initial form factor
    :param ff: Final form factor
    :param Ei: Initial energy
    :param Ef: Final energy
    :param E: Energy of desired form factor
    :return: Form factor at desired energy E
    """
    return ff - (ff-fi)*(Ef-E)/(Ef-Ei)

def form_factor(f,E):
    """
    Purpose: Determines if form factor can be directly taken from database or if linear interpolation needs to be done
    :param f: 2D list that contains the form factors at their respective energies
    :param E: Desired energy
    :return: From factor [real, imaginary]
    """
    from scipy import interpolate

    E_axis = f[:,0]
    tck_real = interpolate.splrep(E_axis,f[:,1])
    f_real = interpolate.splev(E, tck_real)

    tck_imag = interpolate.splrep(E_axis, f[:,2])
    f_imag = interpolate.splev(E, tck_imag)
    """
    idx = 0
    factor = f[0,0]  # First energy
    while (factor < E and factor != E):
        idx = idx + 1
        factor = f[idx, 0]

    # Determines if we need to interpolate between values
    if factor == E:
        f_real = f[idx, 1]  # real component
        f_imag = f[idx, 1]  # imaginary componentd
    else:
        f_real = interpolate(f[idx - 1, 1], f[idx, 1], f[idx - 1, 0], f[idx, 0], E)  # real component
        f_imag = interpolate(f[idx - 1, 2], f[idx, 2], f[idx - 1, 0], f[idx, 0], E)  # imaginary component
    """
    return complex(f_real, -f_imag)

def find_form_factor(element, E):
    """
    Purpose: Finds the form factors at energy 'E' from the form factor database for the specified element
    :param element: The name of the element in abbreviated from (string)
    :param E: Energy of scan
    :return: The complex form factor
    """
    F = 0
    my_dir = os.getcwd() + r'\Scattering_Factor'
    for ifile in os.listdir(my_dir):
        if ifile.endswith(element + '.txt'):
            F = form_factor(np.loadtxt(my_dir +  "\\" + ifile),E)
    #print(element, F)
    return F

def dielectric_constant(rho, sf, E):
    """
    Purpose: Compute the dielectric tensor constant
    :param L: Dictionary {Element 1: [density array], Element 2: [density array],..., Element N: {density array}}
    :param E: Energy of incoming photon (eV)
    :return: Dielectric constant (f_real + i*f_imag)
    """
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]

    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    value = 0
    elements = list(rho.keys())  # get's elements in layer
    F = dict()
    for element in elements:
        F[element] = find_form_factor(sf[element],E)


    if len(elements) == 1:
        value = F[elements[0]]*rho[elements[0]]  # computes alpha and beta values for form factor
    else:
        for element in elements:
            value = value + F[element]*rho[element]  # computes alpha and beta values

    n = 1 - constant*value  # computes the index of refraction
    epsilon = n ** 2  # computes dielectric constant from index of refraction
    #epsilon = n*np.conj(n)

    return epsilon




if __name__ == "__main__":


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
    eps = dielectric_constant(density, E)

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


