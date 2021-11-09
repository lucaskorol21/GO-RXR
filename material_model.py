import cython
import Pythonreflectivity as pr
import numpy
import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
import os

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
    idx = 0
    factor = f[0,0]  # First energy
    while (factor < E and factor != E):
        idx = idx + 1
        factor = f[idx, 0]

    if factor == E:
        f_real = f[idx, 1]
        f_imag = f[idx, 1]
    else:
        f_real = interpolate(f[idx - 1, 1], f[idx, 1], f[idx - 1, 0], f[idx, 0], E)
        f_imag = interpolate(f[idx - 1, 2], f[idx, 2], f[idx - 1, 0], f[idx, 0], E)

    return complex(f_real, f_imag)

def find_form_factor(element, E):
    """
    Purpose: Finds the form factors at energy 'E' from the form factor database for the specified element
    :param element: The name of the element in abbreviated from (string)
    :param E: Energy of scan
    :return: The complex form factor
    """
    for ifile in os.listdir(os.getcwd()):
        if ifile.startswith(element) and (ifile.endswith('.ff') or ifile.endswith('.ffm')):
            F = form_factor(np.loadtxt(ifile),E)

    return F

def dielectric_tensor(L, E):
    # energy [eV]
    # Constants
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e10  # Speed of light in vacuum [cm/s]
    re = 2.817940322719e-13  # Classical electron radius (Thompson scattering length) [cm]
    avocado = 6.02214076e23  # avagoadro's number
    k0 = 2 * pi * E / (h * c)  # photon wavenumber in vacuum [1/cm]
    factor = [0]
    constant = 2 * pi * re * (avocado) / (k0 ** 2)  # constant for density sum

    value = 0
    elements = list(L.keys())

    if len(elements) == 1:
        rho = L[elements[0]][1]
        F = find_form_factor(elements[0],E)
        value = F*L[elements[0]][1]
    else:
        idx = 0
        for element in elements:
            print(element)
            F = find_form_factor(element,E)
            value = value + F*L[element][1]
            idx =idx + 1




    n = 1 - constant*value

    epsilon = n ** 2

    return epsilon




if __name__ == "__main__":


    E = 799.82 # Xray Energy
    h = 4.135667696e-15  # Plank's Constant [eV s]
    c = 2.99792450e18  # Speed of light in vacuum [A/s]

    Theta = np.linspace(0.1, 89.9, 899)  # Angles
    wavelength = (h*c)/E  # Wavelength (same unit as roughness) (Angstroms or nm)
    test = np.loadtxt('test_example.txt')
    L1 = {'Fe': [0,0.028,0]}
    L2 = {'Fe':[38,0.028,0], 'La':[38,0.028,0]}
    epsilon_2 = dielectric_tensor(L2,E)
    epsilon_1 = dielectric_tensor(L1,E)
    A = pr.Generate_structure(2)
    A[0].seteps(epsilon_1)
    A[1].seteps(epsilon_2)
    A[1].setd(38)


    R1 = pr.Reflectivity(A, Theta, wavelength, MultipleScattering=True)

    plt.figure()
    qz = (0.001013546247)*E*sin(Theta*pi/180)
    Sigma, = plt.plot(qz, R1[0], 'k-',label='Python')
    Data, = plt.plot(test[:,0], test[:,1],'y--', label='ReMagX',)
    plt.legend(handles=[Sigma, Data])
    plt.yscale("log")
    plt.xlabel('qz')
    plt.ylabel('Reflectivity')
    plt.title('ReMagX vs. Python Script (800 eV)')
    plt.show()




