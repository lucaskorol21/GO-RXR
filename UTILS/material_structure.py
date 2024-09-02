"""
Library: material_structure
Version: 1.01
Author: Lucas Korol
Institution: University of Saskatchewan
Last Updated: October 11,2023
Python: version 3.7

Purpose: This script contains all functions related to retrieving form factors and calculating the structural and
        magnetic optical constants.

Imported Libraries

numpy (version 1.21.4) - used for array manipulation

material_model (version 0.1) - Part of the 'name' software package and is used to retrieve the form factors and
                               calculate the optical constants

Pythonreflectivity (version 1.0) - Software used to calculate the reflectivity spectra. Currently this software is only
                                   compatible with python 3.7. The function generate_structure, slab.reflectivity, and
                                   slab.energy_scan will need to be altered for newer versions.

numba (version 0.55.2) - Used to speed up certain functions with long executions times and frequent use.

scipy (version 1.7.1) - This python file uses the error function from the scipy library to model the roughness of a
                        material

"""

import matplotlib.pyplot as plt
import numpy as np
from UTILS.material_model import *
import Pythonreflectivity as pr
from numba import *
from collections import OrderedDict
from scipy.special import erf
import warnings
import copy
import os

# Import ROOT_PATH from the __init__.py file
from . import ROOT_PATH

def find_ff(element,E, ff_dict):
    """
    Purpose: Return the form factor from material_model.py
    :param element: The element symbol
    :param E: The energy in units of electronvolts (can be entered as a list of energies)
    :param ff_dict: form factor dictionary
    :return: The real and imaginary component of the form factor {can be a list of tuples (real, imaginary)}
    """
    F = form_factor(ff_dict[element], E)
    return F

@njit()
def ALS(beta, delta, beta_m, delta_m, precision=1e-6):
    """
    Purpose: Return a list of the indices for the adaptive layer segmentation
    :param alpha: numpy array of refractive component
    :param beta: numpy array of absorptive form factor values
    :param precision: precision value used in slicing
    :return: my_slabs - contains indices for slicing
    """

    beta = beta /np.linalg.norm(beta)  # normalize the refractive component
    delta = delta /np.linalg.norm(delta)  # normalize the absorptive component
    beta_m = beta_m/np.linalg.norm(beta_m)
    delta_m = delta_m / np.linalg.norm(delta_m)

    idx_a = 0  # keeps track of surface of previous slab

    n = beta.size
    my_slabs = np.zeros(n) # pre-initialize the slab array to the maximum number of slabs possible

    dsSlab = 1  # counts the number of slices
    for idx_b in range(1,n):

        # retrieves optical values
        f1 = beta[idx_a]
        f2 = beta[idx_b]
        f1d = delta[idx_a]
        f2d = delta[idx_b]

        # retrieved magneto-optical values
        f1m = beta_m[idx_a]
        f2m = beta_m[idx_b]
        f1dm = delta_m[idx_a]
        f2dm = delta_m[idx_b]

        var_b = np.absolute(f2-f1)  # varitation of refractive component
        var_d = np.absolute(f2d-f1d)  # variation of absorptive component
        var_bm = np.absolute(f2m - f1m)  # varitation of the magnetic refractive component
        var_dm = np.absolute(f2dm - f1dm)  # variation of magnetic absorptive component

        # checks if variation meets precision value if not check next value
        if var_b>precision or var_d>precision or var_bm>precision or var_dm>precision:
            my_slabs[dsSlab] = idx_b
            idx_a = idx_b  # change previous slice location
            dsSlab = dsSlab + 1
        elif idx_b == n-1:  # reached the end of the array
            my_slabs[dsSlab] = idx_b
            dsSlab = dsSlab + 1

    my_slabs = my_slabs[:dsSlab]
    return my_slabs

def generate_structure(thickness, structure, my_slabs, epsilon, epsilon_mag, layer_magnetized, transition):
    """
    Purpose: Generate the object structure as defined by pythonreflectivity
    :param thickness: thickness numpy array of length n
    :param structure: Array of length m containing information of each layer
    :param my_slabs: Slab indices for adaptive layer segmentation
    :param epsilon: permittivity numpy array of length n
    :param epsilon_mag: magnetic permittivity numpy array of length n
    :param layer_magnetized: numpy array of length m containing booleans that define if layer is magnetized
    :param transition: array that containes the thickness at which a layer transition occurs
    :return: The object structure as defined by pythonreflectivity
    """
    m = len(my_slabs)  # number of slabs
    # m = len(epsilon)
    A = pr.Generate_structure(m)  # creates object for reflectivity computation
    m_j = 0  # previous slab
    idx = 0  # keeps track of current layer
    layer = 0
    gamma = 90  # pre-initialize magnetization direction
    phi = 90

    # Recent versions of Pythonreflectivity use the susceptibility instead of the dielectric constant
    for m_i in my_slabs:
        d = thickness[m_i] - thickness[m_j]  # computes thickness of slab
        #eps = (epsilon[m_i] + epsilon[m_j]) / 2  # computes the dielectric constant value to use
        eps = epsilon[m_j]  # computes the dielectric constant value to use
        #eps_mag = (epsilon_mag[m_i] + epsilon_mag[m_j]) / 2
        eps_mag = epsilon_mag[m_j]  # computes the magnetic dielectric constant

        # Determines the magnetization direction of the first layer
        if layer == 0:
            if layer_magnetized[0]:
                for ele in structure[layer].keys():
                    if len(structure[0][ele].mag_scattering_factor) != 0:
                        gamma = structure[0][ele].gamma
                        phi = structure[0][ele].phi

        # Determines the magnetization direction of the other layers
        temp_list = [transition[i][layer] for i in range(len(transition))]
        trans = max(temp_list)
        if trans <= thickness[m_j] and layer < len(transition[0]) - 1:
            layer = layer + 1
            if layer_magnetized[layer]:
                for ele in structure[layer].keys():
                    if len(structure[layer][ele].mag_scattering_factor) != 0:
                        gamma = structure[layer][ele].gamma
                        phi = structure[layer][ele].phi

        # sets the magnetization direction based on the input angles
        #  - in future versions we will be able to set the magnetization direction based on phi and theta
        if layer_magnetized[layer]:
            if gamma == 90 and phi == 90:
                A[idx].setmag('y')
            elif gamma == 0 and phi == 90:
                A[idx].setmag('x')
            elif gamma == 0 and phi == 0:
                A[idx].setmag('z')
            else:
                raise ValueError('Values of Gamma and Phi can only be (90,90), (0,90), and (0,0)')

            A[idx].seteps([eps, eps, eps, eps_mag])  # sets dielectric tensor for magnetic layer
        else:
            A[idx].seteps(eps)  # sets dielectric tensor for non-magnetic layer

        if idx != 0:
            A[idx].setd(d)  # sets thickness of layer if and only if not substrate layer

        # move onto the next layer
        m_j = m_i
        idx = idx + 1

    return A

def energy_reflectivity(A, Theta, wavelength, R, E, backS=0, scaleF=1):
    """
    Purpose: Compute the reflectivity of a specific energy for the energy scan
    :param A: Structure object as defined by pythonreflectivity
    :param Theta: Angle of energy scan in degrees
    :param wavelength: Wavelength of photon energy
    :param R: Pre-initialized numpy array
    :param E: Integer that defined energy location for pre-initialized reflectivity dicitonary R
    :param backS: Float value that contains the background shift applied to the scan
    :param scaleF: Float value that contains the scaling factor applied to the scan
    :return: Reflectivity R
    """

    Rtemp = pr.Reflectivity(A, Theta, wavelength, MagneticCutoff=1e-10)  # Computes the reflectivity

    # sets the polarization reflectivity array
    if len(Rtemp) == 2:
        R['S'][E] = Rtemp[0][0]*scaleF + backS  # s-polarized light
        R['P'][E] = Rtemp[1][0]*scaleF + backS  # p-polarized light
        R['AL'][E] = scaleF*(Rtemp[0][0] - Rtemp[1][0]) / (scaleF*(Rtemp[0][0] + Rtemp[1][0])+backS*2)  # Asymmetry linear polarized
    elif len(Rtemp) == 4:
        delta_e = 1e-6
        R['S'][E] = Rtemp[0][0]*scaleF + backS  # s-polarized light
        R['P'][E] = Rtemp[1][0]*scaleF + backS  # p-polarized light
        R['AL'][E] = scaleF*(Rtemp[0][0] - Rtemp[1][0]) / (scaleF*(Rtemp[0][0] + Rtemp[1][0])+backS*2)  # linear asymmetry
        R['LC'][E] = Rtemp[2][0]*scaleF + backS  # left circular polarization
        R['RC'][E] = Rtemp[3][0]*scaleF + backS  # right circular polarization
        R['AC'][E] = scaleF*(Rtemp[2][0] - Rtemp[3][0]) / (scaleF*(Rtemp[2][0] + Rtemp[3][0]))+2*backS   # circular asymmetry
    else:
        raise TypeError('Error in reflectivity computation. Reflection array not expected size.')
    return R


def get_number(string):
    """
    Purpose: Strip successive digits from a string. This function is used to interpret chemical formulas
    Input:
     string - A string containing letters and digits
    Output:
     num - The stripped numbers in integer form
     string - The same string with the successive digits removed

    """

    finish = True  # Boolean used to determine if end of successive digits
    mynumbers = list()  # Stores the successive digits to be joined later
    while len(string)>0 and finish:
        # Determines if left most character is a digit
        if string[0].isdigit():
            mynumbers.append(string[0])  # appends digit to character list
            string = string[1:]  # removes digit from string
        else:
            finish = False
    num = int(''.join(mynumbers))  # joins digits into integer type

    return string, num




def checkstring(formula):
    """
    Purpose: This function identifies the elements and their stoichiometric relations
    :param formula: A string that represents the chemical formula. It is required that the chemical formula always
                    starts with a capital letter. Refer below for documentation description:
                        U - Uppercase letter
                        L - Lowercase letter
                        N - Digit

    :return: Returns a list that contains the chemical formula and it's stoichiometry [symbol, stoich]
    """
    # including dummy variables A and X
    # Used in error check
    symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
               'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se' , 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
               'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po' , 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
               'Rg', 'Cn', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
               'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'A', 'D','E','G','J','L','M','Q','R','T','X','Z']

    # The letters A, D, E, G, J, L, M, Q, R, T, X, and Z represent 'dummy variables'

    info = []  # list to store the element symbol and it's stoichiometry
    n = len(formula)  # Get's the number of characters in the list

    if n <= 2:  # Base Case
        if n == 2:  # Case for UL or UN
            if formula[1].isdigit():  # Case UN
                info.append(formula[0])
                info.append(int(formula[1]))
                formula = ''  # Sets formula for end case
            elif formula[1].islower():  # Case UL
                info.append(formula)
                info.append(1)
                formula = ''  # Sets formula for end case
            elif formula[1].isupper():
                info.append(formula[0])
                info.append(1)
                formula = formula[1]
        else:  # Case U
            info.append(formula)
            info.append(1)
            formula = ''

    elif formula[0].isupper():  # All other cases
        if formula[1].islower():
            if formula[2].isdigit():  # Case ULD+
                info.append(formula[0:2])
                formula, num = get_number(formula[2:])
                info.append(num)
            elif formula[2].islower():
                raise NameError('Unexpected lowercase character ' + formula[0:3])
            else:  # Case UL
                info.append(formula[0:2])
                info.append(1)
                formula = formula[2:]

        elif formula[1].isdigit():  # Case UD+
            info.append(formula[0])
            formula, num = get_number(formula[1:])
            info.append(num)

        else:  # Case U
            info.append(formula[0])
            info.append(1)
            formula = formula[1:]
    else:  # Raises error if unexpected symbol
        raise NameError('Unexpected symbol ' + formula[0])

    if info[0] not in symbols:
        raise NameError('Symbol ' + info[0] + ' does not correspond to an element found on the periodic table.')
    return formula, info


def find_stoichiometry(formula):
    """
    Purpose: Determine the stoichiometry of the formula inputted
    :param formula: String that represents a chemical formula
    :return: dictionary that contains element symbol as key and element class as value
    """
    num_elements = 0  # keeps track of the number of elements in the chemical formula
    ele = []  # keeps track of elements in the chemical formula
    data = []  # keeps track of the data information
    repeat_dict = dict()  # dictionary that keeps track of the repeated elements
    mydict = dict()  # data dictionary

    # Loops through the formula string
    while len(formula) > 0:
        formula, info = checkstring(formula)  # Get's current element
        if info[0] in ele:  # checks if element is repeated
            repeat_dict[info[0]] = 1  # sets the string counter to 1

        ele.append(info[0])  # adds new elements to element list
        data.append([info[0], element(info[0], info[1])])  # adds element information
        num_elements = num_elements + 1  # counts number of elements in chemical formula

    # Loops through all the components to change repeated element keys
    for component in data:
        if component[0] in list(repeat_dict.keys()):  # case of repeated element
            old_key = component[0]  # get's the element symbol
            component[0] = str(old_key) + str(repeat_dict[old_key])  # adds number to string to distinguish between 'different' elements
            repeat_dict[old_key] = repeat_dict[old_key] + 1  # adds number for new distinguishable element


    mydict = dict(OrderedDict(data))


    return mydict, num_elements

def perovskite_density(formula):
    """
    Purpose: Retrieves the density (g/cm^3) for common perovskite materials
    :param formula: Chemical formula of perovskite materials
    :return: Density of perovskite material (g/cm^3)
    """
    density = None
    # file = open("Perovskite_Density.txt", "r")
    file = open(os.path.join(ROOT_PATH, 'DATA', "Perovskite_Density.txt"), "r")
    lines = file.readlines()
    for line in lines:
        if line.split()[0] == formula:
            density = line.split()[1]
    file.close()

    if density == None:
        raise NameError("Inputted formula not found in perovskite density database")

    return float(density)

def atomic_mass(atom):
    """
    Purpose: Returns the molar mass of elements in the periodic table
    :param atom: Element symbol
    :return: Molar mass (g/mol)
    """
    mass = None
    
    # file = open("Atomic_Mass.txt", "r")
    file = open(os.path.join(ROOT_PATH, 'DATA', "Atomic_Mass.txt"), "r")
    lines = file.readlines()
    for line in lines:
        if line.split()[0] == atom:
            mass = line.split()[1]
    file.close()

    if mass == None:
        raise NameError("Inputted formula not found in perovskite density database")

    return float(mass)

def error_function(t, sigma1, sigma2, offset1, offset2):
    """
    Purpose: Calculates the error function
    :param t: thickness
    :param sigma1: interface 1 roughness value
    :param sigma2: interface 2 roughness value
    :param offset1: interface 1 offset thickness value
    :param offset2:  interface 2 offset thickness value
    :return:
    """
    result1 = (erf((t-offset1)/sigma1/np.sqrt(2))+1)/2
    result2 = (erf((offset2-t)/sigma2/np.sqrt(2))+1)/2
    result = result1*result2

    # This function is not used for the error function calculation
    return result

class element:
    def __init__(self, name, stoichiometry):
        """
        Purpose: Used in class slab. Used to keep track of element properties for each layer
        :param name: Element Symbol
        :param stoichiometry: The stoichiometric relation of the chemical formula as input into the class slab
        """
        self.name = name  # Element Symbol
        self.molar_mass = 0  # Total molar mass of all elements in slab layer
        self.density = 0  # Density of the molecule (g/cm^3)
        self.thickness = 0  # Thickness of the layer (Angstrom)
        self.roughness = 0  # Roughness of the surface (Angstrom)
        self.linked_roughness = None
        self.stoichiometry = stoichiometry  # Stoichiometry of the element
        self.poly_ratio = 1  # List that keeps track of polymorphous density ratio
        self.polymorph = []  # A list that contains the 'names' of various forms of the element (e.g. oxidation states)
        self.gamma = 0  #
        self.phi = 0  #
        self.mag_density = []  # The scalling factor we want to multiply our scattering factor by (density is not the correct description)
        self.scattering_factor = name  # Identifies which scattering factor to be used. This parameter will allow us to implement 'scattering functions'
        self.mag_scattering_factor = []
        self.position = None

class slab:
    def __init__(self, num_layers):
        """
        Purpose: Create a slab model of the sample layer. Keeps track of the elemental properties for each layer
        :param num_layers: Initializes the sample with desired number of slabs/layers. Must be inputted as an integer.
        """
        self.link = []  # Describes if elements in current layer have a linked roughness to the next layer. This is a list that contains Boolean type.
        self.myelements = []  # Keeps track of the elements in the material

        # Checks if num_layers is integer type
        if type(num_layers) != int:
            raise NameError('Number of slab layers must be entered as an integer type')
        if num_layers < 0:
            warnings.warn("Number of slab layers should be inputted as a positive integer. "
                          "Code assumed that negative value was intentional and took it's absolute value.")
            num_layers = abs(num_layers)

        self.structure = [dict() for i in range(num_layers)]  # keeps track of the structural properties of each layer
        self.poly_elements = dict()  # Keeps track of the  polymorphous elements
        self.mag_elements = dict()  # Keeps track of the magnetized elements in the material
        self.number_layers = num_layers
        self.find_sf = [dict(), dict()]  # [structural, magnetic]
        self.transition = None
        self.layer_magnetized = [False for i in range(num_layers)]  # keeps track of layers with magnetization
        self.eShift = dict()
        self.mag_eShift = dict()
        self.ff_scale = dict()  # keep track of the structural scaling factors
        self.ffm_scale = dict()  # keep track of the magnetic scaling factors



    def addlayer(self, num_layer, formula, thickness, density=None, roughness=0 , linked_roughness = None, sf = None):
        """
        Purpose: Add layer to sample material
        :param num_layer: Layer number (0 indicates substrate layer)
        :param formula: Chemical formula
        :param thickness: Thickness of the material in Angstrom 10^(-10)m
        :param density: Density of the perovskite material (g/cm^3). No input indicates user wants to use density found in database
        :param roughness: Rougness at the interface
        :param link: Determines if roughness between two of the same sites are linked
        :param gamma: Gamma value for magnetization of entire slab
        :param theta: Theta value for magnetization of entire slab
        """
        # Retrieve the element info
        elements, num_elements = find_stoichiometry(formula)
        molar_mass = 0
        # -------------------------------------- Error Checks ---------------------------------------------------------#

        # Let's user know that layer set is out of total number of layer range
        if num_layer > len(self.structure) - 1:
            raise RuntimeError( 'Layer ' + str(num_layer) + ' selected, but maximum layer is ' + str(len(self.structure)-1))

        # Let's user know layer already initialized
        if len(self.structure[num_layer]) != 0:
            warnings.warn('Layer '+str(num_layer)+' already initialized')

        # Checks Thickness
        if type(thickness) == int or type(thickness) == float:
            temp_thickness = [thickness for i in range(num_elements)]
            thickness = temp_thickness
        if type(thickness) != list and type(thickness) != np.ndarray:
            raise TypeError(
                'Layer ' + str(num_layer) + ': Thickness must be of float, integer, list, or numpy.ndarray type.')
        else:
            for idx in range(len(thickness)):
                if thickness[idx] < 0:
                    thickness[idx] = abs(thickness[idx])
                    warnings.warn('Layer ' + str(
                        num_layer) + ': Thickness should be entered as a positive value. The absolute value was taken as it was assumed the user meant to input a negative value.')
                if thickness[idx] > 100:
                    warnings.warn('Layer ' + str(num_layer) + ': Thickness is much larger than expected')

        # Checks Density
        if density == None:
            pass
        elif type(density) != int and type(density) != float and type(density) != list and type(density) == np.ndarray:
            raise TypeError('Layer ' +str(num_layer)+': Density must be entered in as a float, integer, list, or numpy.ndarray')
        elif type(density) == int or type(density) == float:
            if density < 0:
                warnings.warn('Layer ' + str(num_layer) + ': Density must be positive. The absolute value was taken as it was assumed the user meant to input a negative value.')
                density = abs(density)
            elif density >20:
                warnings.warn('Layer ' + str(num_layer) + ': The density of ' + str(density) + ' g/cm^3 might be too large of a value. Consider double checking your density. ')
        elif type(density) == list or type(density) == np.ndarray:
            for idx,dens in enumerate(density):
                if dens < 0:
                    warnings.warn('Layer ' + str(
                        num_layer) + ': Density must be positive. The absolute value was taken as it was assumed the user meant to input a negative value.')
                    density[idx] = abs(dens)
                elif dens > 20:
                    warnings.warn('Layer ' + str(num_layer) + ': The density of ' + str(
                        density) + ' g/cm^3 might be too large of a value. Consider double checking your density. ')
        elif type(roughness) == list or type(roughness) == np.ndarray:
            if len(density) != num_elements:
                raise ValueError('Density array must have same number of elements as does chemical formula')

            for idx in range(len(density)):
                if density[idx] < 0:
                    density[idx] = abs(roughness[idx])
                    warnings.warn('Layer ' + str(num_layer) + ': Roughness should be entered as a positive value. The absolute value was taken as it was assumed the user meant to input a negative value.')
                if density[idx] > 20:
                    warnings.warn('Layer ' + str(num_layer) + ': Roughness is much larger than expected')



        # Checks Roughness
        if type(roughness) == int or type(roughness) == float:
            temp_roughness = [roughness for i in range(num_elements)]
            roughness = temp_roughness
        if type(roughness) != list and type(roughness) != np.ndarray:
            raise TypeError('Layer ' +str(num_layer)+': Roughness must be of float, integer, list, or numpy.ndarray type.')
        else:
            for idx in range(len(roughness)):
                if roughness[idx] < 0:
                    roughness[idx] = abs(roughness[idx])
                    warnings.warn('Layer ' + str(num_layer) + ': Roughness should be entered as a positive value. The absolute value was taken as it was assumed the user meant to input a negative value.')
                if roughness[idx]>15:
                    warnings.warn('Layer ' + str(num_layer) + ': Roughness is much larger than expected')



        # Determine if current layer is linked to the previous layer
        if linked_roughness == None:  # if no inputs set an array with False booleans with length num_elements
            linked_roughness = [False for i in range(num_elements)]
        elif type(linked_roughness) != list and type(linked_roughness) != np.ndarray and type(linked_roughness) != int and type(linked_roughness) != float:
            raise TypeError('Variable linked_roughness must be of type int, float, list,  or numpy.ndarray.')
        elif type(linked_roughness) ==  float or type(linked_roughness) == int:
            temp_linked_roughness = [linked_roughness for i in range(num_elements)]
            linked_roughness = temp_linked_roughness
        elif len(linked_roughness) != num_elements:
            ValueError('linked_roughness must have the same number of elements as chemical formula.')
        if type(roughness) != list and type(roughness) != np.ndarray:
            raise TypeError('Layer ' + str(num_layer) + ': Roughness must be of float, integer, list, or numpy.ndarray type.')
        else:
            for idx in range(len(linked_roughness)):
                if linked_roughness[idx] < 0:
                    linked_roughness[idx] = abs(linked_roughness[idx])
                    warnings.warn('Layer ' + str(num_layer) + ': Roughness should be entered as a positive value. The absolute value was taken as it was assumed the user meant to input a negative value.')
                if linked_roughness[idx] > 15:
                    warnings.warn('Layer ' + str(num_layer) + ': Roughness is much larger than expected')

        # -------------------------------------------------------------------------------------------------------------#
        # Retrieve the element info
        elements, num_elements = find_stoichiometry(formula)
        molar_mass = 0
        # Computes the total molar mass of perovskite material
        for key in list(elements.keys()):
            molar_mass = molar_mass + atomic_mass(elements[key].name)*elements[key].stoichiometry


        # ------------------------------------------------------------------------------------------------------------ #
        # ------------------------------------------------------------------------------------------------------------ #
        # might also need to check if last layer
        # ------------------------------------------------------------------------------------------------------------ #
        # ------------------------------------------------------------------------------------------------------------ #

        position = 0
        for key in list(elements.keys()):

            if key not in self.myelements:  # Adds element to my element list if not already there
                self.myelements.append(key)

            if density == None:
                density = perovskite_density(formula)
                elements[key].density = density * elements[key].stoichiometry / molar_mass  # sets density  (g/cm^3)
            elif type(density) == float or type(density) == int:
                elements[key].density = density* elements[key].stoichiometry / molar_mass # sets density  (g/cm^3)
            elif type(density) == list or type(density) == np.ndarray:
                elements[key].density = density[position]

            #elements[key].density = density[position]*elements[key].stoichiometry/molar_mass  # sets density  (g/cm^3)
            elements[key].thickness = thickness[position]  # sets thickness  (Angstrom)
            elements[key].roughness = roughness[position]  # Order of Angstrom
            elements[key].linked_roughness = linked_roughness[position]  # sets the linked roughness for each element
            elements[key].molar_mass = molar_mass  # Molar mass of perovskite material
            elements[key].position = position

            # sets scattering factors for scattering factors that are different than their element symbol
            if sf is not None:
                my_sf = sf[position]
                if my_sf is not None:
                    elements[key].scattering_factor = my_sf
            position = position + 1
        self.structure[num_layer] = elements  # sets the layer with the appropriate slab properties


    def polymorphous(self, lay, ele, polymorph, poly_ratio, sf=''):
        """
        Purpose: Allows the user to set polymorphous elements in their material
        :param lay: The layer the polymorphous element is found (integer value)
        :param ele: The symbol of the polymorphous element in the selected layer (string)
        :param polymorph: List of the polymorph 'names' (e.g. ['Mn2+', 'Mn3+'])
        :param poly_ratio: A list of ratios that determine how much of the polymorph makes up the total element
                            - polymorph = ['Mn2+', 'Mn3+'], poly_ratio = [0.2,0.8]
                                - Density(Mn2+) = Density(Mn)*0.2
                                - Density(Mn3+) = Density(Mn)*0.8
                            - Note that the sum of all elements in poly_ratio must equal to 1
        :param sf: List of the scattering factors. Later version will implement scattering function.
        """

        # ------------------------------------------- Error Checks --------------------------------------------------- #

        # Checking that layer exists
        if lay > len(self.structure) - 1:
            raise RuntimeError( 'Layer ' + str(lay) + ' selected, but maximum layer is ' + str(len(self.structure)-1))

        # Checking that the element is found in the material
        if ele not in self.myelements:
            raise RuntimeError('Layer ' + str(lay) + ': Element '+ ele + ' not found in any layer.')

        # Checking that the element is found in the desired layer
        if ele not in list(self.structure[lay].keys()):
            raise RuntimeError('Layer ' + str(lay) + ': Element ' + ele + ' not found in layer ' + str(lay))

        # Checking polymorphous type
        if type(polymorph) != list and type(polymorph) != np.ndarray:
            raise TypeError('Layer ' + str(lay) + ': Input of polymorphous variable must be of type list')


        # Checks that length of list is greater than 2
        if len(polymorph) < 2:
            raise RuntimeError('Layer ' + str(lay) + ': Input more than one different version of the material. Length of list should be greater than two.')

        # Checks that poly_ratio and polymorphous have same length
        if len(poly_ratio) != len(polymorph):
            raise RuntimeError('Layer ' + str(lay) + ': variables poly_ratio and polymorph must be same length.')


        # Chekcs that all poly_ratio types are floating points
        for x in poly_ratio:
            if type(x) != int and type(x) != float:
                raise TypeError('Layer ' +str(lay)+ ': variable poly_ratio can only be type float or interger.')


        if abs(1 - sum(poly_ratio)) > 1e-3:
            print(lay, ':', poly_ratio)
            print(sum(poly_ratio))
            warnings.warn('Addition of all values in poly_ratio should add up to 1.')
        # ------------------------------------------------------------------------------------------------------------ #

        # Finds and sets the polymorphous elements to the right element class
        poly_keys = list(self.poly_elements.keys())
        if ele not in poly_keys:
            if len(polymorph) < 2: # Case where only one polymorph is entered
                self.poly_elements[ele] = [polymorph]
            else: # Multiple entries
                self.poly_elements[ele] = polymorph




        # Sets polymorph values to appropriate element in the correct layer
        self.structure[lay][ele].polymorph = polymorph
        self.structure[lay][ele].poly_ratio = np.array(poly_ratio)

        # Sets scattering factors
        if sf is str:
            self.structure[lay][ele].scattering_factor = polymorph

        else:
            self.structure[lay][ele].scattering_factor = sf

    def magnetization(self, lay, identifier, density, sf, gamma=90, phi=90):

        """
        Purpose: Set magnetization properties
        :param lay: Layer the magnetic material is found
        :param identifier: The symbol of the magnetic element
        :param density: The density of the magnetism
        :param sf: The scattering factor you want to use for the magnetic component
        :param gamma: Azimuthal angle (degrees)
        :param phi: Polar angle (degrees)
        """

        # ---------------------------------------- Error Check ------------------------------------------------------- #
        # Checking that layer exists
        if lay > len(self.structure) - 1:
            raise RuntimeError('Layer ' + str(lay) + ' selected, but maximum layer is ' + str(len(self.structure) - 1))

        self.layer_magnetized[lay] = True
        # Checks to make sure that the variable identifier is either a list or string
        if type(identifier) != str and type(identifier) != list and type(identifier) != np.ndarray:
            raise TypeError('Variable identifier must be a list or string.')

            # Checks to make sure that the variable identifier is either a list or string
        if type(sf) != str and type(sf) != list and type(sf) != np.ndarray:
            raise TypeError('Variable identifier must be a list or string.')

        # Checks to make sure that the variable identifier is either a list or string
        if type(density) != int and type(density) != list and type(density) != float:
            raise TypeError('Density must be entered as a list, integer, or float value')

        # Checks to make sure density is either float or integer value
        if type(density) == list:
            if len(density) != len(identifier) and len(identifier) != len(sf):
                raise RuntimeError('Variables density, identifier, and sf must contain same number of elements.')
            check = set([type(x) for x in density])
            for k in check:
                if k != int and k != float:
                    raise TypeError('Values in density list must be either integers or floating point values')
            for idx, dens in enumerate(density):
                if dens<0:
                    density[idx] = abs(dens)

        # Checks that all identifiers are strings
        if type(identifier) == list:

            if type(identifier[0]) is not str and type(identifier[0]) is not np.str_:
                raise TypeError('All list elements must be strings')
            if len(set([type(x) for x in identifier])) > 1:
                raise TypeError('All list elements must be strings')
            if len(identifier) != len(set(identifier)):
                raise RuntimeError('An element is repeated in identifier list.')


        if type(sf) == list:
            if type(sf[0]) != str:
                raise TypeError('All list elements must be strings')
            if len(set([type(x) for x in sf])) > 1:
                raise TypeError('All list elements must be strings')

        if type(gamma) == list and type(phi)==list:
            if len(gamma) != len(phi):
                raise TypeError('Gamma and Phi must be of same type')
            fg = gamma[0]
            fp = phi[0]
            for g in gamma:
                if g != fg:
                    raise SyntaxError('Current implementation required all elements to have same magnetization direction')
            for p in phi:
                if p != fp:
                    raise SyntaxError('Current implementation required all elements to have same magnetization direction')

        if type(phi) != list and type(gamma) == list:
            raise TypeError('Gamma and Phi must be of same type')

        if type(phi) == list and type(gamma) != list:
            raise TypeError('Gamma and Phi must be of same type')


        # ------------------------------------------------------------------------------------------------------------ #
        my_idx = [i for i in range(len(sf)) if sf[i] != '' and sf[i] != '0' and sf[i] != 0]
        layer = self.structure[lay]
        poly_start = True
        if type(identifier) == list or type(identifier) == np.ndarray:
            for key in list(layer.keys()):
                for idx in range(len(identifier)):
                    if len(layer[key].polymorph) != 0:  # polymorph case
                        if identifier[idx] in layer[key].polymorph:  # checks if identifier in selected polymorph
                            # pre-initializing for polymorph case
                            layer[key].polymorph = list(layer[key].polymorph)
                            poly_idx = layer[key].polymorph.index(identifier[idx])  # determine index of poly

                            # initialization
                            if poly_start:  # first polymorph appearance
                                self.structure[lay][key].mag_scattering_factor = ['' for i in range(len(layer[key].polymorph))]
                                self.structure[lay][key].mag_density = np.zeros(len(layer[key].polymorph))
                                #self.mag_elements[key] = [0 for i in range(len(layer[key].polymorph))]


                                self.mag_elements[key] = []
                                poly_start = False
                                # gamma and phi entered as multiple arrays
                                if type(gamma) == list and type(phi) == list:
                                    self.structure[lay][key].gamma = np.zeros(len(layer[key].polymorph))
                                    self.structure[lay][key].phi = np.zeros(len(layer[key].polymorph))
                            if identifier[idx] not in self.mag_elements[key]:
                                self.mag_elements[key].append(identifier[idx])

                            self.structure[lay][key].mag_scattering_factor[poly_idx] = sf[idx]
                            self.structure[lay][key].mag_density[poly_idx] = density[idx]



                            if type(gamma) == list and type(phi) == list:
                                self.structure[lay][key].gamma = gamma[idx]
                                self.structure[lay][key].phi = phi[idx]
                            else:
                                self.structure[lay][key].gamma = gamma
                                self.structure[lay][key].phi = phi
                    elif key == identifier[idx]:
                        self.structure[lay][identifier[idx]].mag_scattering_factor = [sf[idx]]
                        self.structure[lay][identifier[idx]].mag_density = np.array([density[idx]])

                        self.mag_elements[key] = [key]
                        # case where we want multiple different magnetization angles
                        # this implementation is for futur versions of code
                        if type(gamma) == list and type(phi) == list:
                            self.structure[lay][key].gamma = gamma[idx]
                            self.structure[lay][key].phi = phi[idx]
                        else:
                            self.structure[lay][key].gamma = gamma
                            self.structure[lay][key].phi = phi




    def error_function(self, t, rough, offset, upward):
        """
        Purpose: Computes the roughness using the error function
        :param t: Thickness list
        :param rough: Roughness (angstrom)
        :param offset: Thickness offset (angstrom)
        :param upward: Boolean used to determine if beginning of element or end of element
                        True - Beginning of element, error function slopes upward
                        False - End of element, error function slopes down
        :return: List representing error function
        """

        if rough == 0:  # Use heaviside function for zero roughness
            if upward:
                val = np.heaviside(t - offset, 1) * 2 - 1
            else:
                val = np.heaviside(offset-t, 1) * 2 - 1

        elif upward: # Upward error function
            val = erf((t-offset)/rough/np.sqrt(2))
        else:  # Downward error function
            val = erf((offset-t) / rough / np.sqrt(2))

        return val

    def density_profile(self, step=0.1):
        """
        Purpose: Creates the density profile based on the slab properties
        :return: thickness - thickness array in angstrom
                 density - structural density array in mol/cm^3
                 mag_density - magnetic density array in mol/cm^3
        """

        next_density = 0  # initialization required for algorithm
        n = len(self.structure)  # number of layers
        thickness = np.array([])  # thickness array
        density_struct = {k: np.array([]) for k in self.myelements}  # hold structure density
        density_poly = {k: dict() for k in list(self.poly_elements.keys())}  # hold polymorph elements
        #density_mag = {k: dict() for k in list(self.mag_elements.keys())}  # hold magnetic elements
        density_mag = {k: dict() for k in list(self.find_sf[1].keys())}

        # Pre-initialized density_poly array
        for ele in list(self.poly_elements.keys()):
            density_struct.pop(ele)
            density_poly[ele] = {k: np.array([]) for k in self.poly_elements[ele]}

        # Pre-initializes density_mag array

        for ele in list(self.mag_elements.keys()):
            density_mag[ele] = {k: np.array([]) for k in self.mag_elements[ele]}



        struct_keys = list(density_struct.keys()) # retrieves structure keys
        poly_keys = list(self.poly_elements.keys())  # retrieves polymorphous keys
        mag_keys = list(self.mag_elements.keys())  # retrieves magnetic keys

        # for an arbitrary size of elements
        # note that this relies that all layers in the slab definition have the same number of elements!
        num_ele = len(list(self.structure[0].keys()))  # assumes all layers have the same number of elements
        transition = [[0] for i in range(num_ele)]  # assumes same number elements through entire sample
        thick_array = np.array([0.0 for i in range(len(list(self.structure[0].keys())))])

        for layer in range(1,n):  # loop over all layers
            for i in range(num_ele):  # loop over all elements
                val = transition[i][layer-1] + list(self.structure[layer].values())[i].thickness

                transition[i].append(val)
                thick_array[i] = val

        thick = max(thick_array)

        thickness = np.arange(-25,thick+15+step, step) # Initializes the thickness array

        # Loop through elements in sample
        for ele in self.myelements:

            # structural elements (none polymorphous or magnetic)
            if ele in struct_keys:

                density_struct[ele] = np.zeros(len(thickness))  # sets array same length as thickness array full of zeros

                # Initializations so python is happy
                current_density = 0
                next_density = 0
                position = 0
                sigma = 0
                position = 0
                offset = 0
                offset_list = []

                for layer in range(n):
                    if ele in list(self.structure[layer].keys()):
                        position = self.structure[layer][ele].position  # position of element
                        offset_list = transition[position]  # offset for new implementation
                # Loops through all layers

                for layer in range(n):
                    offset = offset_list[layer]
                    #offset = transition[layer]

                    # Element is found in the current layer (ignore linked roughness)
                    if ele in list(self.structure[layer].keys()):

                        # saves scattering factor to be used in computation
                        if ele not in self.find_sf[0]:
                            self.find_sf[0][ele] = self.structure[layer][ele].scattering_factor
                            name = self.structure[layer][ele].scattering_factor
                            if name not in self.eShift.keys():
                                self.eShift[name] = 0
                                self.ff_scale[name] = 1


                        sigma = self.structure[layer][ele].roughness  # roughness parameterization

                        current_density = self.structure[layer][ele].density  # current density
                        if layer == n - 1:  # Last layer
                            next_density = 0  # density of element in next layer
                        elif ele in list(self.structure[layer+1].keys()):  # element in next layer
                            next_density = self.structure[layer+1][ele].density
                        else:  # element not in the next layer
                            next_density = 0

                        begin = 0
                        if layer == 0:
                            begin = 1

                        const = (next_density - current_density) / 2
                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        density_struct[ele] = density_struct[ele] + const*erf_func + begin*current_density

                    # Element is not found in the current layer (must take care of linked roughness)
                    else:
                        current_density = 0
                        if layer == n - 1:  # Last layer
                            next_density = current_density
                            sigma = 0
                        elif ele in list(self.structure[layer+1].keys()):  # element in next layer
                            if type(self.structure[layer+1][ele].linked_roughness) is float or type(self.structure[layer+1][ele].linked_roughness) is int:  # roughness is NOT linked to the previous site
                                sigma = self.structure[layer+1][ele].linked_roughness
                            else:  # roughness is linked to the previous site
                                position = self.structure[layer + 1][ele].position  # position of element
                                previous_element = list(self.structure[layer].keys())[position]
                                sigma = self.structure[layer][previous_element].roughness

                            next_density = self.structure[layer+1][ele].density # next layer density
                        else:
                            next_density = 0
                            sigma = 0

                        const = (next_density-current_density)/2
                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        density_struct[ele] = density_struct[ele] + const * erf_func

                density_struct[ele][density_struct[ele] < 0] = 0
            # Polymorph element
            pn=0
            if ele in poly_keys:

                # initialization of polymorph density dictionary
                layer = 0
                not_found = True
                while not_found or layer<=n-1:
                    if ele in list(self.structure[layer].keys()):
                        if len(list(self.structure[layer][ele].polymorph)) == 0:
                            raise SyntaxError('Polymorph not defined for ' + str(ele) + ' in layer ' + str(layer))

                        pn = len(list(self.structure[layer][ele].polymorph))  # number of polymorphs for element
                        density_poly[ele] = {k: np.zeros(len(thickness)) for k in list(self.structure[layer][ele].polymorph)}
                        not_found = False
                    layer = layer + 1

                if layer>n:  # element not found in any layers
                    raise RuntimeError(ele + ' defined as a polymorph, but not found in sample.')



            # Polymorphous elements
            if ele in poly_keys:
                for layer in range(n):
                    if ele in list(self.structure[layer].keys()):
                        position = self.structure[layer][ele].position  # position of element
                        offset_list = transition[position]

                for layer in range(n): # loops through all layers
                    offset = offset_list[layer]  # offset for new implementation
                    #offset = transition[layer]

                    # Element found in current layer
                    if ele in list(self.structure[layer].keys()):

                        position = self.structure[layer][ele].position  # position of element
                        sigma = self.structure[layer][ele].roughness  # roughness parameterization
                        current_density = self.structure[layer][ele].density*self.structure[layer][ele].poly_ratio  # current density

                        if layer == n - 1:  # On last layer
                            next_density = np.zeros(pn)  # density of element in next layer
                        elif ele in list(self.structure[layer + 1].keys()):  # element in next layer
                            next_density = self.structure[layer + 1][ele].density* self.structure[layer+1][ele].poly_ratio
                        else:  # element not in the next layer
                            next_density = np.zeros(pn)

                        begin = 0
                        if layer == 0:
                            begin = 1

                        # Implement changes ---------------------------------------------------------------------------
                        erf_func = self.error_function(thickness, sigma, offset, True)+1
                        const = (next_density-current_density)/2

                        po = 0
                        for poly in list(density_poly[ele].keys()):

                            # saves scattering factor of polymorphs
                            if ele not in self.find_sf[0]:
                                self.find_sf[0][poly] = self.structure[layer][ele].scattering_factor[po]
                                name = self.structure[layer][ele].scattering_factor[po]
                                if name not in self.eShift.keys():
                                    self.eShift[name] = 0
                                    self.ff_scale[name] = 1

                            # Density normalization
                            density_poly[ele][poly] = density_poly[ele][poly] + (const[po]*erf_func + begin*current_density[po])

                            po = po + 1

                    else:  # Element not found in current layer


                        current_density = np.zeros(pn)
                        if layer == n - 1:  # Last layer
                            next_density = current_density
                            sigma = 0
                        elif ele in list(self.structure[layer + 1].keys()):
                            if type(self.structure[layer+1][ele].linked_roughness) is float or type(self.structure[layer+1][ele].linked_roughness) is int:
                                sigma = self.structure[layer+1][ele].linked_roughness
                            else:
                                position = self.structure[layer + 1][ele].position  # position of element
                                previous_element = list(self.structure[layer].keys())[position]
                                sigma = self.structure[layer][previous_element].roughness

                            next_density = self.structure[layer + 1][ele].density * self.structure[layer+1][ele].poly_ratio  # next layer density

                        else:
                            next_density = np.zeros(pn)
                            sigma = 0

                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        const = (next_density - current_density) / 2
                        # Loops through all the polymorphs of the selected element
                        po = 0
                        for poly in list(density_poly[ele].keys()):
                            # Density normalization
                            density_poly[ele][poly] = density_poly[ele][poly] + const[po] * erf_func
                            density_poly[ele][poly][density_poly[ele][poly]<0] = 0
                            po = po + 1



            pm=0
            if ele in mag_keys:

                # initialization of magnetization density dictionary
                layer = 0
                not_found = True
                while not_found or layer <= n - 1:
                    if ele in list(self.structure[layer].keys()):
                        density_mag[ele] = {k: np.zeros(len(thickness)) for k in list(self.mag_elements[ele])}
                        if len(self.structure[layer][ele].mag_density) == 0:
                            raise SyntaxError('Magnetization not defined for ' + str(ele) + ' in layer ' + str(layer))
                        pm = len(self.structure[layer][ele].mag_density)
                        not_found = False
                    layer = layer + 1

                if layer > n:  # element not found in any layers
                    raise RuntimeError(ele + ' defined as a polymorph, but not found in sample.')

            # Magnetic elements
            if ele in mag_keys:
                for layer in range(n):
                    if ele in list(self.structure[layer].keys()):
                        position = self.structure[layer][ele].position  # position of element
                        offset_list = transition[position]

                for layer in range(n):  # loops through all layers
                    #offset = transition[layer]
                    offset = offset_list[layer]  # offset for new implementation

                    # Element found in current layer
                    if ele in list(self.structure[layer].keys()):
                        position = self.structure[layer][ele].position  # position of element
                        sigma = self.structure[layer][ele].roughness  # roughness parameterization
                        #current_density = self.structure[layer][ele].stoichiometry * np.array(self.structure[layer][ele].mag_density) * np.array(self.structure[layer][ele].density) / self.structure[layer][ele].molar_mass  # current density
                        current_density = np.array(self.structure[layer][ele].mag_density)
                        if layer == n - 1:  # Last layer
                            next_density = np.zeros(pm)  # density of element in next layer
                        elif ele in list(self.structure[layer + 1].keys()):  # element in next layer
                            #next_density = self.structure[layer + 1][ele].stoichiometry * np.array(self.structure[layer + 1][ele].mag_density) * np.array(self.structure[layer + 1][ele].density) / self.structure[layer + 1][ele].molar_mass
                            next_density = np.array(self.structure[layer + 1][ele].mag_density)
                        else:  # element not in the next layer
                            next_density = np.zeros(pm)

                        begin = 0
                        if layer == 0:
                            begin = 1

                        # Implement changes ---------------------------------------------------------------------------
                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        const = (next_density - current_density) / 2

                        ma = 0
                        for mag in list(density_mag[ele].keys()):
                            if len(self.structure[layer][ele].polymorph) != 0:
                                ma = list(self.structure[layer][ele].polymorph).index(mag)
                            # finds magnetic scattering factors
                            if mag not in self.find_sf[1]:
                                mag_sf = self.structure[layer][ele].mag_scattering_factor[ma]
                                my_check = ['',0,'0']

                                if mag_sf not in my_check:
                                    self.find_sf[1][mag] = self.structure[layer][ele].mag_scattering_factor[ma]
                                    name = self.structure[layer][ele].mag_scattering_factor[ma]


                            # Density normalization

                            density_mag[ele][mag] = density_mag[ele][mag] + (const[ma] * erf_func + begin * current_density[ma])
                            #ma = ma + 1

                    else:  # Element not found in current layer

                        current_density = np.zeros(pm)
                        if layer == n - 1:  # Last layer
                            next_density = current_density
                            sigma = 0
                        elif ele in list(self.structure[layer + 1].keys()):
                            if type(self.structure[layer+1][ele].linked_roughness) is float or type(self.structure[layer+1][ele].linked_roughness) is int:
                                sigma = self.structure[layer+1][ele].linked_roughness
                            else:
                                position = self.structure[layer + 1][ele].position  # position of element
                                previous_element = list(self.structure[layer].keys())[position]
                                sigma = self.structure[layer][previous_element].roughness
                            #next_density = self.structure[layer + 1][ele].stoichiometry * np.array(self.structure[layer + 1][ele].mag_density)*np.array(self.structure[layer + 1][ele].density) / self.structure[layer + 1][ele].molar_mass  # next layer density
                            next_density = np.array(self.structure[layer + 1][ele].mag_density)
                        else:
                            next_density = current_density
                            sigma = 0

                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        const = (next_density - current_density) / 2
                        # Loops through all the polymorphs of the selected element
                        ma = 0
                        for mag in list(density_mag[ele].keys()):

                            # Density normalization
                            density_mag[ele][mag] = density_mag[ele][mag] + const[ma] * erf_func
                            density_mag[ele][mag][density_mag[ele][mag]<0] = 0
                            ma = ma + 1


        # Create single dictionary to use (structural and polymorphs)
        density = density_struct
        for ele in list(density_poly.keys()):
            for poly in list(density_poly[ele].keys()):
                density[poly] = density_poly[ele][poly]

        # Create magnetic dictionary
        density_magnetic = dict()

        for ele in list(density_mag.keys()):
            for mag in list(density_mag[ele].keys()):
                if mag in list(self.find_sf[1].keys()):
                    density_magnetic[mag] = density_mag[ele][mag]

        self.transition = transition  # keeps track of transition
        return thickness, density, density_magnetic

    def _set_form_factors(self, element, ff, mag=False):
        """
        Purpose: Set the form factor
        :param element: Symbol or identifier of the form factor to be set
        :param ff: form factor value
        :param mag: Boolean to determine if magnetic form factor or not
        """
        if not mag:
            self.find_sf[0][element] = ff
        else:
            self.find_sf[1][element] = ff

    def plot_density_profile(self, fig=1, save=False, dir='Plot_Scans'):
        thickness, density, density_magnetic = self.density_profile()
        val = list(density.values())
        mag_val = list(density_magnetic.values())
        check = []
        for key in list(density.keys()):
            if key[-1].isdigit():
                check.append(True)
            else:
                check.append(False)

        plt.figure(fig)
        for idx in range(len(val)):
            if check[idx]:
                plt.plot(thickness, val[idx], ':')
            else:
                plt.plot(thickness, val[idx])

        for idx in range(len(mag_val)):
            plt.plot(thickness, -mag_val[idx], '--')

        center = np.zeros(len(thickness))
        plt.plot(thickness, center, 'k-.', linewidth=2)
        my_legend = list(density.keys())

        for key in list(density_magnetic.keys()):
            my_legend.append('Mag: ' + key)

        plt.legend(my_legend)
        #plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02, 0.5))
        plt.xlabel('Thickness (Angstrom)')
        plt.ylabel('Density (mol/cm^3)')

        if save:
            saveto = dir + '/Density_Profile.png'
            plt.savefig(saveto)



    def reflectivity(self, E, qz, precision=1e-6,s_min = 0.1, bShift=0,sFactor=1, sf_dict={}):

        """
        Purpose: Calculate reflectivity for constant energy using Pythonreflectivity
        :param E: Energy of reflectivity scan (eV)
        :param qi: Starting momentum transfer (A^{-1}) and related to small grazing angle
        :param qf: Ending momentum transfer (A^{-1}) and related to large grazing angle
        :param precision: Precision value for adaptive layer segmentation
        :param s_min: Minimum step size (None indicates using default value)
        :param bShift: float containing the background shift value
        :param sFactor: float containing the scaling factor value
        :return:
            qz - numpy array containing the momentum transfer
            R - dictionary for simulated reflectivity for the different types of x-ray polarizations

        """

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s
        wavelength = h * c / (E * 1e-10)  # wavelength of incoming x-ray

        # computes density profile based on the defined model (depth-dependent concentration)
        thickness, density, density_magnetic = self.density_profile(step=s_min)

        sf = dict()  # scattering factors of non-magnetic components
        sfm = dict()  # scattering factors of magnetic components

        #print(self.find_sf[1])
        if len(sf_dict) == 0:
            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])  # retrieve the energy shift of each scattering factor
                scale = float(self.ff_scale[self.find_sf[0][e]])  # retrieve scaling factor of each scattering factor
                sf[e] = find_form_factor(self.find_sf[0][e], E+dE, False)*scale  # find the scattering factor at energy E + dE
            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em],E + dE,True)*scale
        else:
            # Non-Magnetic Scattering Factor - no need to access original
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])  # retrieve the energy shift of each scattering factor
                scale = float(self.ff_scale[self.find_sf[0][e]])  # retrieve scaling factor of each scattering factor
                sf[e] = find_ff(self.find_sf[0][e],E+dE,sf_dict)

            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em], E + dE, True) * scale


        delta, beta = index_of_refraction(density, sf, E)  # calculates depth-dependent refractive index components
        delta_m, beta_m = magnetic_optical_constant(density_magnetic, sfm, E)   # calculates depth-dependent magnetic components
        if type(delta_m) != list and type(delta_m) != np.ndarray:
            delta_m = np.zeros(len(delta))
        if type(beta_m) != list and type(beta_m) != np.ndarray:
            beta_m = np.zeros(len(beta))

        # definition of magneto-optical constant as described in Lott Dieter Thesis
        n = 1 + np.vectorize(complex)(-delta, beta)  # complex index of refraction
        epsilon = n**2  # dielectric constant computation

        # magneto-optical constant as defined in Lott Dieter Thesis
        Q = np.vectorize(complex)(beta_m, delta_m)
        epsilon_mag = Q*epsilon*2*(-1)

        my_slabs = ALS(epsilon.real, epsilon.imag, Q.real, Q.imag, precision)  # performs the adaptive layer segmentation using Numba

        my_slabs = my_slabs.astype(int)  # sets all values in my_slab to integers

        my_slabs = my_slabs[1:]  # removes first element as it is not needed for structure generation


        m = len(my_slabs)  # number of slabs

        A =pr.Generate_structure(m)  # initializes Pythonreflectivity object class

        m_j=0  # previous slab
        idx = 0  # keeps track of current layer
        layer = 0
        gamma = 90  # pre-initialize magnetization direction
        phi = 90


        # This section will need to be altered for newer version of PythonReflectivity as the newer version uses
        # the definition of chi instead of epsilon
        for m_i in my_slabs:
            d = thickness[m_i] - thickness[m_j]  # computes thickness of slab

            eps = epsilon[m_j]  # non-magnetic dielectric constant

            eps_mag = epsilon_mag[m_j]  # magnetic dielectric constant

            # Determines the magnetization direction of the first layer
            if layer == 0:
                if self.layer_magnetized[0]:
                    for ele in self.structure[layer].keys():
                        if len(self.structure[0][ele].mag_scattering_factor) != 0:
                            gamma = self.structure[0][ele].gamma
                            phi = self.structure[0][ele].phi


            # Determines the magnetization direction of the other layers
            temp_list = [self.transition[i][layer] for i in range(len(self.transition))]
            transition = max(temp_list)

            # makes sure we are properly defining the magnetization direction as desired
            if transition<=thickness[m_j] and layer<len(self.transition[0])-1:
                layer = layer + 1
                if self.layer_magnetized[layer]:
                    for ele in self.structure[layer].keys():
                        if self.structure[layer][ele].mag_scattering_factor is not None:
                            gamma = self.structure[layer][ele].gamma
                            phi = self.structure[layer][ele].phi


            # sets the magnetization direction based on the input angles
            if self.layer_magnetized[layer]:
                if gamma == 90 and phi == 90:
                    A[idx].setmag('y')
                elif gamma == 0 and phi == 90:
                    A[idx].setmag('x')
                elif gamma == 0 and phi == 0:
                    A[idx].setmag('z')
                else:
                    raise ValueError('Values of Gamma and Phi can only be (90,90), (0,90), and (0,0)')

                A[idx].seteps([eps,eps,eps,eps_mag])   # sets the components of the dielectric tensor
            else:
                A[idx].seteps([eps,eps,eps,0])  # non-magnetic case


            if idx != 0:
                A[idx].setd(d)  # sets thickness of layer if and only if not substrate layer


            # move onto the next layer
            m_j = m_i
            idx = idx + 1


        Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
        Rtemp = pr.Reflectivity(A, Theta, wavelength, MagneticCutoff=1e-20)  # Computes the reflectivity
        R = dict()

        # creates the reflectivity dictionary applying the defined scaling factors and background shifts
        if len(Rtemp) == 2:
            R['S'] = Rtemp[0]*sFactor + bShift  # s-polarized light
            R['P'] = Rtemp[1]*sFactor + bShift   # p-polarized light
            R['AL'] = sFactor*(Rtemp[0]-Rtemp[1])/(sFactor*(Rtemp[0]+Rtemp[1])+2*bShift)  # Asymmetry linear polarized
            R['LC'] = Rtemp[0]*sFactor + bShift  # Left circular
            R['RC'] = Rtemp[0]*sFactor + bShift  #  right circular
            R['AC'] = np.zeros(len(Rtemp[0]))  # Asymmetry circular polarized (XMCD)
        elif len(Rtemp) == 4:
            R['S'] = Rtemp[0]*sFactor + bShift
            R['P'] = Rtemp[1]*sFactor + bShift
            R['AL'] = sFactor*(Rtemp[0]-Rtemp[1])/(sFactor*(Rtemp[0]+Rtemp[1])+2*bShift)
            R['LC'] = Rtemp[2]*sFactor + bShift
            R['RC'] = Rtemp[3]*sFactor + bShift
            R['AC'] = sFactor*(Rtemp[2]-Rtemp[3])/(sFactor*(Rtemp[2]+Rtemp[3])+2*bShift)

        else:
            raise TypeError('Error in reflectivity computation. Reflection array not expected size.')

        return qz, R

    def reflectivity_udkm(self, E, qz, precision=1e-6,s_min = 0.1, bShift=0,sFactor=1, sf_dict={}):

        """
        Purpose: Calculates  reflectivity for constant energy using ukdm1Dsim
        :param E: Energy of reflectivity scan (eV)
        :param qi: Starting momentum transfer (A^{-1}) and related to small grazing angle
        :param qf: Ending momentum transfer (A^{-1}) and related to large grazing angle
        :param precision: Precision value for adaptive layer segmentation
        :param s_min: Minimum step size (None indicates using default value)
        :param bShift: float containing the background shift value
        :param sFactor: float containing the scaling factor value
        :return:
            qz - numpy array containing the momentum transfer
            R - dictionary for simulated reflectivity for the different types of x-ray polarizations

        """

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s
        wavelength = h * c / (E * 1e-10)  # wavelength of incoming x-ray

        # computes density profile based on the defined model (depth-dependent concentration)
        thickness, density, density_magnetic = self.density_profile(step=s_min)

        sf = dict()  # scattering factors of non-magnetic components
        sfm = dict()  # scattering factors of magnetic components

        #print(self.find_sf[1])
        if len(sf_dict) == 0:
            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])  # retrieve the energy shift of each scattering factor
                scale = float(self.ff_scale[self.find_sf[0][e]])  # retrieve scaling factor of each scattering factor
                sf[e] = find_form_factor(self.find_sf[0][e], E+dE, False)*scale  # find the scattering factor at energy E + dE
            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em],E + dE,True)*scale
        else:
            # Non-Magnetic Scattering Factor - no need to access original
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])  # retrieve the energy shift of each scattering factor
                scale = float(self.ff_scale[self.find_sf[0][e]])  # retrieve scaling factor of each scattering factor
                sf[e] = find_ff(self.find_sf[0][e],E+dE,sf_dict)

            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em], E + dE, True) * scale


        delta, beta = index_of_refraction(density, sf, E)  # calculates depth-dependent refractive index components
        delta_m, beta_m = magnetic_optical_constant(density_magnetic, sfm, E)   # calculates depth-dependent magnetic components
        if type(delta_m) != list and type(delta_m) != np.ndarray:
            delta_m = np.zeros(len(delta))
        if type(beta_m) != list and type(beta_m) != np.ndarray:
            beta_m = np.zeros(len(beta))

        # definition of magneto-optical constant as described in Lott Dieter Thesis
        n = 1 + np.vectorize(complex)(-delta, beta)  # complex index of refraction
        epsilon = n**2  # dielectric constant computation

        # magneto-optical constant as defined in Lott Dieter Thesis
        Q = np.vectorize(complex)(beta_m, delta_m)
        epsilon_mag = Q*epsilon*2*(-1)

        my_slabs = ALS(epsilon.real, epsilon.imag, Q.real, Q.imag, precision)  # performs the adaptive layer segmentation using Numba

        my_slabs = my_slabs.astype(int)  # sets all values in my_slab to integers

        my_slabs = my_slabs[1:]  # removes first element as it is not needed for structure generation


        m = len(my_slabs)  # number of slabs

        R = dict()
        return qz, R

    def energy_scan(self, Theta, energy, precision=1e-11,s_min = 0.1, bShift=0, sFactor=1, sf_dict={}):
        """
        Purpose: Calculates reflectivity for constant grazing angle using Pythonreflectivity
        :param Theta: Grazing angle in degrees
        :param energy: List or numpy array containing the energies in the energy scan
        :param precision: parameter used in adaptive layer segmentation
        :param s_min: minimum slab slice
        :return: energy, R
                    - energy: exact energy array that user input
                    - R: dictionary containing the reflectivity calculation
                        'S' - s-polarization
                        'P' - p-polarized
                        'AL' - asymmetry linear
                        'RC' - right circular
                        'LC' - left circular
                        'AC' - asymmetry circular

        """
        # initializes the reflectivity array to be returned
        Elen = len(energy)
        R = {'S': np.zeros(Elen),
             'P': np.zeros(Elen),
             'AL': np.zeros(Elen),
             'LC': np.zeros(Elen),
             'RC': np.zeros(Elen),
             'AC': np.zeros(Elen)}

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s

        thickness, density, density_magnetic = self.density_profile(step=s_min)  # Computes the density profile
        # Magnetic Scattering Factor
        sfm = dict()
        sf = dict()

        if len(sf_dict) == 0:
            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])
                scale = float(self.ff_scale[self.find_sf[0][e]])
                sf[e] = find_form_factor(self.find_sf[0][e], energy + dE, False)*scale
            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em], energy + dE, True)*scale
        else:
            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])
                scale = float(self.ff_scale[self.find_sf[0][e]])
                sf[e] = find_ff(self.find_sf[0][e], energy + dE, sf_dict)
            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em], energy + dE, True) * scale


        d_len = len(thickness)
        delta, beta = IoR(density, sf, energy)  # gets absorptive and dispersive components of refractive index

        delta_m, beta_m = MOC(density_magnetic, sfm,energy, d_len)  # absorptive and dispersive components for magnetic components

        epsilon = 1 - 2*delta + 1j*beta*2  # dielectric constant

        # definition as described in Lott Dieter Thesis
        Q = beta_m + 1j*delta_m  # magneto-optical constant
        epsilon_mag = Q * epsilon *(-2)  # magneto-optical permittivity
        # retrieves the slabs at each energy using list comprehension
        all_slabs = [ALS(epsilon[E].real,epsilon_mag[E].imag, Q[E].real, Q[E].imag, precision=precision)[1:].astype(int) for E in range(len(energy))]

        # initializes the object for reflectivity computation using list comprehension

        Alist = [generate_structure(thickness, self.structure, all_slabs[s], epsilon[s], epsilon_mag[s], self.layer_magnetized, self.transition) for s in range(len(all_slabs))]
        wavelength = h * c / (energy * 1e-10)
        # reflectivity computation using list comprehension

        R = [energy_reflectivity(Alist[E],Theta, wavelength[E], R, int(E), backS=bShift, scaleF=sFactor) for E in range(len(all_slabs))]


        R = R[0]

        return energy, R

    def energy_scan_udkm(self, Theta, energy, precision=1e-11,s_min = 0.1, bShift=0, sFactor=1, sf_dict={}):
        """
        Purpose: Calculates reflectivity for constant grazing angle using udkm1Dsim
        :param Theta: Grazing angle in degrees
        :param energy: List or numpy array containing the energies in the energy scan
        :param precision: parameter used in adaptive layer segmentation
        :param s_min: minimum slab slice
        :return: energy, R
                    - energy: exact energy array that user input
                    - R: dictionary containing the reflectivity calculation
                        'S' - s-polarization
                        'P' - p-polarized
                        'AL' - asymmetry linear
                        'RC' - right circular
                        'LC' - left circular
                        'AC' - asymmetry circular

        """
        # initializes the reflectivity array to be returned
        Elen = len(energy)
        R = {'S': np.zeros(Elen),
             'P': np.zeros(Elen),
             'AL': np.zeros(Elen),
             'LC': np.zeros(Elen),
             'RC': np.zeros(Elen),
             'AC': np.zeros(Elen)}

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s

        thickness, density, density_magnetic = self.density_profile(step=s_min)  # Computes the density profile
        # Magnetic Scattering Factor
        sfm = dict()
        sf = dict()

        if len(sf_dict) == 0:
            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])
                scale = float(self.ff_scale[self.find_sf[0][e]])
                sf[e] = find_form_factor(self.find_sf[0][e], energy + dE, False)*scale
            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em], energy + dE, True)*scale
        else:
            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                dE = float(self.eShift[self.find_sf[0][e]])
                scale = float(self.ff_scale[self.find_sf[0][e]])
                sf[e] = find_ff(self.find_sf[0][e], energy + dE, sf_dict)
            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                dE = float(self.mag_eShift[self.find_sf[1][em]])
                scale = float(self.ffm_scale[self.find_sf[1][em]])
                sfm[em] = find_form_factor(self.find_sf[1][em], energy + dE, True) * scale


        d_len = len(thickness)
        delta, beta = IoR(density, sf, energy)  # gets absorptive and dispersive components of refractive index

        delta_m, beta_m = MOC(density_magnetic, sfm,energy, d_len)  # absorptive and dispersive components for magnetic components

        epsilon = 1 - 2*delta + 1j*beta*2  # dielectric constant

        # definition as described in Lott Dieter Thesis
        Q = beta_m + 1j*delta_m  # magneto-optical constant
        epsilon_mag = Q * epsilon *(-2)  # magneto-optical permittivity
        # retrieves the slabs at each energy using list comprehension
        all_slabs = [ALS(epsilon[E].real,epsilon_mag[E].imag, Q[E].real, Q[E].imag, precision=precision)[1:].astype(int) for E in range(len(energy))]

        return energy, R
    def energy_shift(self):
        """
        Purpose: Initialize the energy shift and form factor scaling for the GUI
        """

        self.density_profile()
        #self.eShift = dict()
        #self.mag_eShift = dict()

        key_delete = []
        mag_key_delete = []
        for e in self.find_sf[0].keys():
            if self.find_sf[0][e] == '' or self.find_sf[0][e] == 0 or self.find_sf[0][e] == '0':
                key_delete.append(e)

        for em in self.find_sf[1].keys():
            if self.find_sf[1][em] == '' or self.find_sf[1][em] == 0 or self.find_sf[1][em] == '0':
                mag_key_delete.append(em)

        for key in key_delete:
            del self.find_sf[0][key]
        for key in mag_key_delete:
            del  self.find_sf[1][key]

        for ele in self.mag_elements.keys():
            for key in key_delete:
                if key in self.mag_elements[ele]:
                    self.mag_elements[ele].remove(key)



    def getRoughness(self, layer, identifier):
        """
        Purpose: Retrieve roughness of a specific layer and element
        :param layer: Layer of integer value
        :param identifier: String containing element symbol or dummy variable symbol
        :return: Roughness as a float value
        """
        # gets the roughness of a select element

        keys = list(self.structure[layer].keys())
        sigma = 2
        if identifier in keys:
            sigma = self.structure[layer][identifier].roughness

        return sigma



    def setRoughness(self, layer, identifier, sigma):
        """
        Purpose: Sets roughness of a specific layer and element
        :param layer: Layer of integer value
        :param identifier: String containing element symbol or dummy variable symbol
        :param sigma: Roughness value of type float
        """
        # sets the roughness of a certain element in a certain layer
        keys = list(self.structure[layer].keys())

        if identifier in keys:
            self.structure[layer][identifier].roughness = sigma
        elif identifier.upper() == 'ALL':
            for key in keys:
                self.structure[layer][key].roughness = sigma

    def getDensity(self, layer, identifier):
        """
        Purpose: Retrieve density of a specific layer and element
        :param layer: Layer of integer value
        :param identifier: String containing element symbol or dummy variable symbol
        :return: Density as a float value
        """

        keys = list(self.structure[layer].keys())
        density = 0.028
        if identifier in keys:
            density = self.structure[layer][identifier].density

        return density



    def setDensity(self, layer, identifier, density):
        """
        Purpose: Sets density of a specific layer and element
        :param layer: Layer of integer value
        :param identifier: String containing element symbol or dummy variable symbol
        :param density: Density as float value
        """
        # sets the roughness of a certain element in a certain layer
        keys = list(self.structure[layer].keys())

        if identifier in keys:
            self.structure[layer][identifier].density = density
        elif identifier.upper() == 'ALL':
            for key in keys:
                self.structure[layer][key].density = density



    def setVariationConstant(self, layer, symbol, identifier, val):
        """
        Purpose: Sets element variation of a certain layer to a constant value
        :param layer: Layer as an integer type
        :param symbol: Symbol of element or dummy variable
        :param identifier: Element variation identifier
        :param val: Constant value as float type
        :return:
        """
        # for now this will only work for 3 element variations
        keys = list(self.structure[layer].keys())

        # loops through elements in the layer until symbol is found
        if symbol in keys:
            if len(self.structure[layer][symbol].polymorph) == 3:
                # determines the index of the polymorph that is not identifer 1 or 2
                idx_no = [i for i in range(3) if self.structure[layer][symbol].polymorph[i] != identifier]  # indices of non-desired element variation
                idx = [i for i in range(3) if self.structure[layer][symbol].polymorph[i] == identifier]  # index of desired element variation

                my_total = sum(list(self.structure[layer][symbol].poly_ratio[idx_no] ))  # total value of non-desired element variation values

                # sets non-constant element variation ratio's appropriately
                for i in idx_no:
                    if my_total == 0:
                        r = 1/2
                        self.structure[layer][symbol].poly_ratio[i] = (1-val)*r
                    else:
                        r = self.structure[layer][symbol].poly_ratio[i]
                        self.structure[layer][symbol].poly_ratio[i] = (1 - val) * r / my_total

                self.structure[layer][symbol].poly_ratio[idx] = val  # sets constant ratio value

    def setMultiVarConstant(self, layer, symbol, identifier, value):
        """
        Purpose: Sets multiple element variation ratios constant. This is ideal when optimizing an element with more than four variations.
                 It should be noted that this function only works when two element variations are allowed to vary.
                  - len(identifier) = total number of element variations - 2
        :param layer: integer type that signals the layer
        :param symbol: Symbol of the element or dummy variable
        :param identifier: List of identifiers to hold constant
        :param value: List of the ratio values
        :return:
        """

        gamma = 1-sum(value)  # value of remaining element variations ratio sum
        keys = list(self.structure[layer].keys())  # finds all the elements in the selected layer

        # loops through elements in the layer until symbol is found
        if symbol in keys:
            if len(self.structure[layer][symbol].polymorph) - len(identifier) <= 2:
                # determines the index for the variations to hold constant

                idx = [i for i in range(len(self.structure[layer][symbol].polymorph)) if self.structure[layer][symbol].polymorph[i] not in identifier]

                beta = sum(self.structure[layer][symbol].poly_ratio[idx])  # sum of the polymorphs


                # sets non-constant element variation ratio's appropriately
                for i in idx:
                    if beta == 0:
                        r = 1/2
                        self.structure[layer][symbol].poly_ratio[i] = gamma * r
                    else:
                        r = self.structure[layer][symbol].poly_ratio[i]
                        self.structure[layer][symbol].poly_ratio[i] = gamma * r / beta

                # sets the constant values
                my_i = 0
                # loops through each keys in the polymorphs
                for i,key in enumerate(self.structure[layer][symbol].polymorph):
                    if key in identifier:  # if the key is in the identifiers than set the appropriate value
                        self.structure[layer][symbol].poly_ratio[i] = value[my_i]
                        my_i = my_i + 1



    def setRatio(self,layer, symbol, identifier1, identifier2, ratio):
        """
        Purpose: Sets the ratio between two element variation to be constant
        :param layer: Layer as integer type
        :param symbol: Symbol of element or dummy variable
        :param identifier1: Element variation identifier
        :param identifier2: Element variation identifier
        :param ratio: Constant ratio value as float type
        """
        # ratio = identifier1/identifier2
        # for now this will only work for 3 element variations
        keys = list(self.structure[layer].keys())

        # loops through all elements in layer until symbol is found
        if symbol in keys:
            if len(self.structure[layer][symbol].polymorph) == 3:
                # determines the index of the polymorph that is not identifer 1 or 2
                idx_no = [i for i in range(3) if self.structure[layer][symbol].polymorph[i] != identifier1 and
                       self.structure[layer][symbol].polymorph[i] != identifier2]

                val = self.structure[layer][symbol].poly_ratio[idx_no][0]

                idx1 = list(self.structure[layer][symbol].polymorph).index(identifier1)
                idx2 = list(self.structure[layer][symbol].polymorph).index(identifier2)

                self.structure[layer][symbol].poly_ratio[idx1] = (1-val)/(ratio+1)
                self.structure[layer][symbol].poly_ratio[idx2] = ratio*(1-val)/(ratio+1)


    def getThickness(self, layer, identifier):
        """
        Purpose: Retrieve thickness of a specific layer and element
        :param layer: Layer of integer value
        :param identifier: String containing element symbol or dummy variable symbol
        :return: Thickness as a float value
        """
        # gets the thickness of a certain element in a certain layer
        keys = list(self.structure[layer].keys())
        d = 2
        if identifier in keys:
            d = self.structure[layer][identifier].thickness
        elif identifier.upper() == 'ALL':
            d = self.structure[layer][keys[0]].thickness

        return d

    def getTotalThickness(self,start_layer,end_layer, identifier):
        """
        Purpose: Retrieves total thickness from start_layer to end_layer for specified element or dummy variable
        :param start_layer: Index of layer closest to substrate
        :param end_layer: Index of layer furthest from substrate
        :param identifier: String containing element symbol or dummy variable symbol
        :return: Total thickness as a float value
        """
        # gets the total thickness from start_layer to end_layer of a specific element (identifier)
        d = 0
        for i in range(start_layer, end_layer+1,1):
            keys = list(self.structure[i].keys())

            if identifier in keys:
                d = d + self.structure[i][identifier].thickness
            elif identifier.upper() == 'ALL':  # case where we want the total thickness of all elements
                d = d + self.structure[i][keys[0]].thickness

        return d

    def setThickness(self,layer, identifier, d):
        """
        Purpose: Sets thickness of a specific layer and element
        :param layer: Layer of integer value
        :param identifier: String containing element symbol or dummy variable symbol
        :return: Thickness as a float value
        """
        keys = list(self.structure[layer].keys())

        if identifier in keys:
            self.structure[layer][identifier].thickness = d
        elif identifier.upper() == 'ALL':
            for key in keys:
                self.structure[layer][key].thickness = d

    def setCombinedThickness(self, layer_start, layer_end, identifier, d):
        """
        Purpose: Sets total thickness from start_layer to end_layer for specified element or dummy variable
        :param start_layer: Index of layer closest to substrate
        :param end_layer: Index of layer furthest from substrate
        :param identifier: String containing element symbol or dummy variable symbol
        :param d: Total thickness value as float type
        """
        # makes sure that the combine thickness from layer_start to layer_end is constant (d)
        import copy
        dprime = 0
        for i in range(layer_start, layer_end + 1, 1):
            keys = list(self.structure[i].keys())
            if identifier in keys:
                dprime = dprime + self.structure[i][identifier].thickness
            elif identifier.upper() == 'ALL':
                dprime = dprime + self.structure[i][keys[0]].thickness


        for i in range(layer_start, layer_end+1, 1):
            keys = list(self.structure[i].keys())
            if identifier in keys:
                val = copy.deepcopy(self.structure[i][identifier].thickness)
                self.structure[i][identifier].thickness = val*d/dprime
            elif identifier.upper() == 'ALL':
                val = copy.deepcopy(self.structure[i][keys[0]].thickness)
                for key in keys:
                    self.structure[i][key].thickness = val * d / dprime

    def getMagDensity(self, layer,symbol,variation):
        """
        Purpose: Retrieves the magnetic density
        :param layer: Layer of integer type
        :param symbol: Symbol of element or dummy variable
        :param variation: Element variation identifier, if not an element variation set this to the same as symbol
        :return mag_density: Magnetic density as float type
        """
        # retrieve the magnetic density
        # if there is no element variation just input the same variable

        polymorph = self.structure[layer][symbol].polymorph
        if len(polymorph) != 0:
            idx = [i for i in range(len(polymorph)) if variation == polymorph[i]]
            mag_density = self.structure[layer][symbol].mag_density[idx]

            #if type(mag_density) is list or type(mag_density) is np.ndarray:
            #    mag_density = mag_density
        else:
            mag_density = self.structure[layer][symbol].mag_density[0]

        return mag_density

    def setMagDensity(self, layer, symbol, variation, density):
        """
        Purpose: Sets the magnetic density
        :param layer: Layer of integer type
        :param symbol: Symbol of element or dummy variable
        :param variation: Element variation identifier, if not an element variation set this to the same as symbol
        :param density: Magnetic density value of float type
        """
        # set the magnetic density
        polymorph = self.structure[layer][symbol].polymorph
        if len(polymorph) != 0:
            idx = [i for i in range(len(polymorph)) if variation == polymorph[i]]

            self.structure[layer][symbol].mag_density[idx] = float(density)

        else:
            self.structure[layer][symbol].mag_density[0] = float(density)


    def getEshift(self,ffName):
        """
        Purpose: Retrieves the energy shift of specified non-magnetic form factor
        :param ffName: form factor name without .ff extension
        :return: Form factor energy shift as float type
        """
        # retieve the energy shift
        return self.eShift[ffName]

    def setEshift(self, ffName, dE):
        """
        Purpose: Set the non-magnetic form factor energy shift
        :param ffName: Form factor name without .ff extension
        :param dE: Energy shift as a float type
        """
        self.eShift[ffName] = dE

    def getMagEshift(self,ffmName):
        """
        Purpose: Retrieves the energy shift of specified magnetic form factor
        :param ffmName: form factor name without .ffm extension
        :return: Magnetic form factor energy shift as float type
        """
        # retrieve the magnetic energy shift
        return self.mag_eShift[ffmName]

    def setMagEshift(self, ffmName, dE):
        """
        Purpose: Set magnetic form factor energy shift
        :param ffmName: Magnetic form factor name without .ffm extension
        :param dE: Energy shift as float type
        """
        # set the magnetic energy shift
        self.mag_eShift[ffmName] = dE



if __name__ == "__main__":
    print('Nope')










