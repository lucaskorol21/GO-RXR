import numpy as np
import matplotlib.pyplot as plt
import math
from collections import OrderedDict
from scipy.special import erf
from scipy.integrate import simpson
import scipy.optimize as optimize
import warnings
from KK_And_Merge import *
from material_model import *
import Pythonreflectivity as pr
from tabulate import tabulate
from scipy import interpolate
import multiprocessing
from functools import partial
from time import *
from multiprocessing.pool import ThreadPool
from numba import *


def zero_to_one(func):
    """
    Purpose: Spans the function over a range from 0 to 1 for layer segmentation
    :param func: Input function, in our case it is the dielectric constant function
    :return: The transformed function
    """
    func_min = min(func)  # maximum value in the function
    func_max = max(func)  # minimum value in the function
    amplitude = func_max - func_min  # computes the amplitude from min to max
    if amplitude==0:
        return func
    else:
        return (func-func_min)/amplitude

def total_variation(func):
    """
    Purpose: Calculates the total variation of the input function
    :param func: Input function, in our case the dielectric constant
    :return: Total variation
    """

    tv = np.sum(np.abs(np.diff(np.array(func))))

    return tv


def slice_diff(thickness, eps, idx_a, idx_b, precision, n ):
    """
    Purpose: Determines the steps size to take for the layer segmentation
    :param d: Total thickness
    :param eps: dielectric constant array
    :param idx_a: current index
    :param idx_b: next index
    :param precision: precision value
    :param n: number of elements in dielectric constant array
    :return:
    """

    f1 = eps[idx_a]  # current dielectric constant value
    f2 = eps[idx_b]  # dielectric constant value of next step
    delta = abs(f2-f1)  # difference between dielectric constant values


    # Determines if layer segmentation is too large
    while(delta<precision and not(idx_b>=n-1)):
        idx_b = idx_b + 1  # increases step size until too large of layer segmentation take
        f2 = eps[idx_b]  # retrieves dielectric constant
        delta = abs(f2-f1)  # computes difference


    return idx_b


def slice_diff_new(tck_eps, d_current,s_min, d_next, precision):
    """
    Purpose: Determines the steps size to take for the layer segmentation
    :param d: Total thickness
    :param eps: dielectric constant array
    :param idx_a: current index
    :param idx_b: next index
    :param precision: precision value
    :param n: number of elements in dielectric constant array
    :return:
    """
    delta_d = s_min*0.1
    f1 = interpolate.splev(d_current, tck_eps)  # current dielectric constant value
    f2 = interpolate.splev(d_next,tck_eps)# dielectric constant value of next step
    delta = abs(f2 - f1)  # difference between dielectric constant values

    # Determines if layer segmentation is too large
    while (delta < precision):
        d_next = d_next + delta_d
        f2 = interpolate.splev(d_next, tck_eps)  # retrieves dielectric constant
        delta = abs(f2 - f1)  # computes difference
        if delta < precision and abs(d_next-d_current)<s_min:
            d_next = d_current + s_min


    return d_next


def layer_segmentation(thickness, epsilon, epsilon_mag, precision=1e-6):
    """
    Purpose: Performs the layer segmentation over entire interval based on total variation
    :param thickness: Thickness array of sample (angstrom)
    :param epsilon: Structural dielectric constant array
    :param epsilon_mag: Magnetic dielectric constant array
    :param precision: Precision value
    :return: Indices for slicing
    """
    idx_a = 0  # current index
    idx_b = 1  # index of next slice
    n = len(epsilon)  # total number of elements in dielectric constant array
    my_slabs = list()
    d = thickness[-1]-thickness[0]  # total thickness of the sample
    delta_d = thickness[1]-thickness[0]  # thickness step
    # computes total variation for real and imaginary of structural and magnetic dielectric constants

    #max_1 = max(abs(np.diff(real(epsilon))))
    #max_2 = max(abs(np.diff(imag(epsilon))))
    #max_3 = max(abs(np.diff(real(epsilon_mag))))
    #max_4 = max(abs(np.diff(imag(epsilon_mag))))

    #p_1 = precision * max_1
    #p_2 = precision * max_2
    #p_3 = precision * max_3
    #p_4 = precision * max_4
    p_1 = precision
    p_3 = precision


    while (idx_b < n):

        # Computes the precision values based on the average variance over the interval
        idx_s_r = slice_diff(thickness, real(epsilon), idx_a, idx_b, p_1, n)  # structural real component
        #idx_s_i = slice_diff(thickness, imag(epsilon), idx_a, idx_b, p_2, n)  # structural imaginary component
        idx_m_r = slice_diff(thickness, real(epsilon_mag), idx_a, idx_b, p_3, n)  # magnetic real component
        #idx_m_i = slice_diff(thickness, imag(epsilon_mag), idx_a, idx_b, p_4, n)  # magnetic imaginary component

        #idx_b = min(idx_s_r, idx_s_i, idx_m_r, idx_m_i)  # use the smallest slice value
        idx_b = min(idx_s_r,  idx_m_r)  # use the smallest slice value
        my_slabs.append(idx_b)
        idx_a = idx_b  # step to next slab
        idx_b = idx_b + 1

    return my_slabs

@njit()
def ALS(epsilon, epsilon_mag, precision=1e-6):
    """
    Purpose: Perform the adaptive layer segmentation
    :param epsilon: numpy array of real values
    :param epsilon_mag: numpy array of real values
    :param precision: precision for slicing
    :return: my_slabs - contains indices for slicing
    """
    epsilon = epsilon/np.linalg.norm(epsilon)  # normalizes epsilon
    epsilon_mag = epsilon_mag/np.linalg.norm(epsilon)  # normalizes epsilon_mag
    idx_a = 0  # keeps track of surface of previous slab
    n = epsilon.size
    my_slabs = np.zeros(n) # pre-initialize the slab array

    dsSlab = 1
    for idx_b in range(1,n):

        # retrieves permittivity values
        f1 = epsilon[idx_a]
        f2 = epsilon[idx_b]
        f1m = epsilon_mag[idx_a]
        f2m = epsilon_mag[idx_b]

        delta = np.absolute(f2-f1)  # varitation of epsilon
        delta_m = np.absolute(f2m-f1m)  # variation of epsilon_mag

        # checks if variation is of minimum variation set by 'precision'
        if delta>precision or delta_m>precision:
            #my_slabs = np.append(my_slabs, idx_b)  # append slice
            my_slabs[dsSlab] = idx_b
            idx_a = idx_b  # change previous slice location
            dsSlab = dsSlab + 1

    my_slabs = my_slabs[:dsSlab]
    return my_slabs


class element:
    def __init__(self, name, stoichiometry):
        """
        Purpose: Used in class slab. Used to keep track of elemental properties of each layer
        :param name: Element Symbol
        :param stoichiometry: The stoichiometric relation of the chemical formula as input into the class slab
        """
        self.name = name  # Element Symbol
        self.molar_mass = 0  # Total molar mass of all elements in slab layer
        self.density = 0  # Density of the molecule (g/cm^3)
        self.thickness = 0  # Thickness of the layer (Angstrom)
        self.roughness = 0  # Roughness of the surface (Angstrom)
        self.stoichiometry = stoichiometry  # Stoichiometry of the element
        self.poly_ratio = 1  # List that keeps track of polymorphous density ratio
        self.polymorph = []  # A list that contains the 'names' of various forms of the element (e.g. ions)
        self.gamma = 90  #
        self.phi = 90  #
        self.mag_density = []  # The scalling factor we want to multiply our scattering factor by (density is not the correct description)
        self.scattering_factor = name  # Identifies which scattering factor to be used. This parameter will allow us to implement 'scattering functions'
        self.mag_scattering_factor = None
        self.position = None


def get_number(string):
    """
    Purpose: Strip successive digits from a string
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

    # Used in error check
    symbols = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni',
               'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se' , 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
               'Cs', 'Ba', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po' , 'At', 'Rn', 'Fr', 'Ra', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds',
               'Rg', 'Cn', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk',
               'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr']

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
    file = open("Perovskite_Density.txt", "r")
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
    file = open("Atomic_Mass.txt", "r")
    lines = file.readlines()
    for line in lines:
        if line.split()[0] == atom:
            mass = line.split()[1]
    file.close()

    if mass == None:
        raise NameError("Inputted formula not found in perovskite density database")

    return float(mass)

def error_function(t, sigma1, sigma2, offset1, offset2):
    result1 = (erf((t-offset1)/sigma1/sqrt(2))+1)/2
    result2 = (erf((offset2-t)/sigma2/sqrt(2))+1)/2
    result = result1*result2
    return result

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


    def addlayer(self, num_layer, formula, thickness, density=None, roughness=0, link=None):
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

        # Thickness type check
        if type(thickness) != int and type(thickness) != float:
            raise TypeError('Layer ' +str(num_layer)+': Thickness must be integer or float type')

         # Checks if thickness is in a reasonable range
        if thickness < 0:
            warnings.warn('Layer ' +str(num_layer)+': Thickness must be a positive number. The absolute value was taken as it was assumed the user meant to input a negative value.')
            thickness = abs(thickness)
        elif thickness == 0:
            raise ValueError('Layer ' +str(num_layer)+': Thickness cannot be zero')

        # Checks Density
        if density == None:
            pass
        elif type(density) != int and type(density) != float:
            raise TypeError('Layer ' +str(num_layer)+': Density must be entered in as a float or integer type')
        elif density < 0:
            warnings.warn('Layer ' +str(num_layer)+': Density must be positive. The absolute value was taken as it was assumed the user meant to input a negative value.')
            density = abs(density)
        #elif density == 0:
        #    raise ValueError('Layer ' +str(num_layer)+': The density of a material can not be zero')

        if density == None:
            pass
        elif density > 20:
            warnings.warn('Layer ' +str(num_layer)+': The density of ' + str(density) + ' g/cm^3 might be too large of a value. Consider double checking your density. ')

        # Checks Roughness
        if type(roughness) != int and type(roughness) != float:
            raise TypeError('Layer ' +str(num_layer)+': Roughness must be of float or integer type.')
        elif roughness < 0:
            roughness = abs(roughness)
            warnings.warn('Layer ' +str(num_layer)+': Roughness should be entered as a positive value. The absolute value was taken as it was assumed the user meant to input a negative value.')

        if roughness > 15:
            warnings.warn('Layer ' +str(num_layer)+': Roughness is much larger than expected')

        if link == None:
            pass
        elif type(link) != list:
            raise TypeError('Variable link must be a list')
        elif len(link) != num_elements:
            raise RuntimeError('Layer ' +str(num_layer)+'Length of link must match number of elements in formula')
        elif len(set([type(x) for x in link])) != 1:
            raise TypeError('Layer ' +str(num_layer)+ ': All elements in link list must be booleans')

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

            elements[key].density = density  # sets density  (g/cm^3)
            elements[key].thickness = thickness  # sets thickness  (Angstrom)
            elements[key].roughness = roughness  # Order of Angstrom
            elements[key].molar_mass = molar_mass  # Molar mass of perovskite material
            elements[key].position = position
            position = position + 1

        self.structure[num_layer] = elements  # sets the layer with the appropriate slab properties

        # Determines if A/B site is linked to the A/B site on the next layer
        if self.number_layers == num_layer-1:  # Roughness link not required for last layer
            link = None
        elif link == None:
            mylink = [True for i in range(num_elements)]
            self.link.append(mylink)
        else:
            self.link.append(link)

    def polymorphous(self, lay, ele, polymorph, poly_ratio, sf=None):
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
        if type(polymorph) != list:
            raise TypeError('Layer ' + str(lay) + ': Input of polymorphous variable must be of type list')

        # Checking poly_ratio type
        if type(poly_ratio) != list:
            raise TypeError('Layer ' + str(lay) + ': Input of poly_ratio must be entered as a list')

        # Checks that length of list is greater than 2
        if len(polymorph) < 2:
            raise RuntimeError('Layer ' + str(lay) + ': Input more than one different version of the material. Length of list should be greater than two.')

        # Checks that poly_ratio and polymorphous have same length
        if len(poly_ratio) != len(polymorph):
            raise RuntimeError('Layer ' + str(lay) + ': variables poly_ratio and polymorph must be same length.')

        # Checks that all polymorphous types are strings
        if len(set([type(x) for x in polymorph])) != 1:
            raise TypeError('Layer ' +str(lay)+ ': variable polymorph must only contain strings.')

        # Chekcs that all poly_ratio types are floating points
        if len(set([type(x) for x in poly_ratio])) != 1:
            raise TypeError('Layer ' +str(lay)+ ': variable poly_ratio must only be of type float.')

        if abs(1 - sum(poly_ratio)) > 1e-3:
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
        if sf == None:
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
        if type(identifier) != str and type(identifier) != list:
            raise TypeError('Variable identifier must be a list or string.')

            # Checks to make sure that the variable identifier is either a list or string
        if type(sf) != str and type(sf) != list:
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

        # Checks that all identifiers are strings
        if type(identifier) == list:
            if type(identifier[0]) != str:
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
        layer = self.structure[lay]
        poly_start = True
        if type(identifier) == list:
            for key in list(layer.keys()):
                for idx in range(len(identifier)):
                    if len(layer[key].polymorph) != 0:  # polymorph case
                        if identifier[idx] in layer[key].polymorph:  # checks if identifier in selected polymorph
                            # pre-initializing for polymorph case
                            poly_idx = layer[key].polymorph.index(identifier[idx])  # determine index of poly

                            # initialization
                            if poly_start:  # first polymorph appearance
                                self.structure[lay][key].mag_scattering_factor = [0 for i in range(len(layer[key].polymorph))]
                                self.structure[lay][key].mag_density = np.zeros(len(layer[key].polymorph))
                                self.mag_elements[key] = [0 for i in range(len(layer[key].polymorph))]
                                poly_start = False
                                # gamma and phi entered as multiple arrays
                                if type(gamma) == list and type(phi) == list:
                                    self.structure[lay][key].gamma = np.zeros(len(layer[key].polymorph))
                                    self.structure[lay][key].phi = np.zeros(len(layer[key].polymorph))
                            self.structure[lay][key].mag_scattering_factor[poly_idx] = sf[idx]
                            self.structure[lay][key].mag_density[poly_idx] = density[idx]
                            self.mag_elements[key][poly_idx] = identifier[idx]
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
            val = erf((t-offset)/rough/sqrt(2))
        else:  # Downward error function
            val = erf((offset-t) / rough / sqrt(2))

        return val

    def density_profile(self, step = 0.1):
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
        density_mag = {k: dict() for k in list(self.mag_elements.keys())}  # hold magnetic elements

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

        # Initializing thickness array
        transition = [0]  # array that contains thickness that slab transition occurs
        thick = 0  # contains the total film thickness

        # Find all the slab transitions and computes total thickness
        for layer in range(1, n):
            val = transition[layer-1] + list(self.structure[layer].values())[0].thickness
            transition.append(val)
            thick = thick + list(self.structure[layer].values())[0].thickness

        #step = 0.05  # thickness step size
        thickness = np.arange(-50,thick+15+step, step) # Creates thickness array

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

                # Loops through all layers
                for layer in range(n):
                    offset = transition[layer]
                    # The current layer
                    if ele in list(self.structure[layer].keys()):

                        # saves scattering factor to be used in computation
                        if ele not in self.find_sf[0]:
                            self.find_sf[0][ele] = self.structure[layer][ele].scattering_factor

                        position = self.structure[layer][ele].position  # position of element
                        sigma = self.structure[layer][ele].roughness  # roughness parameterization
                        current_density = self.structure[layer][ele].stoichiometry* self.structure[layer][ele].density/ self.structure[layer][ele].molar_mass  # current density
                        if layer == n - 1:  # Last layer
                            next_density = 0  # density of element in next layer
                        elif ele in list(self.structure[layer+1].keys()):  # element in next layer
                            next_density = self.structure[layer+1][ele].stoichiometry* self.structure[layer+1][ele].density/ self.structure[layer+1][ele].molar_mass
                        else:  # element not in the next layer
                            next_density = 0

                        begin = 0
                        if layer == 0:
                            begin = 1



                        const = (next_density - current_density) / 2
                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        density_struct[ele] = density_struct[ele] + const*erf_func + begin*current_density
                    else:
                        current_density = 0
                        if layer == n - 1:  # Last layer
                            next_density = current_density
                            sigma = 0
                        elif ele in list(self.structure[layer+1].keys()):
                            position = self.structure[layer+1][ele].position  # position of element
                            next_density = self.structure[layer+1][ele].stoichiometry* self.structure[layer+1][ele].density/ self.structure[layer+1][ele].molar_mass # next layer density
                            previous_element = list(self.structure[layer].keys())[position]
                            sigma = self.structure[layer][previous_element].roughness
                        else:
                            next_density = 0
                            sigma = 0

                        const = (next_density-current_density)/2
                        erf_func = self.error_function(thickness, sigma, offset, True) + 1
                        density_struct[ele] = density_struct[ele] + const * erf_func
                        density_struct[ele][density_struct[ele]<0] = 0

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

                for layer in range(n): # loops through all layers

                    offset = transition[layer]
                    # Element found in current layer
                    if ele in list(self.structure[layer].keys()):

                        position = self.structure[layer][ele].position  # position of element
                        sigma = self.structure[layer][ele].roughness  # roughness parameterization
                        current_density = self.structure[layer][ele].stoichiometry * self.structure[layer][ele].density*self.structure[layer][ele].poly_ratio / self.structure[layer][ele].molar_mass  # current density

                        if layer == n - 1:  # On last layer
                            next_density = np.zeros(pn)  # density of element in next layer
                        elif ele in list(self.structure[layer + 1].keys()):  # element in next layer
                            next_density = self.structure[layer + 1][ele].stoichiometry * self.structure[layer + 1][ele].density* self.structure[layer+1][ele].poly_ratio/ self.structure[layer + 1][ele].molar_mass
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

                            # Density normalization
                            density_poly[ele][poly] = density_poly[ele][poly] + (const[po]*erf_func + begin*current_density[po])

                            po = po + 1

                    else:  # Element not found in current layer


                        current_density = np.zeros(pn)
                        if layer == n - 1:  # Last layer
                            next_density = current_density
                            sigma = 0
                        elif ele in list(self.structure[layer + 1].keys()):
                            position = self.structure[layer + 1][ele].position  # position of element
                            previous_element = list(self.structure[layer].keys())[position]
                            next_density = self.structure[layer + 1][ele].stoichiometry * self.structure[layer + 1][ele].density * self.structure[layer+1][ele].poly_ratio / self.structure[layer + 1][ele].molar_mass  # next layer density

                            sigma = self.structure[layer][previous_element].roughness
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
                for layer in range(n):  # loops through all layers
                    offset = transition[layer]
                    # Element found in current layer
                    if ele in list(self.structure[layer].keys()):
                        position = self.structure[layer][ele].position  # position of element
                        sigma = self.structure[layer][ele].roughness  # roughness parameterization
                        #current_density = self.structure[layer][ele].stoichiometry * np.array(self.structure[layer][ele].mag_density) * np.array(self.structure[layer][ele].density) / self.structure[layer][ele].molar_mass  # current density
                        current_density = self.structure[layer][ele].stoichiometry * np.array(self.structure[layer][ele].mag_density)
                        if layer == n - 1:  # Last layer
                            next_density = np.zeros(pm)  # density of element in next layer
                        elif ele in list(self.structure[layer + 1].keys()):  # element in next layer
                            #next_density = self.structure[layer + 1][ele].stoichiometry * np.array(self.structure[layer + 1][ele].mag_density) * np.array(self.structure[layer + 1][ele].density) / self.structure[layer + 1][ele].molar_mass
                            next_density = self.structure[layer + 1][ele].stoichiometry * np.array(self.structure[layer + 1][ele].mag_density)
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

                            # finds magnetic scattering factors
                            if mag not in self.find_sf[1]:
                                self.find_sf[1][mag] = self.structure[layer][ele].mag_scattering_factor[ma]
                            # Density normalization
                            density_mag[ele][mag] = density_mag[ele][mag] + (const[ma] * erf_func + begin * current_density[ma])
                            ma = ma + 1

                    else:  # Element not found in current layer

                        current_density = np.zeros(pm)
                        if layer == n - 1:  # Last layer
                            next_density = current_density
                            sigma = 0
                        elif ele in list(self.structure[layer + 1].keys()):
                            position = self.structure[layer + 1][ele].position  # position of element
                            previous_element = list(self.structure[layer].keys())[position]
                            #next_density = self.structure[layer + 1][ele].stoichiometry * np.array(self.structure[layer + 1][ele].mag_density)*np.array(self.structure[layer + 1][ele].density) / self.structure[layer + 1][ele].molar_mass  # next layer density
                            next_density = self.structure[layer + 1][ele].stoichiometry * np.array(self.structure[layer + 1][ele].mag_density)
                            sigma = self.structure[layer][previous_element].roughness
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
                density_magnetic[mag] = density_mag[ele][mag]

        self.transition = transition  # keeps track of transition
        return thickness, density, density_magnetic

    def plot_density_profile(self, fig=1):
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

        plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02, 0.5))
        plt.xlabel('Thickness (Angstrom)')
        plt.ylabel('Density (mol/cm^3)')



    def reflectivity(self, E, qz, precision=1e-6,s_min = 0.1):

        """
        Purpose: Computes the reflectivity
        :param E: Energy of reflectivity scan (eV)
        :param qi: Starting momentum transfer
        :param qf: End momentum transfer
        :param precision: Precision value
        :param s_min: Minimum step size (None indicates using default value)
        :return:
            qz - momentum transfer
            R1 - reflectivity

        """
        # initialize the reflectivity array

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s
        wavelength = h * c / (E * 1e-10)  # wavelength m
        #wavelength = 19.366478131833802
        # requires angle for reflectivity computation and minimum slab thickness
        #theta_i = arcsin(qi / E / (0.001013546247)) * 180 / pi  # initial angle
        #theta_f = arcsin(qf / E / (0.001013546247)) * 180 / pi  # final angle in interval

        thickness, density, density_magnetic = self.density_profile(step=0.1)  # Computes the density profile

        sf = dict()  # form factors of non-magnetic components
        sfm = dict()  # form factors of magnetic components

        start_new = time_ns()
        # Non-Magnetic Scattering Factor
        for e in self.find_sf[0].keys():
            sf[e] = find_form_factor(self.find_sf[0][e], E, False)
        # Magnetic Scattering Factor
        for em in self.find_sf[1].keys():
            sfm[em] = find_form_factor(self.find_sf[1][em],E,True)

        delta, beta = index_of_refraction(density, sf, E)  # calculates dielectric constant for structural component
        delta_m, beta_m = magnetic_optical_constant(density_magnetic, sfm, E)   # calculates dielectric constant for magnetic component


        # definition as described in Lott Dieter Thesis
        n = 1 + np.vectorize(complex)(-delta, beta)
        #epsilon = 1 + np.vectorize(complex)(-2*delta, 2*beta)
        epsilon = n**2
        #Q = np.vectorize(complex)(delta, beta)
        Q = np.vectorize(complex)(beta_m, delta_m)
        epsilon_mag = Q*epsilon*2*(-1)
        end_new = time_ns()
        my_slabs = ALS(epsilon.real, epsilon_mag.real, precision)  # performs the adaptive layer segmentation using Numba

        my_slabs = my_slabs.astype(int)  # sets all numbers to integers

        #plt.figure()
        #plt.plot(thickness[my_slabs], epsilon.real[my_slabs], '.')
        #plt.legend([len(my_slabs)-1])
        #plt.show()

        my_slabs = my_slabs[1:]  # removes first element

        print("Layer Segmentation: ", end_new-start_new, " ns")
        start = time()
        m = len(my_slabs)  # number of slabs
        #m = len(epsilon)
        A =pr.Generate_structure(m)  # creates object for reflectivity computation
        m_j=0  # previous slab
        idx = 0  # keeps track of current layer
        layer = 0
        gamma = 90  # pre-initialize magnetization direction
        phi = 90

        for m_i in my_slabs:
            d = thickness[m_i] - thickness[m_j]  # computes thickness of slab
            eps = (epsilon[m_i] + epsilon[m_j])/2  # computes the dielectric constant value to use
            eps_mag = (epsilon_mag[m_i] + epsilon_mag[m_j])/2  # computes the magnetic dielectric constant

            # Determines the magnetization direction of the first layer
            if layer == 0:
                if self.layer_magnetized[0]:
                    for ele in self.structure[layer].keys():
                        if self.structure[0][ele].mag_scattering_factor != None:
                            gamma = self.structure[0][ele].gamma
                            phi = self.structure[0][ele].phi


            # Determines the magnetization direction of the other layers
            if self.transition[layer]<=thickness[m_j] and layer<len(self.transition)-1:
                layer = layer + 1
                if self.layer_magnetized[layer]:
                    for ele in self.structure[layer].keys():
                        if self.structure[layer][ele].mag_scattering_factor != None:
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

                A[idx].seteps([eps,eps,eps,eps_mag])  # sets dielectric tensor for magnetic layer
            else:
                A[idx].seteps(eps)  # sets dielectric tensor for non-magnetic layer


            if idx != 0:
                A[idx].setd(d)  # sets thickness of layer if and only if not substrate layer


            # move onto the next layer
            m_j = m_i
            idx = idx + 1


        Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / pi  # initial angle

        end = time()
        print('Initialization: ', end-start)
        start = time()
        Rtemp = pr.Reflectivity(A, Theta, wavelength, MagneticCutoff=1e-10)  # Computes the reflectivity
        end = time()
        print('Reflectivity: ', end-start)
        R = dict()
        print()
        """
        # Used to demonstrate how sample is being segmented
        plot_t = np.insert(thickness[my_slabs],0,thickness[0], axis=0)
        plot_e = np.insert(real(epsilon[my_slabs]),0, real(epsilon[0]), axis=0)
        """
        if len(Rtemp) == 2:
            R['S'] = Rtemp[0]  # s-polarized light
            R['P'] = Rtemp[1]  # p-polarized light
            R['AL'] = (Rtemp[0]-Rtemp[1])/(Rtemp[0]+Rtemp[1])  # Asymmetry linear polarized
            R['LC'] = np.zeros(len(Rtemp[0]))  # Left circular
            R['RC'] = np.zeros(len(Rtemp[0]))  # right circular
            R['AC'] = np.zeros(len(Rtemp[0]))  # Asymmetry circular polarized
        elif len(Rtemp)==4:
            R['S'] = Rtemp[0]
            R['P'] = Rtemp[1]
            R['AL'] = (Rtemp[0]-Rtemp[1])/(Rtemp[0]+Rtemp[1])
            R['LC'] = Rtemp[2]
            R['RC'] = Rtemp[3]
            R['AC'] = (Rtemp[2]-Rtemp[3])/(Rtemp[2]+Rtemp[3])
        else:
            raise TypeError('Error in reflectivity computation. Reflection array not expected size.')

        return qz, R

    def energy_scan(self, Theta, energy, precision=1e-6,s_min = 0.1):
        """
                Purpose: Computes the reflectivity
                :param E: Energy of reflectivity scan (eV)
                :param qi: Starting momentum transfer
                :param qf: End momentum transfer
                :param precision: Precision value
                :param s_min: Minimum step size (None indicates using default value)
                :return:
                    qz - momentum transfer
                    R1 - reflectivity

                """
        Elen = len(energy)
        R = {'S': np.zeros(Elen),
             'P': np.zeros(Elen),
             'AL': np.zeros(Elen),
             'LC': np.zeros(Elen),
             'RC': np.zeros(Elen),
             'AC': np.zeros(Elen)}

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s

        # wavelength = 19.366478131833802
        # requires angle for reflectivity computation and minimum slab thickness
        # theta_i = arcsin(qi / E / (0.001013546247)) * 180 / pi  # initial angle
        # theta_f = arcsin(qf / E / (0.001013546247)) * 180 / pi  # final angle in interval

        thickness, density, density_magnetic = self.density_profile(step=s_min)  # Computes the density profile



        for E in range(len(energy)):

            wavelength = h * c / (energy[E] * 1e-10)  # wavelength m
            sf = dict()  # form factors of non-magnetic components
            sfm = dict()  # form factors of magnetic components

            # Non-Magnetic Scattering Factor
            for e in self.find_sf[0].keys():
                sf[e] = find_form_factor(self.find_sf[0][e], energy[E], False)

            # Magnetic Scattering Factor
            for em in self.find_sf[1].keys():
                sfm[em] = find_form_factor(self.find_sf[1][em], energy[E], True)


            delta, beta = index_of_refraction(density, sf, energy[E])  # calculates dielectric constant for structural component
            delta_m, beta_m = magnetic_optical_constant(density_magnetic, sfm, energy[E])   # calculates dielectric constant for magnetic component



            # definition as described in Lott Dieter Thesis
            n = 1 + np.vectorize(complex)(-delta, beta)
            # epsilon = 1 + np.vectorize(complex)(-2*delta, 2*beta)
            epsilon = n ** 2
            # Q = np.vectorize(complex)(delta, beta)
            Q = np.vectorize(complex)(beta_m, delta_m)
            epsilon_mag = Q * epsilon * 2 * (-1)

            my_slabs = ALS(epsilon.real,epsilon_mag.real, precision=1e-6)
            my_slabs = my_slabs[1:]
            my_slabs = my_slabs.astype(int)




            m = len(my_slabs)  # number of slabs

            A = pr.Generate_structure(m)  # creates object for reflectivity computation
            m_j = 0  # previous slab
            idx = 0  # keeps track of current layer
            layer = 0
            gamma = 90  # pre-initialize magnetization direction
            phi = 90

            for m_i in my_slabs:
                d = thickness[m_i] - thickness[m_j]  # computes thickness of slab
                eps = (epsilon[m_i] + epsilon[m_j]) / 2  # computes the dielectric constant value to use
                eps_mag = (epsilon_mag[m_i] + epsilon_mag[m_j]) / 2  # computes the magnetic dielectric constant

                # Determines the magnetization direction of the first layer
                if layer == 0:
                    if self.layer_magnetized[0]:
                        for ele in self.structure[layer].keys():
                            if self.structure[0][ele].mag_scattering_factor != None:
                                gamma = self.structure[0][ele].gamma
                                phi = self.structure[0][ele].phi

                # Determines the magnetization direction of the other layers
                if self.transition[layer] <= thickness[m_j] and layer < len(self.transition) - 1:
                    layer = layer + 1
                    if self.layer_magnetized[layer]:
                        for ele in self.structure[layer].keys():
                            if self.structure[layer][ele].mag_scattering_factor != None:
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

                    A[idx].seteps([eps, eps, eps, eps_mag])  # sets dielectric tensor for magnetic layer
                else:
                    A[idx].seteps(eps)  # sets dielectric tensor for non-magnetic layer

                if idx != 0:
                    A[idx].setd(d)  # sets thickness of layer if and only if not substrate layer

                # move onto the next layer
                m_j = m_i
                idx = idx + 1


            #Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / pi  # initial angle
            Rtemp = pr.Reflectivity(A, Theta, wavelength, MagneticCutoff=1e-10)  # Computes the reflectivity

            if len(Rtemp) == 2:
                R['S'][E] = Rtemp[0][0]  # s-polarized light
                R['P'][E] = Rtemp[1][0]  # p-polarized light
                R['AL'][E] = (Rtemp[0][0] - Rtemp[1][0]) / (Rtemp[0][0] + Rtemp[1][0])  # Asymmetry linear polarized
            elif len(Rtemp) == 4:
                R['S'][E] = Rtemp[0][0]
                R['P'][E] = Rtemp[1][0]
                R['AL'][E] = (Rtemp[0][0] - Rtemp[1][0]) / (Rtemp[0][0] + Rtemp[1][0])
                R['LC'][E] = Rtemp[2][0]
                R['RC'][E] = Rtemp[3][0]
                R['AC'][E] = (Rtemp[2][0] - Rtemp[3][0]) / (Rtemp[2][0] + Rtemp[3][0])
            else:
                raise TypeError('Error in reflectivity computation. Reflection array not expected size.')


        return energy, R

    def energy_scan_multi(self, Theta, energy, precision=0.5,s_min = 0.1):
        """
                        Purpose: Computes the reflectivity
                        :param E: Energy of reflectivity scan (eV)
                        :param qi: Starting momentum transfer
                        :param qf: End momentum transfer
                        :param precision: Precision value
                        :param s_min: Minimum step size (None indicates using default value)
                        :return:
                            qz - momentum transfer
                            R1 - reflectivity

                        """

        h = 4.135667696e-15  # Plank's constant eV*s
        c = 2.99792458e8  # speed of light m/s

        # wavelength = 19.366478131833802
        # requires angle for reflectivity computation and minimum slab thickness
        # theta_i = arcsin(qi / E / (0.001013546247)) * 180 / pi  # initial angle
        # theta_f = arcsin(qf / E / (0.001013546247)) * 180 / pi  # final angle in interval

        thickness, density, density_magnetic = self.density_profile(step=s_min)  # Computes the density profile

        fsf = self.find_sf
        st = self.structure
        lm = self.layer_magnetized
        tran = self.transition



        prod = partial(multi_energy_calc,thickness, density, density_magnetic, fsf, st, lm, tran, Theta, precision)
        cores = multiprocessing.cpu_count()
        with multiprocessing.Pool(cores) as pool:
            result_list = pool.map(prod, energy)
        pool.join()



        return energy, result_list

def multi_energy_calc(thickness, density, density_magnetic, find_sf, structure, layer_magnetized, transition, theta, prec , E):

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    wavelength = h * c / (E * 1e-10)  # wavelength m
    sf = dict()  # form factors of non-magnetic components
    sfm = dict()  # form factors of magnetic components

    # Non-Magnetic Scattering Factor
    for e in find_sf[0].keys():
        sf[e] = find_form_factor(find_sf[0][e], E, False)
    # Magnetic Scattering Factor
    for em in find_sf[1].keys():
        sfm[em] = find_form_factor(find_sf[1][em], E, True)

    delta, beta = index_of_refraction(density, sf, E)  # calculates dielectric constant for structural component
    delta_m, beta_m = magnetic_optical_constant(density_magnetic, sfm, E)  # calculates dielectric constant for magnetic component

    # definition as described in Lott Dieter Thesis
    n = 1 + np.vectorize(complex)(-delta, beta)
    # epsilon = 1 + np.vectorize(complex)(-2*delta, 2*beta)
    epsilon = n ** 2
    # Q = np.vectorize(complex)(delta, beta)
    Q = np.vectorize(complex)(beta_m, delta_m)
    epsilon_mag = Q * epsilon * 2 * (-1)

    my_slabs = layer_segmentation(thickness, epsilon, epsilon_mag, prec)  # computes the layer segmentation

    m = len(my_slabs)  # number of slabs
    A = pr.Generate_structure(m)  # creates object for reflectivity computation
    m_j = 0  # previous slab
    idx = 0  # keeps track of current layer
    layer = 0
    gamma = 90  # pre-initialize magnetization direction
    phi = 90

    for m_i in my_slabs:
        d = thickness[m_i] - thickness[m_j]  # computes thickness of slab
        eps = (epsilon[m_i] + epsilon[m_j]) / 2  # computes the dielectric constant value to use
        eps_mag = (epsilon_mag[m_i] + epsilon_mag[m_j]) / 2  # computes the magnetic dielectric constant

        # Determines the magnetization direction of the first layer
        if layer == 0:
            if layer_magnetized[0]:
                for ele in structure[layer].keys():
                    if structure[0][ele].mag_scattering_factor != None:
                        gamma = structure[0][ele].gamma
                        phi = structure[0][ele].phi

        # Determines the magnetization direction of the other layers
        if transition[layer] <= thickness[m_j] and layer < len(transition) - 1:
            layer = layer + 1
            if layer_magnetized[layer]:
                for ele in structure[layer].keys():
                    if structure[layer][ele].mag_scattering_factor != None:
                        gamma = structure[layer][ele].gamma
                        phi = structure[layer][ele].phi

        # sets the magnetization direction based on the input angles
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

    # Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / pi  # initial angle

    R = pr.Reflectivity(A, theta, wavelength, MagneticCutoff=1e-10)  # Computes the reflectivity


    return R


if __name__ == "__main__":
    """
    # Example 1: Simple sample creation
    sample = slab(5)  # Initializing four layers
    s = 0.1

    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, roughness=2, link=[False,True,True])  # substrate layer
    sample.polymorphous(0,'Ti',['Ti3+','Ti4+'], [0.25,0.75])

    # First Layer
    # Link: La-->La and O-->O
    # * denotes optional parameter
    sample.addlayer(1, 'LaMnO3', 25, roughness=1.5, link=[True,False,True])  # Film 1 on top of substrate
    sample.polymorphous(1,'Mn',['Mn2+','Mn3+'], [0.9,0.1], sf=['Mn','Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [6.1,7.1], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)



    sample.addlayer(2, 'LaMnO3',5, roughness=2, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [0.1 , 0.9], sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [1, 2], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    sample.addlayer(3, 'LaMnO3', 5, roughness= 2, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.1, 0.9], sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.5, 5], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    sample.addlayer(4, 'COC2', 5, roughness= 1, density=0.5, link=[True, False, True])  # Film 1 on top of substrate
    """

    """
    # Example 2: Simple sample creation
    sample = slab(6)  # Initializing four layers
    s = 0.1

    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, roughness=2, link=[False, True, True])  # substrate layer
    sample.polymorphous(0, 'Ti', ['Ti3+', 'Ti4+'], [0.25, 0.75])

    # First Layer
    # Link: La-->La and O-->O
    # * denotes optional parameter
    sample.addlayer(1, 'LaMnO3', 25, roughness=1.5, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(1, 'Mn', ['Mn2+', 'Mn3+'], [0.9, 0.1],
                        sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [6.1, 7.1],
                         ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    # Second Layer
    # Link: La-->C and O-->C
    sample.addlayer(2, 'LaAlO3', 16, density=5, roughness=2, link=[True, False, True])  # Film 2 on top film 1
    sample.magnetization(2, 'Al', 5, 'Co')  # mag_type is preset to 'isotropic

    sample.addlayer(3, 'LaMnO3', 5, roughness=2, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.1, 0.9],
                        sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [1, 2],
                         ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    sample.addlayer(4, 'LaMnO3', 5, roughness=2, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [0.1, 0.9],
                        sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.5, 5],
                         ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    sample.addlayer(5, 'COC2', 5, roughness=1, density=0.5, link=[True, False, True])  # Film 1 on top of substrate

    # Impurity on surface
    # Impurity takes the form 'CCC'
    # The exact implementation of this has not quite been determined yet
    #sample.addlayer(3, 'CCC', 10, density=1) #  Density initialized to 5g/cm^3
    hi
    """


    # Example 2: Simple sample creation
    sample = slab(6)  # Initializing four layers
    s = 0.1
    mag_dense = 0.1
    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, density = 5.120891853,roughness=2, link=[False, True, True])  # substrate layer
    sample.addlayer(1, 'SrTiO3', 4, density=5.120891853, roughness=2, link=[False, True, True])  # substrate layer


    sample.addlayer(2,'LaMnO3', 4, density = 6.8195658, roughness=2)
    sample.polymorphous(2,'Mn', ['Mn2+','Mn3+'], [1,0], sf=['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+','Mn3+'], [0,0],['Co','Ni'])

    sample.addlayer(3, 'LaMnO3', 30, density = 6.8195658, roughness=2)
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [mag_dense, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 4, density = 6.8195658, roughness=2)
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 4, density=6.8195658, roughness=2)
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0, 0], ['Co', 'Ni'])

    #sample.addlayer(4, 'CCC', 4, density = 0, roughness = 2)

    """
    result = sample.structure[2]

    e = 'O'
    print('STRUCTURE')
    print('Name: ', result[e].name)
    print('Scattering Factor: ', result[e].scattering_factor)
    print('Stoichiometry: ',result[e].stoichiometry)
    print('Polymorph: ',result[e].polymorph)
    print('Polymorph Ratio : ', result[e].poly_ratio)
    print('Link: ', sample.link)
    print('')
    print('MAGNETIZATION')
    print('Scattering Factors: ', result[e].mag_scattering_factor)
    print('Density: ', result[e].mag_density)
    print('Phi: ', result[e].phi, ' degrees')
    print('Theta: ', result[e].theta, ' degrees')
    print('Position: ', result[e].position)

    #  print(sample.myelements)
    #  print(sample.poly_elements)
    #print(sample.mag_elements)
    
    
    # Example 2: Simple sample creation
    sample = slab(2)  # Initializing four layers
    s = 0.1
    mag_dense = 0.1
    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, density=0, roughness=2, link=[False, True, True])  # substrate layer


    sample.addlayer(1, 'LaMnO3', 30, density=6.8195658, roughness=2)
    sample.polymorphous(1, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [mag_dense, 0], ['Co', 'Ni'])

    """
    sample.plot_density_profile()
    thickness, density, density_magnetic = sample.density_profile(step=0.1)



    E = 642.2 # eV

    elements = density.keys()
    sf = dict()
    elements_mag = density_magnetic.keys()
    sfm = dict()
    # Non-Magnetic Scattering Factor
    for e in sample.find_sf[0].keys():
        sf[e] = find_form_factor(sample.find_sf[0][e], E, False)
    # Magnetic Scattering Factor
    for em in sample.find_sf[1].keys():
        sfm[em] = find_form_factor(sample.find_sf[1][em], E, True)

    delta, beta = index_of_refraction(density, sf, E)


    plt.figure(2)
    plt.plot(thickness, delta, thickness, beta)
    plt.suptitle('Optical Profile')
    plt.xlabel('Thickness')
    plt.ylabel('Profile')
    plt.legend(['delta','beta'])

    h = 4.1257e-15 # Plank's constant eV*s
    c = 2.99792458e8 # speed of light m/s
    y = h*c/(E*1e-10)  # wavelength m

    F = np.loadtxt('test_example.txt')
    qz = F[:,0]
    p1 = 1e-5
    p2 = 1e-5
    start = time()
    qz, Rtemp =  sample.reflectivity(E, qz,precision=0,s_min=0.1)  # baseline
    end = time()
    print(end-start)
    start = time()
    qz1, R1temp = sample.reflectivity(E, qz, precision = p1, s_min=0.1)
    end = time()
    print(end-start)
    qz2, R2temp = sample.reflectivity(E,qz, precision = p2, s_min=0.1)

    #R = np.log10(Rtemp['LC'])
    #R1 = np.log10(R1temp['LC'])
    #R2 = np.log10(R2temp['LC'])

    R = Rtemp['AC']
    R1 = R1temp['AC']
    R2 = R2temp['AC']


    plt.figure(4)
    plt.plot(qz, R, 'k-')
    plt.plot(qz1, R1, 'b--')
    plt.plot(qz2, R2, 'r--')
    plt.legend(['baseline',str(p1), str(p2)])
    #plt.yscale("log")
    plt.xlabel('qz')
    plt.ylabel('Reflectivity ' + "$(log_{10}(R))$")
    plt.title('ReMagX vs. Python Script (800 eV)')

    diff_1 = abs(R-R1)/abs(R+R1)
    figure(5)
    plt.plot(qz1, diff_1)
    plt.suptitle("Difference in Spectra: " + str(p1))
    plt.xlabel("Thickness (Angstroms)")
    plt.ylabel("Percent Difference")

    diff_2 = abs(R-R2)/abs(R+R2)
    figure(6)
    plt.plot(qz1, diff_2)
    plt.suptitle("Difference in Spectra: " + str(p2))
    plt.xlabel("Thickness (Angstrom)")
    plt.ylabel("Percent Difference")


    max1 = max(abs(R-R1))
    max2 = max(abs(R-R2))
    A1 = simpson(abs(R-R1), qz)
    A2 = simpson(abs(R-R2), qz)

    print()
    print()
    print(tabulate([[p1, max1, A1], [p2, max2, A2]], headers=['Precision', 'Maximum', 'Total Area']))



    q = F[:,0]
    I = F[:,1]
    # I = np.log10(I)
    plt.figure(55)
    plt.plot(q,I,'k')
    plt.plot(qz, R, 'r--')
    plt.suptitle('Zak Formalism: Asymmetry 642.2 eV (rho=0.5) ')
    plt.xlabel('qz')
    plt.ylabel('Reflectivity ' + "$(log_{10}(R))$")
    plt.legend(['ReMagX','Lucas'])

    diff_3 = abs(I-R)/abs(I+R)
    plt.figure(56)
    plt.plot(q, diff_3, 'k')
    plt.suptitle('Percent Difference: precision = ' + str(0))
    plt.xlabel('qz')
    plt.ylabel('Percent Difference')
    plt.legend(['ReMagX', 'Lucas'])

    plt.show()










