import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.special import erf
import warnings
from KK_And_Merge import *

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
        self.theta = 0  #
        self.phi = 0  #
        self.mag_density = None  # The scalling factor we want to multiply our scattering factor by (density is not the correct description)
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
    mydict = dict()
    while len(formula) > 0:
        formula, info = checkstring(formula)
        mydict[info[0]] = element(info[0], info[1])
        num_elements = num_elements + 1

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

        self.structure = [dict() for i in range(num_layers)]  # keeps track of the structural properties of each layer
        self.poly_elements = dict()  # Keeps track of the  polymorphous elements
        self.mag_elements = dict()  # Keeps track of the magnetized elements in the material

    def addlayer(self, num_layer, formula, thickness, density=None, roughness=0, link=None):
        """
        Purpose: Add layer to sample material
        :param num_layer: Layer number (0 indicates substrate layer)
        :param formula: Chemical formula
        :param thickness: Thickness of the material in Angstrom 10^(-10)m
        :param density: Density of the perovskite material (g/cm^3). No input indicates user wants to use density found in database
        :param roughness: Rougness at the interface
        :param link: Determines if roughness between two of the same sites are linked
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
        elif density == 0:
            raise ValueError('Layer ' +str(num_layer)+': The density of a material can not be zero')

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
            molar_mass = molar_mass + atomic_mass(key)*elements[key].stoichiometry

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
        n = len(elements)
        if link == None:
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
        self.structure[lay][ele].poly_ratio = poly_ratio

        # Sets scattering factors
        if sf == None:
            self.structure[lay][ele].scattering_factor = polymorph
        else:
            self.structure[lay][ele].scattering_factor = sf

    def magnetization(self, lay, identifier, density, sf, phi=0, theta=0):
        """
        Purpose: Set magnetization properties
        :param lay: Layer the magnetic material is found
        :param identifier: The symbol of the magnetic element
        :param density: The density of the magnetism
        :param sf: The scattering factor you want to use for the magnetic component
        :param phi: Azimuthal angle (degrees)
        :param theta: Polar angle (degrees)
        """

        # ---------------------------------------- Error Check ------------------------------------------------------- #

        # Checking that layer exists
        if lay > len(self.structure) - 1:
            raise RuntimeError('Layer ' + str(lay) + ' selected, but maximum layer is ' + str(len(self.structure) - 1))

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


        # ------------------------------------------------------------------------------------------------------------ #
        layer = self.structure[lay]
        if type(identifier) == list:

            # Loop through all elements in selected layer
            for key in list(layer.keys()):
                if set(layer[key].polymorph) == set(identifier): # Checks if elements in identifier and polymorph list the same
                    if layer[key].polymorph == identifier:  # Checks if they are in the same order, if not switch order
                        self.structure[lay][key].mag_scattering_factor = sf
                        self.structure[lay][key].mag_density = density
                        self.structure[lay][key].phi = phi
                        self.structure[lay][key].theta = theta
                        self.mag_elements[key] = identifier
                    else:
                        temp_sf = ['X' for i in range(len(layer[key].polymorph))]
                        temp_density = [0 for i in range(len(layer[key].polymorph))]
                        for n in range(len(layer[key].polymorph)):
                            idx = identifier.index(layer[key].polymorph[n])
                            temp_sf[idx] = sf[n]
                            temp_density[idx] = density[n]
                        self.structure[lay][key].mag_scattering_factor = temp_sf
                        self.structure[lay][key].mag_density = temp_density
                        self.structure[lay][key].phi = phi
                        self.structure[lay][key].theta = theta
                        self.mag_elements[key] = identifier
                else:
                    # Plolymorph names do not match magnetic names
                    raise RuntimeError('Polymorph names '+str(layer[key].polymorph)+' do not match magnetic names ' + str(identifier))
        else:
            if identifier not in list(layer.keys()):
                raise NameError('Variable identifier ' +str(identifier)+' does not match any elements in selected layer.')
            # Case where magnetization does not belong to polymorphs
            self.structure[lay][identifier].mag_scattering_factor = [sf]
            self.structure[lay][identifier].mag_density = [density]
            self.structure[lay][identifier].phi = phi
            self.structure[lay][identifier].phi = theta
            self.mag_elements[identifier] = [identifier]


    def convert_to_slice(self):
            print('Time to convert')

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

    def density_profile(self):
        """
        Purpose: Creates the density profile based on the slab properties
        :return: thickness - thickness array in angstrom
                 density - structural density array in mol/cm^3
                 mag_density - magnetic density array in mol/cm^3
        """

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

        step = 1e-4  # thickness step size
        thickness = np.arange(-25,thick+15+step, step) # Creates thickness array

        # Loop through elements in sample
        for ele in self.myelements:

            # structural elements (none polymorphous or magnetic)
            if ele in struct_keys:
                density_struct[ele] = np.zeros(len(thickness))  # sets array same length as thickness array full of zeros

                # Loops through all layers
                for layer in range(n):
                    if ele in list(self.structure[layer].keys()): # determines if desired element is in current layer
                        if layer == 0:  # first layer
                            offset = transition[0]  # offset value for error function
                            rough = self.structure[layer][ele].roughness  # surface roughness
                            val = (self.error_function(thickness,rough, offset, False) + 1)/2 # calculate error function
                            ratio = self.structure[layer][ele].stoichiometry / self.structure[layer][ele].molar_mass # ratio calculation
                            density_struct[ele] = density_struct[ele] + val*self.structure[layer][ele].density*ratio  # Compute density

                        else:  # not first layer
                            position = self.structure[layer][ele].position  # position of element
                            offset1 = transition[layer-1]  # thickness offset of upward error function
                            offset2 = transition[layer]  # thickness offset for downward error function

                            # gets element of previous layer that is in the same site as current element (e.g. A-site for ABO3)
                            previous_element = list(self.structure[layer-1].keys())[position]
                            rough1 = self.structure[layer-1][previous_element].roughness  # roughness of previous element
                            rough2 = self.structure[layer][ele].roughness # roughness of surface
                            val1 = self.error_function(thickness, rough1, offset1, True)  # computes error function
                            val2 = self.error_function(thickness, rough2, offset2, False)  # computes error function

                            val = val1 + val2  # adds both error functions together
                            ratio = self.structure[layer][ele].stoichiometry / self.structure[layer][ele].molar_mass  # computes ratio
                            density_struct[ele] = density_struct[ele] + val * self.structure[layer][ele].density * ratio # computes density

            # Polymorphous elements
            if ele in poly_keys:
                first = True
                for layer in range(n):
                    if ele in list(self.structure[layer].keys()):
                        if first:
                            density_poly[ele] = {k:np.zeros(len(thickness)) for k in list(self.structure[layer][ele].polymorph)}
                            first = False
                        if layer == 0:
                            offset = transition[0]  # offset value for error function
                            rough = self.structure[layer][ele].roughness  # surface roughness
                            val = (self.error_function(thickness, rough, offset,False) + 1) / 2  # calculate error function
                            ratio = self.structure[layer][ele].stoichiometry / self.structure[layer][ele].molar_mass

                            po = 0
                            for poly in list(density_poly[ele].keys()):
                                density_poly[ele][poly] = density_poly[ele][poly] + val * self.structure[layer][ele].density * ratio * self.structure[layer][ele].poly_ratio[po]
                                po = po + 1
                        else:
                            position = self.structure[layer][ele].position  # position of element
                            offset1 = transition[layer - 1]
                            offset2 = transition[layer]

                            previous_element = list(self.structure[layer - 1].keys())[position]
                            rough1 = self.structure[layer - 1][previous_element].roughness
                            rough2 = self.structure[layer][ele].roughness
                            val1 = self.error_function(thickness, rough1, offset1, True)
                            val2 = self.error_function(thickness, rough2, offset2, False)
                            val = val1 + val2
                            ratio = self.structure[layer][ele].stoichiometry / self.structure[layer][ele].molar_mass
                            po = 0
                            for poly in list(density_poly[ele].keys()):
                                density_poly[ele][poly] = density_poly[ele][poly] + val * self.structure[layer][ele].density * ratio * self.structure[layer][ele].poly_ratio[po]
                                po = po + 1


            if ele in mag_keys:
                first = True
                for layer in range(n):
                    if ele in list(self.structure[layer].keys()):
                        if first:
                            density_mag[ele] = {k:np.zeros(len(thickness)) for k in list(self.mag_elements[ele])}
                            first = False

                        if layer == 0:
                            offset = transition[0]  # offset value for error function
                            rough = self.structure[layer][ele].roughness  # surface roughness
                            val = (self.error_function(thickness, rough, offset,False) + 1) / 2  # calculate error function
                            ratio = self.structure[layer][ele].stoichiometry / self.structure[layer][ele].molar_mass

                            ma = 0
                            for mag in self.mag_elements[ele]:
                                density_mag[ele][mag] = density_mag[ele][mag] + val * ratio * self.structure[layer][ele].mag_density[ma]
                                ma = ma + 1
                        else:
                            position = self.structure[layer][ele].position  # position of element
                            offset1 = transition[layer - 1]
                            offset2 = transition[layer]

                            previous_element = list(self.structure[layer - 1].keys())[position]
                            rough1 = self.structure[layer - 1][previous_element].roughness
                            rough2 = self.structure[layer][ele].roughness
                            val1 = self.error_function(thickness, rough1, offset1, True)
                            val2 = self.error_function(thickness, rough2, offset2, False)
                            val = val1 + val2
                            ratio = self.structure[layer][ele].stoichiometry / self.structure[layer][ele].molar_mass
                            ma = 0

                            for mag in self.mag_elements[ele]:
                                density_mag[ele][mag] = density_mag[ele][mag] + val  * ratio * self.structure[layer][ele].mag_density[ma]
                                ma = ma + 1

        # Create single dictionary to use
        density = density_struct
        for ele in list(density_poly.keys()):
            for poly in list(density_poly[ele].keys()):
                density[poly] = density_poly[ele][poly]

        density_magnetic = dict()
        for ele in list(density_mag.keys()):
            for mag in list(density_mag[ele].keys()):
                density_magnetic[mag] = density_mag[ele][mag]
        return thickness, density, density_magnetic






if __name__ == "__main__":

    # Example: Simple sample creation
    sample = slab(6)  # Initializing four layers

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

    # Second Layer
    # Link: La-->C and O-->C
    sample.addlayer(2, 'LaAlO3', 16, density=5, roughness=2, link=[True, False,True])   # Film 2 on top film 1
    sample.magnetization(2,'Al', 5, 'Co') # mag_type is preset to 'isotropic

    sample.addlayer(3, 'LaMnO3', 15, roughness=0.3, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.1 , 0.9], sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [1, 2], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    sample.addlayer(4, 'LaMnO3', 2, roughness= 0.3, link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [0.1, 0.9], sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.5, 8], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    sample.addlayer(5, 'LaMnO3', 2, roughness=2 , link=[True, False, True])  # Film 1 on top of substrate
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [0.01, 0.99], sf=['Mn', 'Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0, 0], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)
    # Impurity on surface
    # Impurity takes the form 'CCC'
    # The exact implementation of this has not quite been determined yet
    #sample.addlayer(3, 'CCC', 10, density=1) #  Density initialized to 5g/cm^3

    result = sample.structure[1]

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
    print(sample.mag_elements)
    print()
    thickness, density, density_magnetic = sample.density_profile()

    # thickness, density, density_magnetic = sample.density_profile()
    val = list(density.values())

    mag_val = list(density_magnetic.values())
    plt.figure()
    for idx in range(len(val)):
        plt.plot(thickness, val[idx])


    for idx in range(len(mag_val)):
        plt.plot(thickness, -mag_val[idx],'--')

    center = np.zeros(len(thickness))
    plt.plot(thickness, center, 'k-.', linewidth=2)
    my_legend = list(density.keys())



    for key in list(density_magnetic.keys()):
        my_legend.append('Mag: ' + key)

    plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02,0.5))
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('Density (mol/cm^3)')
    plt.show()

    t = np.arange(-20,20,0.01)




