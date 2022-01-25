import numpy as np
import matplotlib.pyplot as plt
import cmath
from scipy import special
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
        self.mag_dir = None  # The type of magnetization (isotropic and anisotropic)
        self.mag_density = None  # The scalling factor we want to multiply our scattering factor by (density is not the correct description)
        self.scattering_factor = name  # Identifies which scattering factor to be used. This parameter will allow us to implement 'scattering functions'
        self.mag_scattering_factor = None

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
        raise NameError(formula)

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
        if len(self.structure[num_layer]) != 0:
            warnings.warn('Layer '+str(num_layer)+' already initialized')

        # Thickness type check
        if type(thickness) != int and type(thickness) != float:
            raise TypeError('Thickness must be integer or float type')

         # Checks if thickness is in a reasonable range
        if thickness < 0:
            warnings.warn('Thickness must be a positive number')
            thickness = abs(thickness)
        elif thickness == 0:
            raise ValueError('Thickness cannot be zero')

        # Checks Density
        if density == None:
            pass
        elif type(density) != int and type(density) != float:
            raise TypeError('Density must be entered in as a float or integer type')
        elif density < 0:
            warnings.warn('Density must be positive')
            density = abs(density)
        elif density == 0:
            raise ValueError('The density of a material can not be zero')

        if density == None:
            pass
        elif density > 20:
            warnings.warn('The density of ' + str(density) + ' g/cm^3 might be too large of a value. Consider double checking your density. ')

        # Checks Roughness
        if type(roughness) != int and type(roughness) != float:
            raise TypeError('Roughness must be of float or integer type.')
        elif roughness < 0:
            roughness = abs(roughness)
            warnings.warn('Roughness should be entered as a positive value')

        if roughness > 15:
            warnings.warn('Roughness is much larger than expected')

        if link == None:
            pass
        elif type(link) != list:
            raise TypeError('Variable link must be a list')
        elif len(link) != num_elements:
            raise RuntimeError("Length of link must match number of elements in formula")
        else:
            print('y')

        # -------------------------------------------------------------------------------------------------------------#
        # Retrieve the element info
        elements, num_elements = find_stoichiometry(formula)
        molar_mass = 0
        # Computes the total molar mass of perovskite material
        for key in list(elements.keys()):
            molar_mass = molar_mass + atomic_mass(key)*elements[key].stoichiometry

        for key in list(elements.keys()):

            if key not in self.myelements:  # Adds element to my element list if not already there
                self.myelements.append(key)

            if density == None:
                density = perovskite_density(formula)

            elements[key].density = density  # sets density  (g/cm^3)
            elements[key].thickness = thickness  # sets thickness  (Angstrom)
            elements[key].roughness = roughness  # Order of Angstrom
            elements[key].molar_mass = molar_mass  # Molar mass of perovskite material

        self.structure[num_layer] = elements  # sets the layer with the appropriate slab properties

        # Determines if A/B site is linked to the A/B site on the next layer
        n = len(elements)
        if link == None:
            mylink = [True for i in range(num_elements)]
            self.link.append(mylink)
        else:
            self.link.append(link)

    def polymorphous(self, lay, ele, polymorph, poly_ratio, sf=None):

        poly_keys = list(self.poly_elements.keys())
        if ele not in poly_keys:
            self.poly_elements[ele] = polymorph

        self.structure[lay][ele].polymorph = polymorph
        self.structure[lay][ele].poly_ratio = poly_ratio
        if sf == None:
            self.structure[lay][ele].scattering_factor = polymorph
        else:
            self.structure[lay][ele].scattering_factor = sf

    def magnetization(self, lay, identifier, density, sf, mag_dir='z'):
        layer = self.structure[lay]
        if type(identifier) == list:



            for key in list(layer.keys()):
                if set(layer[key].polymorph) == set(identifier):
                    if layer[key].polymorph == identifier:
                        self.structure[lay][key].mag_scattering_factor = sf
                        self.structure[lay][key].mag_density = density
                        self.structure[lay][key].mag_dire = mag_dir
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
                        self.structure[lay][key].mag_type = mag_dir
                        self.mag_elements[key] = identifier
        else:

            self.structure[lay][identifier].mag_scattering_factor = sf
            self.structure[lay][identifier].mag_density = [density]
            self.structure[lay][identifier].mag_type = mag_dir
            self.mag_elements[identifier] = [identifier]


    def convert_to_slice(self):
            print('Time to convert')


    def density_profile(self):

        thick = np.array([])  # thickness array
        density_struct = {k:np.array([]) for k in self.myelements}
        density_poly = {k:dict() for k in list(self.poly_elements.keys())}
        density_mag = {k:dict() for k in list(self.mag_elements.keys())}

        for ele in list(self.poly_elements.keys()):
            density_struct.pop(ele)
            density_poly[ele] = {k:np.array([]) for k in self.poly_elements[ele]}

        for ele in list(self.mag_elements.keys()):
            density_mag[ele] = {k:np.array([]) for k in self.mag_elements[ele]}

        poly_keys = list(self.poly_elements.keys())

        first_layer = True  # Boolean that determines if using first layer
        step = 0.01  # Thickness step size (Angstrom)
        last = 0 + step  # Start of next slab

        for layer in self.structure:
            mykeys = list(layer.keys())  # retrieve elements in layer

            if first_layer:  # substrate layer
                temp = np.arange(-50,0+step,step)  # temporay thickness
                n = len(temp)  # length of the temporary thickness array
                thick = np.concatenate((thick, temp))  # Add slab thickness to material thickness
                first_layer = False

                # Takes care of the none polymorphous and none magnetic case
                for ele in list(density_struct.keys()):  # Create density array for each element
                    if ele in mykeys:  # Element in substrate layer
                        ratio = layer[ele].stoichiometry/layer[ele].molar_mass
                        density_struct[ele] = np.concatenate((density_struct[ele], np.ones(n)*layer[ele].density*ratio))
                    else:  # Element in film layer
                        density_struct[ele] = np.concatenate((density_struct[ele], np.zeros(n)))

                # Takes care of polymorphous case
                for ele in list(density_poly.keys()):
                    if ele in mykeys:  #Poly in substrate material
                        idx = 0
                        for poly in list(density_poly[ele].keys()):
                            ratio = layer[ele].stoichiometry*layer[ele].poly_ratio[idx]/layer[ele].molar_mass
                            density_poly[ele][poly] = np.concatenate((density_poly[ele][poly], np.ones(n)*layer[ele].density*ratio))
                            idx = idx + 1
                    else:
                        for poly in list(density_poly[ele].keys()):
                            density_poly[ele][poly] = np.concatenate((density_poly[ele][poly], np.zeros(n)))

                for ele in list(density_mag.keys()):
                    if ele in mykeys:  # Magnet in substrate material
                        idx = 0
                        for mag in list(density_mag[ele].keys()):
                            density_mag[ele][mag] = np.concatenate((density_mag[ele][mag], (-1)*np.ones(n)*layer[ele].mag_density[idx]))
                            idx = idx + 1
                    else:
                        for mag in list(density_mag[ele]):
                            density_mag[ele][mag] = np.concatenate((density_mag[ele][mag], np.zeros(n)))

            else:  # Film Layers
                first_element = list(layer.keys())[0]
                temp = np.arange(last, last + layer[first_element].thickness+step,step)
                n = len(temp)
                thick = np.concatenate((thick, temp))  # Add slab thickness to material thickness
                last = last + layer[first_element].thickness  # Update start of next layer

                # Takes care of the none polymorphous and none magnetic case
                for ele in list(density_struct.keys()):  # Create density array for each element
                    if ele in mykeys:  # Element in substrate layer
                        ratio = layer[ele].stoichiometry / layer[ele].molar_mass

                        density_struct[ele] = np.concatenate((density_struct[ele], np.ones(n) * layer[ele].density * ratio))
                    else:  # Element in film layer
                        density_struct[ele] = np.concatenate((density_struct[ele], np.zeros(n)))

                # Takes care of polymorphous case
                for ele in list(density_poly.keys()):
                    if ele in mykeys:  # Poly in substrate material
                        idx = 0
                        for poly in list(density_poly[ele].keys()):
                            ratio = layer[ele].stoichiometry * layer[ele].poly_ratio[idx]/layer[ele].molar_mass
                            density_poly[ele][poly] = np.concatenate((density_poly[ele][poly], np.ones(n) * layer[ele].density * ratio))
                            idx = idx + 1
                    else:
                        for poly in list(density_poly[ele].keys()):
                            density_poly[ele][poly] = np.concatenate((density_poly[ele][poly], np.zeros(n)))

                # Takes care of magnetic component
                for ele in list(density_mag.keys()):
                    if ele in mykeys:  # Magnet in substrate material
                        idx = 0
                        for mag in list(density_mag[ele].keys()):
                            density_mag[ele][mag] = np.concatenate((density_mag[ele][mag], (-1)*np.ones(n)*layer[ele].mag_density[idx]))
                            idx = idx + 1
                    else:
                        for mag in list(density_mag[ele]):
                            density_mag[ele][mag] = np.concatenate((density_mag[ele][mag], np.zeros(n)))

        # Create single dictionary to use
        density = density_struct
        for ele in list(density_poly.keys()):
            for poly in list(density_poly[ele].keys()):
                density[poly] = density_poly[ele][poly]

        density_magnetic = dict()
        for ele in list(density_mag.keys()):
            for mag in list(density_mag[ele].keys()):
                density_magnetic[mag] = density_mag[ele][mag]
        return thick, density, density_magnetic



if __name__ == "__main__":

    # Example: Simple sample creation
    sample = slab(4)  # Initializing four layers

    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, link=[False,True,True])  # substrate layer
    sample.polymorphous(0,'Ti',['Ti3+','Ti4+'], [0.25,0.75])

    # First Layer
    # Link: La-->La and O-->O
    # * denotes optional parameter
    sample.addlayer(1, 'LaMnO3', 25, link=[True,False,True])  # Film 1 on top of substrate
    sample.polymorphous(1,'Mn',['Mn3+','Mn2+'], [0.9,0.1], sf=['Mn','Fe'])  # (Layer, Element, Polymorph Symbols, Ratios, Scattering Factor)
    sample.magnetization(1, ['Mn3+', 'Mn2+'], [0.01,0.03], ['Ni', 'Co'])  # (Layer, Polymorph/Element, density, Scattering Factor, type*)

    # Second Layer
    # Link: La-->C and O-->C
    sample.addlayer(2, 'LaAlO3', 16, density=5, roughness=2, link=[True, False,True])   # Film 2 on top film 1
    sample.magnetization(2,'Al', 0.01, 'Co', mag_type='anisotropic') # mag_type is preset to 'isotropic


    # Impurity on surface
    # Impurity takes the form 'CCC'
    # The exact implementation of this has not quite been determined yet
    sample.addlayer(3, 'CCC', 10, density=1) #  Density initialized to 5g/cm^3


    result = sample.structure[1]

    e = 'Mn'
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
    print('Type: ', result[e].mag_type)

    #  print(sample.myelements)
    #  print(sample.poly_elements)
    print(sample.mag_elements)

    thickness, density, density_magnetic = sample.density_profile()
    val = list(density.values())
    mag_val = list(density_magnetic.values())
    plt.figure()
    for idx in range(len(val)):
        plt.plot(thickness, val[idx])
    for idx in range(len(mag_val)):
        plt.plot(thickness, mag_val[idx],'--')

    center = np.zeros(len(thickness))
    plt.plot(thickness, center, 'k-.', linewidth=2)
    my_legend = list(density.keys())
    for key in list(density_magnetic.keys()):
        my_legend.append('Mag: ' + key)
    plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02,0.5))
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('Density (mol/cm^3)')
    plt.show()



