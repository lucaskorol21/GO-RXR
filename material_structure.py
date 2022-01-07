import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *

class element:
    def __init__(self, name, stoichiometry):
        """
        Purpose: Class that keeps track of elemental properties
        :param name: Elemental symbol
        :param stoichiometry: Stoichiometry of each element
        :param polymorph:
        :param poly_ratio: Amount of ion that makes of element ---> (ion density) = (ionratio)*(element)
        """
        self.name = name  # Elemental Symbol
        self.molar_mass = 0
        self.density = 0  # Density of the molecule (g/cm^3)
        self.thickness = 0  # Thickness of the layer (Angstrom)
        self.roughness = 0  # Roughness of the surface (Angstrom)
        self.stoichiometry = stoichiometry  # Stoichiometry of the element
        self.poly_ratio = 1  # Ratio of the total density that this polymorphous form makes up of the total element.
        self.polymorph = []  # A list that contains the various forms (e.g. ions)
        self.scattering_factor = name  # Identifies which scattering factor to be used. This parameter will allow us to implement 'scattering functions'

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
            else:  # Case UL
                info.append(formula)
                info.append(1)

            formula = ''  # Sets formula for end case
        else:  # Case U
            info[0].append(formula)
            info[1].append(1)

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
    mydict = dict()
    while len(formula) > 0:
        formula, info = checkstring(formula)
        mydict[info[0]] = element(info[0], info[1])

    return mydict

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

        self.structure = [dict() for i in range(num_layers)]  # Physical structure
        self.link = [[True, True] for i in range(num_layers)]  # [A-site, B-site]]
        self.myelements = []  # Keeps track of the elements in the material

    def addlayer(self, num_layer, formula, thickness, density=None, roughness=0, A_site=True, B_site=True):
        """
        Purpose: Add layer to sample material
        :param num_layer: Layer number (0 indicates substrate layer)
        :param formula: Perovskite chemical formula
        :param thickness: Thickness of the material in Angstrom 10^(-10)m
        :param density: Density of the perovskite material (g/cm^3). No input indicates user wants to use density found in database
        :param roughness: Rougness of at the interface
        :param A_site: Boolean value
                        True - Indicates roughness of A-site from current layer is linked to roughness of A-site in next layer
                        False - Indicates roughness of A-site from current layer is NOT linked to roughness of A-site in next layer
        :param B_site: Boolean value
                        True - Indicates roughness of B-site from current layer is linked to roughness of B-site in next layer
                        False - Indicates roughness of B-site from current layer is NOT linked to roughness of B-site in next layer
        """


        # Retrieve the element info
        elements = find_stoichiometry(formula)
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
        if not(A_site):
            self.link[num_layer][0] = False
        if not(B_site):
            self.link[num_layer][1] = False

    def polymorphous(self, lay, ele, polymorph, poly_ratio, sf=None):
        self.structure[lay][ele].polymorph = polymorph
        self.structure[lay][ele].poly_ratio = poly_ratio
        if sf == None:
            self.structure[lay][ele].scattering_factor = polymorph
        else:
            self.structure[lay][ele].scattering_factor = sf


    def setsigma(self, el, lay, rough):
       
        self.structure[lay][el].roughness = rough

    def setd(self, el, lay, thick):
        self.structure[lay][element].thickness = thick

    def setdensity(self, el, lay, rho):
        
        self.structure[lay][el].density = rho

    def convert_to_slice(self):
            print('Time to convert')


    def showprofile(self):
        thick = np.array([])  # thickness array
        dense = {k:np.array([]) for k in self.myelements}  # Density dictionary
        first_layer = True  # Boolean that determines if using first layer
        step = 0.1  # Thickness step size (Angstrom)
        last = 0 + step  # Start of next slab

        for layer in self.structure:
            mykeys = list(layer.keys())  # retrieve elements in layer

            if first_layer:  # substrate layer
                temp = np.arange(-50,0+step,step)
                n = len(temp)
                thick = np.concatenate((thick, temp))  # Add slab thickness to material thickness
                first_layer = False
                for ele in self.myelements:  # Create density array for each element
                    if ele in mykeys:  # Element in substrate layer
                        dense[ele] = np.concatenate((dense[ele], np.ones(n)*layer[ele].density*layer[ele].stoichiometry/layer[ele].molar_mass))
                    else:  # Element in film layer
                        dense[ele] = np.concatenate((dense[ele], np.zeros(n)))
            else:  # Film Layers
                first_element = list(layer.keys())[0]
                temp = np.arange(last, last + layer[first_element].thickness+step,step)
                print(temp)
                n = len(temp)
                thick = np.concatenate((thick, temp))  # Add slab thickness to material thickness
                last = last + layer[first_element].thickness  # Update start of next layer
                for ele in self.myelements:  # Create density array for each element
                    if ele in mykeys:  # Element in current layer
                        dense[ele] = np.concatenate((dense[ele], np.ones(n)*layer[ele].density*layer[ele].stoichiometry/layer[ele].molar_mass))
                    else:  # Element not found in current layer
                        dense[ele] = np.concatenate((dense[ele], np.zeros(n)))

        # Create plot
        mydense = np.transpose(list(dense.values()))
        plt.figure()
        plt.plot(thick, mydense)
        plt.xlabel('Thickness (Angstrom)')
        plt.ylabel('Density (mol/cm^3)')
        plt.suptitle('Density Profile')
        plt.show()



if __name__ == "__main__":

    # Example: Simple sample creation
    sample = slab(3)  # Initializing three layers

    # Sr roughness not linked to La
    sample.addlayer(0, 'SrTiO3', 50, A_site=False)  # substrate layer

    # Mn roughness not linked to Al
    sample.addlayer(1, 'LaMnO3', 25, B_site=False)  # Film 1 on top of substrate
    sample.polymorphous(1,'Mn',['Mn3+','Mn2+'], [0.5,0.5], sf=['Mn','Fe'])  # [Layer, Element, Polymorph Symbols, Ratios, Scattering Factor]

    # Showing other layer properties
    sample.addlayer(2, 'LaAlO3', 16, density=8.08, roughness=2)   # Film 2 on top film 1

    sample.showprofile()  # Showing the density profile

    result = sample.structure[1]
    e = 'La'
    print('Name: ', result[e].name)
    print('Scattering Factor: ', result[e].scattering_factor)
    print('Stoichiometry: ',result[e].stoichiometry)
    print('Polymorph: ',result[e].polymorph)
    print('Polymorph Ratio : ', result[e].poly_ratio)




