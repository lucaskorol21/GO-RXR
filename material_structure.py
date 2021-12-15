import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
import regex as re
import os
from tabulate import tabulate
from prettytable import PrettyTable

class element:
    def __init__(self, name, stoich, ion=None, ionratio=1):
        """
        Purpose: Class that keeps track of elemental properties
        :param name: Elemental symbol
        :param stoich: Stoichiometry of each element
        :param ion:
        :param ionratio: Amount of ion that makes of element ---> (ion density) = (ionratio)*(element)
        """
        self.name = name  # Elemental Symbol
        self.density = 0  # Density of the molecule (mol/cm^3)
        self.thickness = 0  # Thickness of the layer (Angstrom)
        self.roughness = 0  # Roughness of the surface (Angstrom)
        self.stoich = stoich  # Stoichiometry of the element
        self.ionratio = ionratio  # Ion Ratio
        self.ion = ion  # The different ions of the same element
        self.sf = name  # The scattering factor name. This parameter will allow us to implement 'scattering functions'
        self.relations = None  # Can set relations if need be

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
            if formula[2].isdigit():  # Case ULD
                info.append(formula[0:2])
                info.append(int(formula[2]))
                formula = formula[3:]
            else:  # Case UL
                info.append(formula[0:2])
                info.append(1)
                formula = formula[2:]

        elif formula[1].isdigit():  # Case UD
            info.append(formula[0])
            info.append(int(formula[1]))
            formula = formula[2:]

        else:  # Case U
            info.append(formula[0])
            info.append(1)
            formula = formula[1:]
    else:  # Raises error if unexpected symbol
        raise NameError(formula)

    return formula, info

def getElements(formula):
    """
    Purpose: The purpose of this function is to identify all information from an inputted string.

    :param formula: The inputted string
                    has the form of a chemical formula with some extra notation. The element symbols identifies which
                    elements scattering factor you want to use, but can later be adjusted to use a function later.
                    Capability to use your own scattering factor textfile will also be provided, but how this will be
                    done has not yet been determined. The notation to describe the string is shown below:
                        U     - Uppercase letter
                        l     - Lowercase letter
                        N     - Number of atoms
                        (...) - Brackets is the notation we use to denote ions.
                        X     - Dummy variable that states if there is a linked roughness

                    Example:
                                formula = "AlMnO3"      ---> UlUlUN
                                formula = "La(MnFe)O3"  ---> Ul(UlUl)UN

                    Further description of notation.
                        Ions: The brackets state a single element with multiple ions,  where sum(ion density) = element density.
                            The elements in the inside corresponds to the scattering factor that best represents the ion. The
                            format of the ions is shown below:

                                                        (A_1N_1 A_2N_2 ... A_mN_m)

                            N_i - represents the ratio of the total density that is the ion A_i, for i = 1,2,...,m

    :return: A dictionary that has the element symbol as the key and the element class as the value {Symbol: element}
    """
    bracket = False  # Boolean that determines if dealing with bracket case
    mydict = dict()  # Dictionary [Symbol: element]
    ion_list = list()  # Keeps track of the ions in an element and their info
    myions = list()  # Keeps track of the ions in an element
    total_ion = 0  # total ion value
    n = len(formula)  # number of characters is 'formula'

    while n > 0:

        if formula[0] == '(':
            formula = formula[1:]
            ion_list = list()
            myions = list()
            bracket = True

        if bracket:
            if formula[0] == ')':
                if len(formula) == 1:
                    digit = 1
                    formula = ''
                else:
                    if formula[1].isdigit():
                        digit = formula[1]
                        formula = formula[2:]
                    else:
                        digit = 1
                        formula = formula[1:]

                for ele in ion_list:
                    mydict[ele[0]] = element(ele[0], digit, ion=myions, ionratio=ele[1]/total_ion)

                bracket = False
            else:
                formula, info = checkstring(formula)
                ion_list.append(info)
                myions.append(info[0])
                total_ion = total_ion + info[1]

        else:
            formula, info = checkstring(formula)
            mydict[info[0]] = element(info[0],info[1])

        n = len(formula)
    return mydict





"""
class slab:
    def __init__(self, num_layers):
       
        self.structure = [dict() for i in range(num_layers)]  # Physical structure
        self.myelements = []  # Keeps track of the elements in the material

    def addlayer(self, num_layer, layer_info):
        

        # Retrieve the element info
        elements = get_elements(layer_info[0])
        thickness = layer_info[1]
        density = layer_info[2]

        # set the layer properties
        info = dict()
        for key in list(elements.keys()):
            info[key] = element(key,thickness, density*float(elements[key]))
            if not(key in self.myelements):
                self.myelements.append(key)

        self.structure[num_layer] = info

    def addIon(self, identifier, ion1, ion2, roughness=None, relation=None):
        

        for layer in self.structure:
            if identifier in layer:  # determine if ion is found in the layer
                density = layer[identifier].density
                thickness = layer[identifier].thickness
                rough = layer[identifier].roughness
                layer.pop(identifier)  # remove identifier from layer
                if relation == None:
                    layer[ion1] = element(ion1,thickness,density*0.5,roughness=rough)  # add ion1 to layer
                    layer[ion2] = element(ion2, thickness, density*0.5, roughness=rough)  # add ion2 to layer
                else:
                    layer[ion1] = element(ion1, thickness, density * relation, roughness=rough)  # add ion1 to layer
                    layer[ion2] = element(ion2, thickness, density * (1-relation), roughness=rough)  # add ion2 to layer

        self.myelements.pop(self.myelements.index(identifier))  # remove identifier from list containing elements in material
        self.myelements.append(ion1)  # add ion1 to list
        self.myelements.append(ion2)  # add ion2 to list

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
                        dense[ele] = np.concatenate((dense[ele], np.ones(n)*layer[ele].density))
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
                        dense[ele] = np.concatenate((dense[ele], np.ones(n)*layer[ele].density))
                    else:  # Element not found in current layer
                        dense[ele] = np.concatenate((dense[ele], np.zeros(n)))

        # Create plot
        mydense = np.transpose(list(dense.values()))
        plt.figure()
        plt.plot(thick, mydense)
        plt.show()

"""


if __name__ == "__main__":

    """
    # Example showing how to create slab model of sample
    sample = slab(3)  # Initializing three layers
    sample.addlayer(0, ['SrTiO3', 50, 0.028])  # substrate layer
    sample.addlayer(1, ['LaMnO3', 25, 0.005])  # Film 1 on top of substrate
    sample.addlayer(2, ['LaMnO3', 16, 0.02])   # Film 2 on top film 1
    sample.addIon('Mn','Mn','Fe')  # Create ions for Manganese
    sample.showprofile()  # Showing the density profile
    """
    molecule = 'La(MnFe)7O3'
    result = getElements(molecule)

    e = 'Fe'
    print('Name: ', result[e].name)
    print('Scattering Factor: ', result[e].sf)
    print('Stoichiometry: ',result[e].stoich)
    print('Ion: ',result[e].ion)
    print('Ion Ratio: ', result[e].ionratio)

