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

        if formula[0] == '(' and formula[1] == ')':  # Error Check
            raise NameError('Element symbol must be found inside parentheses')

        if formula[0] == '(':  # Checks for the start of parenthesis
            formula = formula[1:]  # Removes '(' from the list.
            ion_list = list()  # temporary list that contains the ion properties
            myions = list()  # List that contains all the related ions
            bracket = True  # Booleans that determines if the bracket is still open

        if bracket:  # Bracket is open
            if formula[0] == ')':  # Checks for a close bracket
                if len(formula) == 1:  # Case closed bracket is last character
                    digit = 1  # sets ion stoichiometry
                    formula = ''
                else:  # Case closed bracket in middle of string
                    if formula[1].isdigit():  # Checks stoichiometry of ion
                        digit = formula[1]
                        formula = formula[2:]
                    else:  # Stocihiometry is 1
                        digit = 1
                        formula = formula[1:]

                for ele in ion_list:  # Loops through all the ions and sets their properties to element class
                    mydict[ele[0]] = element(ele[0], digit, ion=myions, ionratio=ele[1]/total_ion)

                bracket = False  # closes bracket
            else:  # Bracket still open
                formula, info = checkstring(formula)  # Retrieves info
                ion_list.append(info)
                myions.append(info[0])
                total_ion = total_ion + info[1]

        else:  # None bracket case
            formula, info = checkstring(formula)
            mydict[info[0]] = element(info[0],info[1])

        n = len(formula)
    return mydict






class slab:
    def __init__(self, num_layers):
       
        self.structure = [dict() for i in range(num_layers)]  # Physical structure
        self.myelements = []  # Keeps track of the elements in the material

    def addlayer(self, num_layer, formula, thickness, density, roughness=0):
        

        # Retrieve the element info
        elements = getElements(formula)

        for key in list(elements.keys()):

            if key not in self.myelements:  # Adds element to my element list if not already there
                self.myelements.append(key)

            elements[key].density = density  # sets density  (mol/cm^3)
            elements[key].thickness = thickness  # sets thickness  (Angstrom)
            elements[key].roughness = roughness  # Order of Angstrom
        
        self.structure[num_layer] = elements  # sets the layer with the appropriate slab properties
        
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
                        dense[ele] = np.concatenate((dense[ele], np.ones(n)*layer[ele].density*layer[ele].stoich*layer[ele].ionratio))
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
                        dense[ele] = np.concatenate((dense[ele], np.ones(n)*layer[ele].density*layer[ele].stoich*layer[ele].ionratio))
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


    # Example showing how to create slab model of sample
    sample = slab(3)  # Initializing three layers
    sample.addlayer(0, 'SrTiO3', 50, 0.028)  # substrate layer
    sample.addlayer(1, 'La(MnFe)O3', 25, 0.01)  # Film 1 on top of substrate
    sample.addlayer(2, 'La(MnFe)O3', 16, 0.02)   # Film 2 on top film 1

    sample.showprofile()  # Showing the density profile

    molecule = 'La(MnFe)7O3'
    result = getElements(molecule)


    e = 'Fe'
    print('Name: ', result[e].name)
    print('Scattering Factor: ', result[e].sf)
    print('Stoichiometry: ',result[e].stoich)
    print('Ion: ',result[e].ion)
    print('Ion Ratio: ', result[e].ionratio)
    

