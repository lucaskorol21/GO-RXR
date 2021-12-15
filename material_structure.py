import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
import regex as re
import os
from tabulate import tabulate
from prettytable import PrettyTable

class element:
    def __init__(self, name, stoich, ion=None):
        """
        Purpose: Keep track of element properties in slab layer
        :param name: Element symbol
        :param thickness: Thickness of the slab (Angstrom)
        :param density: Density of the A/B-Site [Oxygen = 3X(A/B-Site)] mol/cm^3
        :param roughness: Roughness of the surface layer
        :param relation: N/A
        """
        self.name = name
        self.stoich = stoich
        self.density = 0
        self.thickness = 0
        self.roughness = 0
        self.ion = ion
        self.sf = name
        self.relations = None  # May not need


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

    :return: A dictionary that contains all the information
    """


    ele = None
    digit = 0
    n = len(formula)  # number of characters is 'formula'

    while n > 0:
        if n <= 2:  # Base Case
            if n == 2:  #  Case for Ul or UN
                if formula[1].isdigit():
                    ele = formula[0]
                    digit = formula[1]
                else:
                    ele = formula
                    digit = 1

                formula = ''
            else:
                ele = formula
                digit = 1
                formula = ''
            print(ele, digit)
        elif formula[0].isupper():
            if formula[1].islower():
                if formula[2].isdigit():
                    ele = formula[0:2]
                    digit = formula[2]
                    formula = formula[3:]
                    print(ele, digit)
                else:
                    ele = formula[0:2]
                    digit = 1
                    formula = formula[2:]
                    print(ele, digit)
            elif formula[1].isdigit():
                ele = formula[0]
                digit = formula[1]
                formula = formula[2:]
                print(ele, digit)
            else:
                ele = formula[0]
                digit = 1
                formula = formula[1:]
                print(ele, digit)
        else:
            raise NameError(formula)

        n = len(formula)
    return formula





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
    molecule = 'LaMn2Fe3O3'
    result = getElements(molecule)
    #print(result)

