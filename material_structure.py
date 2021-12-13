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
    # Finds all the elements in the chemical formula
    myelement = re.findall(r'[A-Z][a-z]*|\d+', re.sub('[A-Z][a-z]*(?![\da-z])', r'\g<0>1', formula))

    # Find the ions in the formula which are identified with parentheses.
    ions = re.findall('[A-Z][a-z]?\d*|\([^()]*(?:\(.*\))?[^()]*\)\d+', formula)

    print(myelement)


    """
    # Separate the ions into their appropriate element symbol
    for ion in range(len(ions)):
        ions[ion] = re.findall('[A-Z][a-z]', ions[ion])

    stoich = dict()
    n = len(myelement)
    for idx in range(0,n,2):
        stoich[myelement[idx]] = element(myelement[idx],myelement[idx+1])

    for ionset in ions:
        for e in ionset:
            m = stoich[e].stoich
            print(m)
    """
    return ions




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
    molecule = 'La(Mn2Fe3)O3'
    result = getElements(molecule)
    print(result)

