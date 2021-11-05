import numpy
import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
from prettytable import PrettyTable
import scipy.special as special


def continuous(element):
    thick = [i[0] for i in element]
    dense = [i[1] for i in element]
    rough = [i[2] for i in element]

    n = 100000

    thickness = []
    density = []

    extra = 10
    total_thick = sum(thick) + 30 + extra
    interval = n / total_thick
    begin = -30
    end = 0
    for idx in range(len(thick)):
        if idx == 0:
            begin = -30
            end = 0
        else:
            begin = end
            end = end + thick[idx]

        steps = int((end - begin) * interval)
        temp_thick = list(np.linspace(begin, end, steps))
        temp_dense = list(np.ones(steps) * dense[idx])

        thickness = thickness + temp_thick
        density = density + temp_dense

    steps = int((extra) * interval)
    temp_thick = list(np.linspace(end, end + extra, steps))
    temp_dense = list(np.zeros(steps))

    thickness = thickness + temp_thick
    density = density + temp_dense
    return thickness, density

class material:
    def __init__(self,elements_list):
        """
        Purpose: Initialize the material class
        :param elements_list: A list of strings that contain the element symbols for all elements in the material
        """
        # Creates a dictionary with all the elements in the list
        element_dict = dict()  # dictionary that contains the element and layer information
        profile_dict = dict()  # contain the continuous depth profile of each element

        for el in elements_list:
            element_dict[el] = self.layer()  # for each element create a layer class type
            profile_dict[el] = list()

        self.elements = element_dict
        self.profile = profile_dict

    def addbulk(self,element, thickness, density, roughness):
        """
        Purpose: Add a bulk layer to the layer class for a specified element
        :param element: The element symbol
        :param thickness: Thickness of the layer (Angstrom)
        :param density: Density of the layer (mol/cm^3)
        :param roughness: Roughness at the layers surface (Angstrom)
        """
        self.elements[element].layers[0] = [thickness,density, roughness]

    def addlayer(self, element, thickness, density, roughness):
        """
        Purpose: Add a layer on the surface for a specified material
        :param element: Element symbol
        :param thickness: Layer thickness (Angstrom)
        :param density: Layer density (mol/cm^3)
        :param roughness: Layer roughness (Angstrom)
        :return:
        """
        self.elements[element].layers.append([thickness,density,roughness])
        self.elements[element].num_layers = self.elements[element].num_layers + 1

    def create_element_profile(self,element):
        """
        Purpose: Creates the continuous density profile of an element
        :param element: Element symbol
        """
        layer_profile = self.elements[element].layers  # gets all the layers of the element
        thickness, density = continuous(self.elements[element].layers)  # creates the continuous profile from the layers

        #  Plotting the results
        plt.figure()
        plt.plot(thickness,density)
        plt.show()

    def create_material_profile(self):
        """
        Purpose: Show the material density profile for all elements
        """
        plt.figure()
        for element in list(self.elements.keys()):
            thickness, density = continuous(self.elements[element].layers)
            plt.plot(thickness, density)

        plt.legend(list(self.elements.keys()))
        plt.show()



    def print_element(self,element):
        """
        Purpose: Prints out the depth profile parameters (thickness, density, roughness)
        :param element: The desired elements symbol
        """
        my_layers = self.elements[element].layers
        t = PrettyTable(["Layer","Thickness (A)","Density (mol/cm^3)","Roughness (A)"])
        idx = 0
        for layer in my_layers:
            if idx == 0:
                t.add_row(["Bulk " + str(idx), layer[0], layer[1], layer[2]])
            else:
                t.add_row(["Layer " + str(idx), layer[0], layer[1], layer[2]])
            idx = idx + 1
        print("                             " + element)
        print(t)

    def print_elements(self):
        """
        Purpose: Prints the depth profile parameters of all elements
        """
        elements = list(self.elements.keys())
        for element in elements:
            self.print_element(element)
            print(" ")

    class layer:
        def __init__(self):
            self.layers = [[0,0,0]]
            self.num_layers = 0




if __name__ == "__main__":

    # Example Material
    elements = ['Fe','O','La']  # Input all the elements in the material in a list
    mat = material(elements)  # create the material class

    # Input the layer parameters (thickness, density, roughness)
    mat.addbulk('Fe',0,10,10)  # In general the bulk will not have a thickness
    mat.addlayer('Fe', 15,15,15)  # Layer on top of bulk
    mat.addlayer('Fe',7,7,7)  # Layer on top of previous layer

    # Input the layer parameters
    mat.addbulk('O', 0, 11, 10)
    mat.addlayer('O', 5, 16, 15)
    mat.addlayer('O', 7, 4, 7)

    # Input the layer parameters
    mat.addbulk('La', 0, 4, 10)
    mat.addlayer('La', 16, 6, 15)
    mat.addlayer('La', 6, 9, 7)

    mat.create_material_profile()

