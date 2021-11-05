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
        # Creates a dictionary with all the elements in the list
        element_dict = dict()
        profile_dict = dict()

        for el in elements_list:
            element_dict[el] = self.layer()
            profile_dict[el] = list()

        self.elements = element_dict
        self.profile = profile_dict

    def addbulk(self,element, thickness, density, roughness):
        self.elements[element].layers[0] = [thickness,density, roughness]

    def addlayer(self, element, thickness, density, roughness):
        self.elements[element].layers.append([thickness,density,roughness])
        self.elements[element].num_layers = self.elements[element].num_layers + 1

    def create_element_profile(self,element):
        layer_profile = self.elements[element].layers
        thickness, density = continuous(self.elements[element].layers)

        plt.figure()
        plt.plot(thickness,density)
        plt.show()

    def create_material_profile(self):
        plt.figure()
        for element in list(self.elements.keys()):
            thickness, density = continuous(self.elements[element].layers)
            plt.plot(thickness, density)

        plt.legend(list(self.elements.keys()))
        plt.show()



    def print_element(self,element):
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
        elements = list(self.elements.keys())
        for element in elements:
            self.print_element(element)
            print(" ")

    class layer:
        def __init__(self):
            self.layers = [[0,0,0]]
            self.num_layers = 0




if __name__ == "__main__":
    elements = ['Fe','O','La']
    mat = material(elements)

    mat.addbulk('Fe',0,10,10)
    mat.addlayer('Fe', 15,15,15)
    mat.addlayer('Fe',7,7,7)

    mat.addbulk('O', 0, 11, 10)
    mat.addlayer('O', 5, 16, 15)
    mat.addlayer('O', 7, 4, 7)

    mat.addbulk('La', 0, 4, 10)
    mat.addlayer('La', 16, 6, 15)
    mat.addlayer('La', 6, 9, 7)

    mat.create_material_profile()

