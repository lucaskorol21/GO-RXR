import numpy
import numpy as np
import matplotlib.pyplot as plt
import cmath
from KK_And_Merge import *
import os
from tabulate import tabulate
from prettytable import PrettyTable


class material:
    def __init__(self,elements_list):
        # Creates a dictionary with all the elements in the list
        element_dict = dict()
        profile_dict = dict()
        print('helloIniti')

        for el in elements_list:
            element_dict[el] = self.layer()
            profile_dict[el] = list()

        self.elements = element_dict
        self.profile = profile_dict

    def addbulk(self,element, thickness, density, roughness):
        self.elements[element].bulk = [thickness,density, roughness]

    def addlayer(self, element, thickness, density, roughness):
        self.elements[element].layers.append([thickness,density,roughness])
        self.elements[element].num_layers = self.elements[element].num_layers + 1

    def create_density_profile(self):
        print('hello')

    def print_element(self,element):
        my_layers = self.elements[element].layers
        t = PrettyTable(["Layer","Thickness (A)","Density (mol/cm^3)","Roughness (A)"])
        idx = len(my_layers)
        for layer in reversed(my_layers):
            t.add_row(["Layer " + str(idx), layer[0], layer[1], layer[2]])
            idx = idx - 1

        bul = self.elements[element].bulk
        t.add_row(["Bulk", bul[0], bul[1], bul[2]])
        print("                             " + element)
        print(t)

    def print_elements(self):
        elements = list(self.elements.keys())
        for element in elements:
            self.print_element(element)
            print(" ")

    class layer:
        def __init__(self):
            self.bulk = [0,0,0]
            self.layers = list()
            self.num_layers = 0




if __name__ == "__main__":
    elements = ['Fe','O','La']
    mat = material(elements)

    mat.addbulk('Fe',0,0.028,2)
    mat.addlayer('Fe', 15,0.03,2)
    mat.addlayer('Fe',2,0.01,2)

    mat.addbulk('O', 0, 0, 2)
    mat.addlayer('O', 30, 0.05, 2)
    mat.addlayer('O', 2, 0.05, 2)

    mat.addbulk('La', 0, 0, 0)
    mat.addlayer('La', 10, 0.028, 2)
    mat.addlayer('La', 7, 0.028, 2)
    mat.print_elements()

