import tkinter
from math import sin
import numpy as np
from scipy import interpolate
from scipy import optimize
import Pythonreflectivity as pr
from material_model import *


import tkinter as tk
from tkinter import ttk


def showSampleParameters(sample):

    parameters = ['Thickness', 'Density', 'Roughness', 'Linked Roughness',
                  'Stoichiometry', 'Polymorph', 'Gamma', 'Phi', 'Magnetic Density',
                  'Scattering Factor', 'Magnetic Scattering Factor']
    root = tk.Tk()
    root.title('Sample Information')
    root.geometry("900x900")
    tree = ttk.Treeview(root)
    ttk.Style().configure('Treeview', rowheight=30)

    structure = sample.structure

    # Loop through all the different elements
    for layer in range(len(structure)):
        if layer == 0:
            layer_name = 'Substrate'
        else:
            layer_name = 'Layer ' + str(layer)
        tree.insert("",layer, layer_name, text=layer_name)  # adding the layer to the tree

        elements = list(structure[layer].keys())

        for e in range(len(elements)):
            element = elements[e]
            element_name = element + " " + str(layer)
            tree.insert(layer_name,e,element_name,text=element)

            for param in parameters:
                if param == 'Thickness':
                    thick_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name,1,thick_name,text=param)
                    thick_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(thick_name, 1, thick_data_name, text=str(structure[layer][element].thickness) + " (A)")
                elif param == 'Density':
                    density_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 2, density_name, text=param)
                    density_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(density_name, 1, density_data_name, text=str(structure[layer][element].density) + " (mol/cm^3)")
                elif param == 'Roughness':
                    rough_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 3, rough_name, text=param)
                    rough_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(rough_name, 1, rough_data_name, text=str(structure[layer][element].roughness)) + " (A)"
                elif param == 'Linked Roughness':
                    link_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 4, link_name, text=param)
                    link_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(link_name, 1, link_data_name, text=str(structure[layer][element].linked_roughness)) + " (A)"
                elif param == 'Stoichiometry':
                    stoich_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 5, stoich_name, text=param)
                    stoich_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(stoich_name, 1, stoich_data_name, text=structure[layer][element].stoichiometry)
                elif param == 'Polymorph':
                    polymorphs = structure[layer][element].polymorph
                    if len(polymorphs) != 0:
                        polymorphs_name = param + " " + element + " " + str(layer)
                        tree.insert(element_name, 6, polymorphs_name, text=param)
                        for poly in range(len(polymorphs)):
                            poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                            tree.insert(polymorphs_name, poly, poly_name, text=polymorphs[poly])
                            poly_density_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                            poly_data = 'Ratio = ' + str(structure[layer][element].poly_ratio[poly])
                            tree.insert(poly_name, 1, poly_density_name, text=poly_data)
                elif param == 'Gamma':
                    gamma_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 7, gamma_name, text=param)
                    gamma_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(gamma_name, 1, gamma_data_name, text=str(structure[layer][element].gamma) + " degrees")
                elif param == 'Phi':
                    phi_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 8, phi_name, text=param)
                    phi_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(phi_name, 1, phi_data_name, text=str(structure[layer][element].phi) + " degrees")
                elif param == 'Magnetic Density':
                    polymorphs = structure[layer][element].polymorph
                    if len(polymorphs) != 0:
                        magnetic_name = param + " " + element + " " + str(layer)
                        tree.insert(element_name, 9, magnetic_name, text=param)
                        for poly in range(len(polymorphs)):
                            poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                            tree.insert(magnetic_name, poly, poly_name, text=polymorphs[poly])
                            mag_density_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                            poly_data = str(structure[layer][element].mag_density[poly]) + " mol/cm^3"
                            tree.insert(poly_name, 1, mag_density_name, text=poly_data)
                elif param == 'Scattering Factor':
                    polymorphs = structure[layer][element].polymorph
                    if len(polymorphs) != 0:
                        scat_name = param + " " + element + " " + str(layer)
                        tree.insert(element_name, 9, scat_name, text=param)
                        for poly in range(len(polymorphs)):
                            poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                            tree.insert(scat_name, poly, poly_name, text=polymorphs[poly])
                            scat_fact_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                            poly_data = str(structure[layer][element].scattering_factor[poly]) + " mol/cm^3"
                            tree.insert(poly_name, 1, scat_fact_name, text=poly_data)
                    else:
                        print('bye')
                elif param == 'Magnetic Scattering Factor':
                    scat_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 9, scat_name, text=param)
                    #scat_fact_name = param + " " + element +

    tree.pack(expand=True, fill='both')
    root.mainloop()

    """
    root = tk.Tk()
    root.title('Sample Information')
    root.geometry("900x900")
    tree = ttk.Treeview(root)
    ttk.Style().configure('Treeview', rowheight=30)

    def add_node(k, v):
        if isinstance(v, dict):
            for i, j in v.items():
                if tree.exists(j):
                    tree.insert(k, 1, i, text=i)
                    add_node(i, j)
                else:
                    tree.insert(k, 1, i, text=i)
                    add_node(i, j)

    for k, v in hierarchy.items():
        tree.insert("", 1, k, text=k)
        add_node(k, v)

    tree.pack(expand=True, fill='both')
    root.mainloop()
    """







if __name__ == "__main__":
    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density=[0.027904, 0.027904, 0.083712], roughness=[7.58207, False, 5.77093])
    sample.addlayer(1, 'SrTiO3', 6, density=[0, 0.027904, 0], roughness=[7.58207, 4.03102, 5.77093])

    sample.addlayer(2, 'LaMnO3', 4, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2],
                    linked_roughness=[False, 0.5, False])
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 17.8, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 9, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 2.5, density=[0.025, 0.024, 0.05], roughness=[0.25, 0.25, 2])
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 4.5, density=[0.025, 0.042, 0.04], roughness=[0.25, 0.25, 2])
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.0053, 0], ['Co', 'Ni'])

    sample.addlayer(7, 'CCO', 11.1, density=[0.05, 0.05, 0.01], roughness=2, linked_roughness=[3, 1.5, False])


    showSampleParameters(sample)




