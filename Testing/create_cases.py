import os
import pickle
import material_structure as ms
import matplotlib.pyplot as plt
import material_model as mm
import numpy as np
import data_structure as ds


def saveTest(data, fname):

    my_dir = os.getcwd()

    filename = my_dir + '/test_data/' + fname
    with open(filename, 'wb') as file:
        pickle.dump(data, file)


def simple_model():
    sample = ms.slab(5)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(1, 'LaMnO3', 10, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])
    sample.addlayer(2, 'CaTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[5, 2, 0.5])
    sample.addlayer(3, 'LaAlO3', 30, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(4, 'AlYO2', 50, density=4.5, roughness=[3, 2, 0.77])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'simple_density.pkl')

def var_model():
    sample = ms.slab(2)
    sample.addlayer(0, 'ABO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(1, 'ABO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])

    sample.polymorphous(0, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(0, 'B', ['Ti', 'Mn2', 'Mn3'], [1, 0, 0], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(1, 'A', ['Sr', 'La'], [0.25, 0.75], sf=['Sr', 'La'])
    sample.polymorphous(1, 'B', ['Ti', 'Mn2', 'Mn3'], [0.2, 0.35, 0.45], sf=['Ti', 'Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'var_density.pkl')

def surfaceImpurity_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(1, 'LaMnO3', 10, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])
    sample.addlayer(2, 'CCO', [15, 10, 4], density=[0, 0.08, 0.05], roughness=[2, 3, 0.75])

    thickness, density, mag_density = sample.density_profile()

    test_case = density
    test_case['Thickness'] = thickness

    saveTest(test_case, 'surfaceImpurity_density.pkl')

def unitCell_model():
    sample = ms.slab(6)
    sample.addlayer(0, 'ABO3', 50, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(1, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(2, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(3, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(4, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(5, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)

    sample.polymorphous(0, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(0, 'B', ['Ti', 'Mn2', 'Mn3'], [1, 0, 0], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(1, 'A', ['Sr', 'La'], [0.75, 0.25], sf=['Sr', 'La'])
    sample.polymorphous(1, 'B', ['Ti', 'Mn2', 'Mn3'], [0.2, 0.35, 0.45], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(2, 'A', ['Sr', 'La'], [0.5, 0.5], sf=['Sr', 'La'])
    sample.polymorphous(2, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0.25, 0.75], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(3, 'A', ['Sr', 'La'], [0.25, 0.75], sf=['Sr', 'La'])
    sample.polymorphous(3, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0, 1], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(4, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(4, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0, 1], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(5, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(5, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0, 1], sf=['Ti', 'Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    test_case = density
    test_case['Thickness'] = thickness

    saveTest(test_case, 'unitCell_density.pkl')

def linkedRough_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],linked_roughness=[False,0.25,False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5], linked_roughness=[4,4.5,False])

    sample.polymorphous(1, 'Mn', ['Mn2','Mn3'], [0.5,0.5], sf=['Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'linked_density.pkl')

def magnetic_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, 4.5, False])

    sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

    sample.magnetization(1, ['Mn2','Mn3'], [0.015,0.01], ['Co','Ni'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness
    for key in mag_density:
        data['Mag:'+key] = mag_density[key]

    saveTest(data, 'mag_density.pkl')

def dummy_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnX0', 15, density=4, roughness=[0, 1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiQ0', 20, density=3.5, roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, 4.5, False])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'dummy_density.pkl')

def negative_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[-0.028, 0.028, 0.084], roughness=[1.5, 5, -2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, -0.084], roughness=[0, -1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, -0.028, 0.084], roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, -4.5, False])

    sample.magnetization(1, ['Mn'], [-0.015], ['Co'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'negative_density.pkl')


if __name__ == "__main__":

    filename = '7uc_sample_test.h5'
    if os.getcwd().split('\\')[-1] == 'Testing':
        my_path = os.getcwd() + '/test_data/' + filename
    else:
        my_path = os.getcwd() + '/Testing/test_data/' + filename

    sample = ds.ReadSampleHDF5(my_path)


    parameters = [[0, 'B', 'Ti', 'Mn2', 0.5], [1, 'B', 'Ti', 'Mn3', 0.1], [2, 'B', 'Mn3', 'Mn2', 0.9],
                  [3, 'B', 'Mn2', 'Mn3', 0.75], [4, 'B', 'Ti', 'Mn2', 0.55], [5, 'B', 'Mn2', 'Ti', 0.925],
                  [6, 'B', 'Mn3', 'Ti', 0.4], [7, 'B', 'Ti', 'Mn2', 0.33333337], [8, 'B', 'Ti', 'Mn3', 0.84],
                  [9, 'B', 'Mn3', 'Mn2', 0.65315315131], [00, 'B', 'Mn2', 'Mn3', 0.22232]]

    #layer, symbol, identifier1, identifier2, ratio
    # Ti, Mn2, Mn3

    my_list = []
    my_dict = {}
    for params in parameters:
        my_dict = {}
        sample.setRatio(params[0], params[1], params[2], params[3], params[4])
        for i, key in enumerate(sample.structure[params[0]]['B'].polymorph):
            my_dict[key] = sample.structure[params[0]]['B'].poly_ratio[i]
        my_list.append(my_dict)
    print(my_list)