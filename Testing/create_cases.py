import os
import pickle
import material_structure as ms
import matplotlib.pyplot as plt


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
    sample.polymorphous(5, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0.1, 0.9], sf=['Ti', 'Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    test_case = density
    test_case['Thickness'] = thickness

    saveTest(test_case, 'unitCell_density.pkl')

if __name__ == "__main__":
    unitCell_model()
    print('DONE!')