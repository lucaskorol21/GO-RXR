import os
import sys
# Add the parent directory of GUI_GO.py to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import importlib
from collections import Counter
import UTILS.material_structure as ms
import pickle
import unittest
importlib.reload(ms)
import numpy as np

# This test script can be executed by inputting
#  ->  python -m unittest -v Testing/density_testing.py
# into the terminal

class TestDensityProfile(unittest.TestCase):
    
    def test_simple_profile(self):
        # Test a simple density profile

        # create the test model
        sample = ms.slab(5)
        sample.addlayer(0,'SrTiO3', 50, density=[0.028,0.028,0.084], roughness=[1.5,2,2.5])
        sample.addlayer(1, 'LaMnO3', 10, density=[0.028,0.028,0.084], roughness=[0, 1, 3])
        sample.addlayer(2, 'CaTiO3', 50, density=[0.028,0.028,0.084], roughness=[5, 2, 0.5])
        sample.addlayer(3, 'LaAlO3', 30, density=[0.028,0.028,0.084], roughness=[1.5, 2, 2.5])
        sample.addlayer(4, 'AlYO2', 50, density=4.5, roughness=[3, 2, 0.77])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness

        filename = 'simple_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)
            
        print('test_case', test_case)
        print('solution', solution)
        print(len(test_case), len(solution))

        # len(loaded_data)
        #Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())
        print('total', total)

        # self.assertEqual(total, 0)  # checks to make sure that the values are correct
        # Assert that the total is close to 0 with an epsilon difference
        self.assertAlmostEqual(total, 0, delta=np.finfo(float).eps)

    def test_elementVar(self):
        # Testing models that include element variations (defined as polymorphous in the code)

        # create the test model
        sample = ms.slab(2)
        sample.addlayer(0, 'ABO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
        sample.addlayer(1, 'ABO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])

        sample.polymorphous(0, 'A', ['Sr','La'], [0,1], sf=['Sr','La'])
        sample.polymorphous(0, 'B', ['Ti', 'Mn2', 'Mn3'], [1,0,0], sf=['Ti','Fe','Mn'])
        sample.polymorphous(1, 'A', ['Sr', 'La'], [0.25, 0.75], sf=['Sr','La'])
        sample.polymorphous(1, 'B', ['Ti', 'Mn2', 'Mn3'], [0.2, 0.35, 0.45], sf=['Ti', 'Fe', 'Mn'])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness

        filename = 'var_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

    def test_surfaceImpurity(self):
        # testing to make sure that the surface layer is being defined properly
        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
        sample.addlayer(1, 'LaMnO3', 10, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])
        sample.addlayer(2, 'CCO', [15, 10, 4], density=[0, 0.08, 0.05], roughness=[2, 3, 0.75])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness

        filename = 'surfaceImpurity_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

    def test_unitCell(self):
        # create the test model
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

        filename = 'unitCell_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

    def test_linkedRoughness(self):
        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness

        filename = 'linked_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

    def test_magnetic(self):
        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

        sample.magnetization(1, ['Mn2', 'Mn3'], [0.015, 0.01], ['Co', 'Ni'])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness
        for key in mag_density:
            test_case['Mag:' + key] = mag_density[key]

        filename = 'mag_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

    def test_dummy(self):
        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnX0', 15, density=4, roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiQ0', 20, density=3.5, roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness

        filename = 'dummy_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

    def test_negative_test(self):
        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[-0.028, 0.028, 0.084], roughness=[1.5, 5, -2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, -0.084], roughness=[0, -1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, -0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, -4.5, False])

        sample.magnetization(1, ['Mn'], [-0.015], ['Co'])

        thickness, density, mag_density = sample.density_profile()

        test_case = density
        test_case['Thickness'] = thickness

        filename = 'negative_density.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        with open(my_path, 'rb') as file:
            solution = pickle.load(file)

        # len(loaded_data)
        # Counter(loaded_data)
        self.assertEqual(len(test_case), len(solution))  # make sure there are the same number of elements
        self.assertEqual(Counter(test_case.keys()), Counter(solution.keys()))  # checks that all keys are the same

        total = sum(sum(abs(test_case[key]-solution[key])) for key in test_case.keys())

        self.assertEqual(total,0)  # checks to make sure that the values are correct

if __name__ == '__main__':
    unittest.main()
