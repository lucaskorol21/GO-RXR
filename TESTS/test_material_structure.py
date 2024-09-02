import os
import sys

# Get the parent directory of the current script's directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# Add the parent directory to the system path
sys.path.append(parent_dir)

import unittest
import numpy as np
import UTILS.material_structure as ms
import pickle
import UTILS.data_structure as ds
from UTILS import TESTS_DIR

# This test script can be executed by inputting
#  ->  python -m unittest -v test_material_structure.py
# into the terminal

class TestMaterialStructure(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        self.tests_dir = TESTS_DIR

        self.filename1 = 'optical_energy.txt'
        self.my_path1 = self.tests_dir + '/test_data/' + self.filename1 

        self.filename2 = 'ALS_test.pkl'
        self.my_path2 = self.tests_dir + '/test_data/' + self.filename2

        self.filename3 = '7uc_sample_test.h5'
        self.my_path3 = self.tests_dir + '/test_data/' + self.filename3

        self.filename4 = 'LSMO_test.h5'
        self.my_path4 = self.tests_dir + '/test_data/' + self.filename4


    def test_adaptive_layer_segmentation(self):
        
        # Test the adaptive layer segmentation implementation
        optical = np.loadtxt(self.my_path1)

        with open(self.my_path2, 'rb') as file:
            solution = pickle.load(file)

        precision = [1e-1, 1e-2, 1e-5, 1e-10, 1e-20]


        for i in range(4):
            delta = optical[i, :]
            beta = optical[i + 4, :]
            delta_m = optical[i + 8, :]
            beta_m = optical[i + 12, :]

            n = 1 + np.vectorize(complex)(-delta, beta)  # complex index of refraction
            epsilon = n ** 2  # dielectric constant computation

            # definition as described in Lott Dieter Thesis
            Q = beta_m + 1j * delta_m  # magneto-optical constant

            for prec in precision:
                my_slabs = ms.ALS(epsilon.real, epsilon.imag, Q.real, Q.imag, prec).astype(int)
                my_slabs = [idx for idx in my_slabs]  # transforms numpy array to list
                compare_list = [idx for idx in solution[prec][i]]  # transforms numpy array to list
                self.assertListEqual(my_slabs, compare_list)

    def test_checkstring(self):
        # Testing the function that retrieves the next element an its stoichiometric relation
        #  - this function can only deal with integer stoichiometries

        tests = ['LaMnO3', 'C7C10C50C100LaMnO3','He132','J1K2E3','M84J8L5', 'La5ULaMnOKBe','Cr1Mn1O3']
        stoich = [['La',1],['C',7],['He',132],['J',1],['M',84], ['La', 5], ['Cr',1]]
        remaining_characters = ['MnO3', 'C10C50C100LaMnO3','','K2E3','J8L5','ULaMnOKBe','Mn1O3']

        for idx,test in enumerate(tests):
            test_rs,test_stoich = ms.checkstring(test)

            self.assertListEqual(stoich[idx],test_stoich)
            self.assertEqual(remaining_characters[idx],test_rs)

    def test_find_stoichiometry(self):
        # Testing to make sure that the proper stoichiometry is set
        tests = ['LaMnO3', 'C7C10C50C100LaMnO3', 'He132', 'J1K2E3', 'M84J8L5', 'La5ULaMnOKBe', 'Cr1Mn1O3']
        number_elements = [3,7,1,3,3,7,3]
        formula_keys = [['La','Mn','O'],['C1','C2','C3','C4','La','Mn','O'], ['He'],['J','K','E'],['M','J','L'],
                        ['La1','U','La2','Mn','O','K','Be'],['Cr','Mn','O']]
        stoichiometry = [[1,1,3],[7,10,50,100,1,1,3],[132],[1,2,3],[84,8,5],[5,1,1,1,1,1,1],[1,1,3]]



        for idx, test in enumerate(tests):
            formula, info = ms.find_stoichiometry(test)

            self.assertEqual(info,number_elements[idx])
            self.assertEqual(list(formula.keys()), formula_keys[idx])

            for i,ele in enumerate(formula_keys[idx]):
                self.assertEqual(formula[ele].stoichiometry,stoichiometry[idx][i])

    def test_get_roughness(self):

        # Testing the getRoughness function used in the scrip
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[10, 'A'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [11, 'all']]
        roughness = [0, 0.5, 3.726885630758733, 2.0653473900870867, 0.0, 0.0, 0.0, 2]

        for idx, params in enumerate(parameters):
            sigma = sample.getRoughness(params[0],params[1])
            self.assertEqual(sigma, roughness[idx])

    def test_set_roughness(self):

        # Testing the setRoughness function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        # retrieves initial configuration of the sample definition
        info = []
        for layer in sample.structure:
            my_info = {}
            for ele in layer.keys():
                my_info[ele] = layer[ele].roughness
            info.append(my_info)

        # new roughness values
        parameters = [[11, 'C2'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [9, 'all']]
        roughness = [2.5, 1, 3.25, 4.8118, 1.2, 5, 0.5, 0.984]

        for idx, params in enumerate(parameters):
            sample.setRoughness(params[0], params[1], roughness[idx])

            # manually changes the parameters
            if params[1] == 'all':
                for key in info[params[0]].keys():
                    info[params[0]][key] = roughness[idx]
            else:
                info[params[0]][params[1]] = roughness[idx]

            # checks to make sure that all roughness parameters have been altered properly
            for i, layer in enumerate(info):
                for ele in layer.keys():
                    self.assertEqual(layer[ele], sample.structure[i][ele].roughness)

    def test_get_density(self):
        
        # Testing the setDensity function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[10, 'A'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [8, 'all']]
        density = [0.0020259688666223116, 0.084, 0.08597982010497454, 0.046124560151674036, 0.028, 0.028, 0.084, 0.028]

        for idx, params in enumerate(parameters):
            rho = sample.getDensity(params[0], params[1])
            self.assertEqual(rho, density[idx])

    def test_set_density(self):

        # Testing the setDensity function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        # retrieves initial configuration of the sample definition
        info = []
        for layer in sample.structure:
            my_info = {}
            for ele in layer.keys():
                my_info[ele] = layer[ele].density
            info.append(my_info)

        # new roughness values
        parameters = [[11, 'C2'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [9, 'all']]
        density = [0.045, 0.08656, 0.0014, 0.096, 0.056, 0.045, 0.05, 0.025]

        for idx, params in enumerate(parameters):
            sample.setDensity(params[0], params[1], density[idx])

            # manually changes the parameters
            if params[1] == 'all':
                for key in info[params[0]].keys():
                    info[params[0]][key] = density[idx]
            else:
                info[params[0]][params[1]] = density[idx]

            # checks to make sure that all roughness parameters have been altered properly
            for i, layer in enumerate(info):
                for ele in layer.keys():
                    self.assertEqual(layer[ele], sample.structure[i][ele].density)

    def test_get_thickness(self):
        
        # Testing the getThickness function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[10, 'A'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [7, 'all'], [3,'all']]
        thickness = [3.97, 3.97, 10.862294545882191, 3.9972874500963402, 3.97, 3.97, 3.97, 3.97,3.94]

        for idx, params in enumerate(parameters):
            d = sample.getThickness(params[0], params[1])
            self.assertEqual(d, thickness[idx])

    def test_set_thickness(self):

        # Testing the setDensity function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        # retrieves initial configuration of the sample definition
        info = []
        for layer in sample.structure:
            my_info = {}
            for ele in layer.keys():
                my_info[ele] = layer[ele].thickness
            info.append(my_info)

        # new roughness values
        parameters = [[11, 'C2'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [9, 'all']]
        thickness = [25.5, 16.525, 0.55, 5, 3.9, 4.0, 6.5, 25]

        for idx, params in enumerate(parameters):
            sample.setThickness(params[0], params[1], thickness[idx])

            # manually changes the parameters
            if params[1] == 'all':
                for key in info[params[0]].keys():
                    info[params[0]][key] = thickness[idx]
            else:
                info[params[0]][params[1]] = thickness[idx]

            # checks to make sure that all roughness parameters have been altered properly
            for i, layer in enumerate(info):
                for ele in layer.keys():
                    self.assertEqual(layer[ele], sample.structure[i][ele].thickness)

    def test_get_totalThickness(self):

        # Testing the getCombinedThickness function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[0,1,'A'],[1,2, 'B'],[1,3, 'all'],[1,4,'A'],[1,5,'B'],[1,6, 'O'],[2,3, 'B'],[1,10, 'all'],
                      [5,7, 'A'],[7,10, 'B'],[0,11, 'O']]

        solutions = [53.905, 7.81, 11.75, 15.72, 19.69, 23.66, 7.845, 39.54, 11.91, 15.88, 93.53728745009633]

        for i,params in enumerate(parameters):
            dtot = sample.getTotalThickness(params[0],params[1],params[2])
            self.assertEqual(dtot,solutions[i])

    def test_set_combinedThickness(self):

        # Testing the setCombinedThickness function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[1, 2, 'B'], [1, 3, 'all'], [1, 4, 'A'], [1, 5, 'B'], [1, 6, 'O'], [2, 3, 'B'],
                      [1, 10, 'all'],
                      [5, 7, 'A'], [7, 10, 'B'], [0, 11, 'O']]

        total_thickness = [17, 4.5, 2.5, 5.55555, 19, 4, 30.54, 20.5, 17.5, 100.6]

        total_thickness_solutions = [17.0, 4.5, 2.5, 5.555549999999999, 18.999999999999996, 4.0, 30.54, 20.5, 17.5, 100.59999999999995]

        test_list = []
        for i, params in enumerate(parameters):
            sample.setCombinedThickness(params[0], params[1], params[2], total_thickness[i])
            test_list.append(sample.getTotalThickness(params[0], params[1], params[2]))

        self.assertListEqual(test_list, total_thickness_solutions)

    def test_setVariationConstant(self):

        # Tests setVariationConstant function for script feature
        sample = ds.ReadSampleHDF5(self.my_path3)

        # layers 1 to 10
        # Ti, Mn2, Mn3
        # self, layer, symbol, identifier, val

        parameters = [[1, 'Ti', 0.05], [2, 'Ti', 0.75], [3, 'Ti', 0.61], [4, 'Ti', 0.002], [5, 'Ti', 0.085],
                      [6, 'Ti', 0.99999], [7, 'Ti', 0.00001], [8, 'Ti', 0.5], [9, 'Ti', 0.25], [10, 'Ti', 0.84],
                      [1, 'Mn2', 0.05], [2, 'Mn2', 0.75], [3, 'Mn2', 0.61], [4, 'Mn2', 0.002], [5, 'Mn2', 0.085],
                      [6, 'Mn2', 0.99999], [7, 'Mn2', 0.00001], [8, 'Mn2', 0.5], [9, 'Mn2', 0.25], [10, 'Mn2', 0.84],
                      [1, 'Mn3', 0.05], [2, 'Mn3', 0.75], [3, 'Mn3', 0.61], [4, 'Mn3', 0.002], [5, 'Mn3', 0.085],
                      [6, 'Mn3', 0.99999], [7, 'Mn3', 0.00001], [8, 'Mn3', 0.5], [9, 'Mn3', 0.25], [10, 'Mn3', 0.84]]

        solutions = [{'Ti': 0.05, 'Mn2': 0.475, 'Mn3': 0.475},
                     {'Ti': 0.75, 'Mn2': 0.24792774892792938, 'Mn3': 0.0020722510720706192},
                     {'Ti': 0.61, 'Mn2': 0.1842952550116507, 'Mn3': 0.20570474498834931},
                     {'Ti': 0.002, 'Mn2': 0.27131489943446113, 'Mn3': 0.7266851005655389},
                     {'Ti': 0.085, 'Mn2': 0.0538550233401543, 'Mn3': 0.8611449766598457},
                     {'Ti': 0.99999, 'Mn2': 4.7245953881851626e-07, 'Mn3': 9.527540461135974e-06},
                     {'Ti': 1e-05, 'Mn2': 0.0, 'Mn3': 0.99999}, {'Ti': 0.5, 'Mn2': 0.0, 'Mn3': 0.5},
                     {'Ti': 0.25, 'Mn2': 1.21406976988786e-05, 'Mn3': 0.7499878593023012},
                     {'Ti': 0.84, 'Mn2': 0.06588935386026525, 'Mn3': 0.09411064613973479},
                     {'Ti': 0.09047619047619047, 'Mn2': 0.05, 'Mn3': 0.8595238095238095},
                     {'Ti': 0.24931115292808748, 'Mn2': 0.75, 'Mn3': 0.0006888470719125218},
                     {'Ti': 0.2916496458573352, 'Mn2': 0.61, 'Mn3': 0.0983503541426648},
                     {'Ti': 0.002739180475147477, 'Mn2': 0.002, 'Mn3': 0.9952608195248525},
                     {'Ti': 0.08220199009518324, 'Mn2': 0.085, 'Mn3': 0.8327980099048169},
                     {'Ti': 9.999904724504866e-06, 'Mn2': 0.99999, 'Mn3': 9.527544962472112e-11},
                     {'Ti': 9.9999e-06, 'Mn2': 1e-05, 'Mn3': 0.9999800001000001},
                     {'Ti': 0.25, 'Mn2': 0.5, 'Mn3': 0.25},
                     {'Ti': 0.1875022764084557, 'Mn2': 0.25, 'Mn3': 0.5624977235915443},
                     {'Ti': 0.14388017153579788, 'Mn2': 0.84, 'Mn3': 0.016119828464202162},
                     {'Ti': 0.611864406779661, 'Mn2': 0.338135593220339, 'Mn3': 0.05},
                     {'Ti': 0.06237075214200787, 'Mn2': 0.18762924785799212, 'Mn3': 0.75},
                     {'Ti': 0.12615028731720695, 'Mn2': 0.26384971268279306, 'Mn3': 0.61},
                     {'Ti': 0.5768301351959196, 'Mn2': 0.4211698648040803, 'Mn3': 0.002},
                     {'Ti': 0.44984405325723126, 'Mn2': 0.4651559467427688, 'Mn3': 0.085},
                     {'Ti': 9.999904725412102e-11, 'Mn2': 9.999900000907236e-06, 'Mn3': 0.99999},
                     {'Ti': 0.49999250001250006, 'Mn2': 0.4999974999874999, 'Mn3': 1e-05},
                     {'Ti': 0.16666666666666666, 'Mn2': 0.3333333333333333, 'Mn3': 0.5},
                     {'Ti': 0.3214308013681088, 'Mn2': 0.42856919863189114, 'Mn3': 0.25},
                     {'Ti': 0.023397999178897025, 'Mn2': 0.136602000821103, 'Mn3': 0.84}]

        for i,params in enumerate(parameters):
            sample.setVariationConstant(params[0],'B',params[1],params[2])
            for j, key in enumerate(sample.structure[params[0]]['B'].polymorph):
                self.assertEqual(sample.structure[params[0]]['B'].poly_ratio[j], solutions[i][key])

    def test_setMultiVarConstant(self):

        # Tests setMultiVarConstant function for script feature
        sample = ds.ReadSampleHDF5(self.my_path4)

        parameters = [[2, 'B', ['Ti', 'Mn2+'], [0.95, 0]], [3, 'B', ['Ti', 'Mn2+'], [0.44, 0.25]],
                  [4, 'A', ['Sr'], [0.777]], [5, 'B', ['Ti', 'Mn3+'], [0.22, 0.22]],
                  [6, 'B', ['Ti', 'Mn4+'], [0.55, 0.1]], [7, 'A', ['La'], [0.22]],
                  [8, 'B', ['Mn3+', 'Mn2+'], [0.1, 0.01]], [9, 'B', ['Mn3+', 'Mn4+'], [0.01, 0.7]],
                  [12, 'B', ['Mn4+', 'Mn2+'], [0.548, 0.22]],[10, 'A', ['Sr'], [0.0001]]]

        solutions = [{'Ti': 0.95, 'Mn2+': 0.0, 'Mn3+': 8.852355544947185e-29, 'Mn4+': 0.050000000000000044},
                     {'Ti': 0.44, 'Mn2+': 0.25, 'Mn3+': 0.0567722767814002, 'Mn4+': 0.2532277232185999},
                     {'Sr': 0.777, 'La': 0.22299999999999998},
                     {'Ti': 0.22, 'Mn2+': 0.4410992818927896, 'Mn3+': 0.22, 'Mn4+': 0.11890071810721055},
                     {'Ti': 0.55, 'Mn2+': 0.0135181896388328, 'Mn3+': 0.3364818103611672, 'Mn4+': 0.1},
                     {'Sr': 0.78, 'La': 0.22}, {'Ti': 0.0, 'Mn2+': 0.1, 'Mn3+': 0.01, 'Mn4+': 0.8900000000000001},
                     {'Ti': 0.14500000000000002, 'Mn2+': 0.14500000000000002, 'Mn3+': 0.01, 'Mn4+': 0.7},
                     {'Ti': 0.0, 'Mn2+': 0.548, 'Mn3+': 0.232, 'Mn4+': 0.22}, {'Sr': 0.0001, 'La': 0.9999}]

        for i,params in enumerate(parameters):
            sample.setMultiVarConstant(params[0],params[1],params[2],params[3])
            for j, key in enumerate(sample.structure[params[0]][params[1]].polymorph):
                self.assertEqual(sample.structure[params[0]][params[1]].poly_ratio[j], solutions[i][key])

    def test_setRatio(self):

        # Tests setRatio function for script feature
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[0, 'B', 'Ti', 'Mn2', 0.5], [1, 'B', 'Ti', 'Mn3', 0.1], [2, 'B', 'Mn3', 'Mn2', 0.9],
                      [3, 'B', 'Mn2', 'Mn3', 0.75], [4, 'B', 'Ti', 'Mn2', 0.55], [5, 'B', 'Mn2', 'Ti', 0.925],
                      [6, 'B', 'Mn3', 'Ti', 0.4], [7, 'B', 'Ti', 'Mn2', 0.33333337], [8, 'B', 'Ti', 'Mn3', 0.84],
                      [9, 'B', 'Mn3', 'Mn2', 0.65315315131], [00, 'B', 'Mn2', 'Mn3', 0.22232]]

        solutions = [{'Ti': 0.6666666666666666, 'Mn2': 0.3333333333333333, 'Mn3': 0.0},
                     {'Ti': 0.9090909090909091, 'Mn2': 0.0, 'Mn3': 0.09090909090909091},
                     {'Ti': 0.8764781545769192, 'Mn2': 0.058510347831985624, 'Mn3': 0.06501149759109513},
                     {'Ti': 0.7221884925576368, 'Mn2': 0.15874943282420753, 'Mn3': 0.11906207461815564},
                     {'Ti': 0.32497901795260065, 'Mn2': 0.17873845987393036, 'Mn3': 0.496282522173469},
                     {'Ti': 0.0357005070581337, 'Mn2': 0.03859514276554994, 'Mn3': 0.9257043501763164},
                     {'Ti': 0.27225198485360036, 'Mn2': 0.0471180530123988, 'Mn3': 0.680629962134001},
                     {'Ti': 0.0, 'Mn2': 0.0, 'Mn3': 1.0},
                     {'Ti': 0.5434782608695653, 'Mn2': 0.0, 'Mn3': 0.4565217391304348},
                     {'Ti': 0.0, 'Mn2': 0.395095367172984, 'Mn3': 0.604904632827016},
                     {'Ti': 0.6666666666666666, 'Mn2': 0.2727054562907695, 'Mn3': 0.060627877042563866}]

        for i,params in enumerate(parameters):
            sample.setRatio(params[0],params[1],params[2],params[3], params[4])
            for j, key in enumerate(sample.structure[params[0]][params[1]].polymorph):
                self.assertEqual(sample.structure[params[0]][params[1]].poly_ratio[j], solutions[i][key])

    def test_get_magDensity(self):

        # Testing the getMagDensity function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        parameters = [[2, 'B','Mn3'], [3, 'B','Mn3'], [4, 'B','Mn3'], [5, 'B','Mn3'], [6, 'B','Mn3'],
                      [7, 'B', 'Mn3'], [8, 'B', 'Mn3'], [9, 'B', 'Mn3'], [10, 'B', 'Mn3']]
        mag_density = [0.0003066750126945097, 0.0002068731239481137, 0.016932058948664633, 0.02655782307499808,
                     0.0314329491242263, 0.0314329491242263, 0.0314329491242263, 0.0314329491242263,
                     0.00010936761100236586]

        for idx, params in enumerate(parameters):
            rho = sample.getMagDensity(params[0], params[1], params[2])
            self.assertEqual(rho, mag_density[idx])

    def test_set_magDensity(self):

        # Testing the setDensity function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        # retrieves initial configuration of the sample definition

        # new roughness values
        parameters = [[2, 'B', 'Mn3'], [3, 'B', 'Mn3'], [4, 'B', 'Mn3'], [5, 'B', 'Mn3'], [6, 'B', 'Mn3'],
                      [7, 'B', 'Mn3'], [8, 'B', 'Mn3'], [9, 'B', 'Mn3'], [10, 'B', 'Mn3']]

        mag_eval = [0.0003066750126945097, 0.0002068731239481137, 0.016932058948664633, 0.02655782307499808,
                       0.0314329491242263, 0.0314329491242263, 0.0314329491242263, 0.0314329491242263,
                       0.00010936761100236586]

        mag_new = [0.01,0.25,0.0002,0.3,5,0.000095, 0.08, 0.004, 0.00025]

        for idx, params in enumerate(parameters):
            mag = mag_new[idx]
            sample.setMagDensity(params[0], params[1], params[2], mag)

            mag_eval[idx] = mag
            # checks to make sure that all roughness parameters have been altered properly

            # retrieves all information for changed parameters
            test_list = []
            for i in range(2, 11):
                my_idx = list(sample.structure[i]['B'].polymorph).index('Mn3')
                mag_density = sample.structure[i]['B'].mag_density[my_idx]
                test_list.append(mag_density)

            self.assertListEqual(test_list, mag_eval)

    def test_get_eShift(self):

        # Testing the getEshift function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        eShift = {'O': 0, 'C': 0, 'Sr': 0, 'La': -0.1, 'A': 0, 'B': 0, 'Ti': 0, 'Mn2': -0.95, 'Mn3': -0.95}

        for key in eShift.keys():
            eShift_func = sample.getEshift(key)
            self.assertEqual(eShift_func, eShift[key])

    def test_set_eShift(self):

        # Testing the getEshift function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        eShift_new = {'O': 1.5, 'C': -1, 'Sr': -0.0001, 'La': -1.1, 'A': 0.652, 'B': 2.11, 'Ti': 0.05, 'Mn2': 0, 'Mn3': 0.55}
        eShift = {'O': 0, 'C': 0, 'Sr': 0, 'La': -0.1, 'A': 0, 'B': 0, 'Ti': 0, 'Mn2': -0.95, 'Mn3': -0.95}

        for key in eShift_new.keys():
            sample.setEshift(key, eShift_new[key])
            eShift[key] = eShift_new[key]
            for test_key in eShift_new.keys():
                value = sample.eShift[test_key]
                self.assertEqual(value, eShift[test_key])

    def test_get_eShift_mag(self):

        # Testing the getMagEshift function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        mag_eShift = {'Ni': -0.95, 'Co': -0.95}

        for key in mag_eShift.keys():
            eShift_func = sample.getMagEshift(key)
            self.assertEqual(eShift_func, mag_eShift[key])

    def test_set_eShift_mag(self):

        # Testing the getMagEshift function used in the script
        sample = ds.ReadSampleHDF5(self.my_path3)

        mag_eShift = {'Ni': -0.95, 'Co': -0.95}
        mag_eShift_new = {'Ni': -0.55, 'Co': 1.2}

        for key in mag_eShift_new.keys():
            sample.setMagEshift(key, mag_eShift_new[key])
            mag_eShift[key] = mag_eShift_new[key]
            for test_key in mag_eShift_new.keys():
                value = sample.mag_eShift[test_key]
                self.assertEqual(value, mag_eShift[test_key])


if __name__ == "__main__":
    unittest.main()