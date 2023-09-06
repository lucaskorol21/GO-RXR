import unittest
import os
import numpy as np
import material_structure as ms
import pickle
import data_structure as ds

class TestMaterialStructure(unittest.TestCase):
    def test_adaptive_layer_segmentation(self):
        # Test the adaptive layer segmentation implementation
        filename = 'optical_energy.txt'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
            path = os.getcwd() + '/test_data/'
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename
            path = os.getcwd() + '/Testing/test_data/'

        optical = np.loadtxt(my_path)

        with open(path+'ALS_test.pkl', 'rb') as file:
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

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

        parameters = [[10, 'A'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [11, 'all']]
        roughness = [0, 0.5, 3.726885630758733, 2.0653473900870867, 0.0, 0.0, 0.0, 2]

        for idx, params in enumerate(parameters):
            sigma = sample.getRoughness(params[0],params[1])
            self.assertEqual(sigma, roughness[idx])


    def test_set_roughness(self):
        # Testing the setRoughness function used in the script

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

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

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

        parameters = [[10, 'A'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [8, 'all']]
        roughness = [0.0020259688666223116, 0.084, 0.08597982010497454, 0.046124560151674036, 0.028, 0.028, 0.084, 0.028]

        for idx, params in enumerate(parameters):
            rho = sample.getDensity(params[0], params[1])
            self.assertEqual(rho, roughness[idx])

    def test_set_density(self):
        # Testing the setDensity function used in the script
        pass

    def test_get_thickness(self):
        # Testing the getThickness function used in the script

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

        parameters = [[10, 'A'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [7, 'all'], [3,'all']]
        roughness = [3.97, 3.97, 10.862294545882191, 3.9972874500963402, 3.97, 3.97, 3.97, 3.97,3.94]

        for idx, params in enumerate(parameters):
            d = sample.getThickness(params[0], params[1])
            self.assertEqual(d, roughness[idx])

    def test_set_thickness(self):
        pass

    def test_get_combinedThickness(self):
        # Testing the getCombinedThickness function used in the script

        pass

    def test_set_combinedThickness(self):
        pass

    def test_setVariationConstant(self):
        pass

    def test_setMultiVarConstant(self):
        pass

    def test_setRatio(self):
        pass

    def test_get_magDensity(self):
        # Testing the getMagDensity function used in the script

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

        parameters = [[10, 'C2'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [11, 'all']]
        roughness = [2, 0.5, 3.726885630758733, 2.0653473900870867, 0.0, 0.0, 0.0, 2]

        for idx, params in enumerate(parameters):
            sigma = sample.getRoughness(params[0], params[1])
            self.assertEqual(sigma, roughness[idx])

    def test_set_magDensity(self):
        pass

    def test_get_eShift(self):
        # Testing the getEshift function used in the script

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

        parameters = [[10, 'C2'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [11, 'all']]
        roughness = [2, 0.5, 3.726885630758733, 2.0653473900870867, 0.0, 0.0, 0.0, 2]

        for idx, params in enumerate(parameters):
            sigma = sample.getRoughness(params[0], params[1])
            self.assertEqual(sigma, roughness[idx])

    def test_set_eShift(self):
        pass

    def test_get_eShift_mag(self):
        # Testing the getMagEshift function used in the script

        filename = '7uc_sample_test.h5'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename

        sample = ds.ReadSampleHDF5(my_path)

        parameters = [[10, 'C2'], [10, 'O'], [11, 'C2'], [11, 'O'], [4, 'A'], [4, 'B'], [4, 'O'], [11, 'all']]
        roughness = [2, 0.5, 3.726885630758733, 2.0653473900870867, 0.0, 0.0, 0.0, 2]

        for idx, params in enumerate(parameters):
            sigma = sample.getRoughness(params[0], params[1])
            self.assertEqual(sigma, roughness[idx])

    def test_set_eShift_mag(self):
        pass

if __name__ == "__main__":
    unittest.main()