import os
import copy
import unittest
import sys 

# Get the parent directory of the current script's directory
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
# Add the parent directory to the system path
sys.path.append(parent_dir)

from UTILS import (
    material_structure as ms,
    data_structure as ds,
    global_optimization as go,
    ROOT_DIR, TESTS_DIR,
)
from UTILS.global_optimization import changeSampleParams
from GUI_GO import checkscript

# This test script can be executed by inputting
#  ->  python -m unittest -v test_data_fitting.py
# into the terminal

class TestDataFitting(unittest.TestCase):

    def __init__(self, methodName: str = "runTest") -> None:
        super().__init__(methodName)

        self.tests_dir = TESTS_DIR
        self.filename = 'Pim4uc_test.h5'
        self.my_path = self.tests_dir + '/test_data/' + self.filename
        self.script_path = self.tests_dir + '/test_data/test_script.txt'

    def test_ChangeSampleParams_element(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [[1, 'STRUCTURAL', 'ELEMENT', 'A', 'THICKNESS'], [3, 'STRUCTURAL', 'ELEMENT', 'A', 'THICKNESS'],
                      [2, 'STRUCTURAL', 'ELEMENT', 'Ti', 'DENSITY'], [4, 'STRUCTURAL', 'ELEMENT', 'O', 'DENSITY'],
                      [5, 'STRUCTURAL', 'ELEMENT', 'O', 'ROUGHNESS'], [6, 'STRUCTURAL', 'ELEMENT', 'Mn', 'ROUGHNESS'],
                      [7, 'STRUCTURAL', 'ELEMENT', 'C2', 'LINKED ROUGHNESS']]

        x = [10, 20, 0.025, 0.086, 0.5, 5, 2.75]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        # my_script, problem, my_error = checkscript(sample, fname=script_name)  # load in the script
        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF, my_script, orbitals)

        # test sample parameters
        sample_solution = copy.deepcopy(sample)
        sample_solution.structure[1]['A'].thickness = 10
        sample_solution.structure[3]['A'].thickness = 20
        sample_solution.structure[2]['Ti'].density = 0.025
        sample_solution.structure[4]['O'].density = 0.086
        sample_solution.structure[5]['O'].roughness = 0.5
        sample_solution.structure[6]['Mn'].roughness = 5
        sample_solution.structure[7]['C2'].linked_roughness = 2.75

        for i in range(0,8):
            for key in list(sample_solution.structure[i].keys()):
                test_thickness = sample_new.structure[i][key].thickness
                sol_thickness = sample_solution.structure[i][key].thickness
                test_density = sample_new.structure[i][key].density
                sol_density = sample_solution.structure[i][key].density
                test_roughness = sample_new.structure[i][key].roughness
                sol_roughness = sample_solution.structure[i][key].roughness
                test_linked = sample_new.structure[i][key].linked_roughness
                sol_linked = sample_solution.structure[i][key].linked_roughness

                self.assertEqual(test_thickness, sol_thickness)
                self.assertEqual(test_density, sol_density)
                self.assertEqual(test_roughness, sol_roughness)
                self.assertEqual(test_linked, sol_linked)

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

    def test_ChangeSampleParams_compound(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [[0, 'STRUCTURAL', 'COMPOUND','ROUGHNESS',0],[1, 'STRUCTURAL', 'COMPOUND','THICKNESS',1],
                      [2, 'STRUCTURAL', 'COMPOUND','DENSITY',2], [3, 'STRUCTURAL', 'COMPOUND','ROUGHNESS',0],
                      [4, 'STRUCTURAL', 'COMPOUND','THICKNESS',1], [5, 'STRUCTURAL', 'COMPOUND','DENSITY',1],
                      [6, 'STRUCTURAL', 'COMPOUND','ROUGHNESS',0], [7, 'STRUCTURAL', 'COMPOUND','THICKNESS',1]]



        x = [1.5, 1.95, 0.074, 4, 3.5, 0.025, 2.15, 12]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        # test sample parameters
        sample_solution = copy.deepcopy(sample)
        sample_solution.structure[0]['A'].roughness = 1.5
        sample_solution.structure[0]['Ti'].roughness = 3.295897580593646
        sample_solution.structure[0]['O'].roughness = 3.295897580593646
        sample_solution.structure[1]['A'].thickness = 1.95
        sample_solution.structure[1]['Ti'].thickness = 1.95
        sample_solution.structure[1]['O'].thickness = 1.95
        sample_solution.structure[2]['A'].density = 0.01799999999999999
        sample_solution.structure[2]['Ti'].density = 0.01799999999999999
        sample_solution.structure[2]['O'].density = 0.074
        sample_solution.structure[3]['A'].roughness = 4
        sample_solution.structure[3]['Mn'].roughness = 4.0
        sample_solution.structure[3]['O'].roughness = 4.75
        sample_solution.structure[4]['A'].thickness = 3.5
        sample_solution.structure[4]['Mn'].thickness = 3.5
        sample_solution.structure[4]['O'].thickness = 3.5
        sample_solution.structure[5]['A'].density = 0.025
        sample_solution.structure[5]['Mn'].density = 0.025
        sample_solution.structure[5]['O'].density = 0.07500000000000001
        sample_solution.structure[6]['A'].roughness = 2.15
        sample_solution.structure[6]['Mn'].roughness = 1.9
        sample_solution.structure[6]['O'].roughness = 1.9
        sample_solution.structure[7]['C1'].thickness = 7.364428710924312
        sample_solution.structure[7]['C2'].thickness = 12
        sample_solution.structure[7]['O'].thickness = 5.203439593060714

        for i in range(0, 8):
            for key in list(sample_solution.structure[i].keys()):
                test_thickness = sample_new.structure[i][key].thickness
                sol_thickness = sample_solution.structure[i][key].thickness
                test_density = sample_new.structure[i][key].density
                sol_density = sample_solution.structure[i][key].density
                test_roughness = sample_new.structure[i][key].roughness
                sol_roughness = sample_solution.structure[i][key].roughness
                test_linked = sample_new.structure[i][key].linked_roughness
                sol_linked = sample_solution.structure[i][key].linked_roughness

                self.assertEqual(test_thickness, sol_thickness)
                self.assertEqual(test_density, sol_density)
                self.assertEqual(test_roughness, sol_roughness)
                self.assertEqual(test_linked, sol_linked)

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

    def test_ChangeSampleParams_variation(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [[0, 'POLYMORPHOUS', 'A', 'Sr'], [1, 'POLYMORPHOUS', 'A', 'La'],
                      [2, 'POLYMORPHOUS', 'A', 'Sr'], [3, 'POLYMORPHOUS', 'Mn', 'Mn2+'],
                      [4, 'POLYMORPHOUS', 'Mn', 'Mn3+'], [5, 'POLYMORPHOUS', 'Mn', 'Mn2+'],
                      [6, 'POLYMORPHOUS', 'Mn', 'Mn3+']]

        x = [0.5, 0.15468, 0.7, 0.641, 0.9999, 0.0001, 0.123456789]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF, my_script, orbitals)

        # test sample parameters
        sample_solution = copy.deepcopy(sample)

        sample_solution = copy.deepcopy(sample)

        idx = list(sample_solution.structure[0]['A'].polymorph).index('Sr')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[0]['A'].poly_ratio[idx] = 0.5
        sample_solution.structure[0]['A'].poly_ratio[idx_no] = 1 - 0.5

        idx = list(sample_solution.structure[1]['A'].polymorph).index('La')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[1]['A'].poly_ratio[idx] = 0.15468
        sample_solution.structure[1]['A'].poly_ratio[idx_no] = 1 - 0.15468

        idx = list(sample_solution.structure[2]['A'].polymorph).index('Sr')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[2]['A'].poly_ratio[idx] = 0.7
        sample_solution.structure[2]['A'].poly_ratio[idx_no] = 1 - 0.7

        idx = list(sample_solution.structure[3]['Mn'].polymorph).index('Mn2+')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[3]['Mn'].poly_ratio[idx] = 0.641
        sample_solution.structure[3]['Mn'].poly_ratio[idx_no] = 1 - 0.641

        idx = list(sample_solution.structure[4]['Mn'].polymorph).index('Mn3+')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[4]['Mn'].poly_ratio[idx] = 0.9999
        sample_solution.structure[4]['Mn'].poly_ratio[idx_no] = 1 - 0.9999

        idx = list(sample_solution.structure[5]['Mn'].polymorph).index('Mn2+')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[5]['Mn'].poly_ratio[idx] = 0.0001
        sample_solution.structure[5]['Mn'].poly_ratio[idx_no] = 1 - 0.0001

        idx = list(sample_solution.structure[6]['Mn'].polymorph).index('Mn3+')
        idx_no = 0 if idx == 1 else 1
        sample_solution.structure[6]['Mn'].poly_ratio[idx] = 0.123456789
        sample_solution.structure[6]['Mn'].poly_ratio[idx_no] = 1 - 0.123456789

        for j, param in enumerate(parameters):
            idx_s = list(sample_solution.structure[param[0]][param[2]].polymorph).index(param[3])
            idx_s_no = 0 if idx_s == 1 else 1
            solution1 = sample_solution.structure[param[0]][param[2]].poly_ratio[idx_s]
            solution2 = sample_solution.structure[param[0]][param[2]].poly_ratio[idx_s_no]

            idx_t = list(sample_new.structure[param[0]][param[2]].polymorph).index(param[3])
            idx_t_no = 0 if idx_s == 1 else 1
            test1 = sample_new.structure[param[0]][param[2]].poly_ratio[idx_t]
            test2 = sample_new.structure[param[0]][param[2]].poly_ratio[idx_t_no]

            self.assertEqual(solution1, test1)
            self.assertEqual(solution2, test2)


        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

    def test_ChangeSampleParams_magnetic(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [[3, 'MAGNETIC', 'Mn', 'Mn2+'], [4, 'MAGNETIC', 'Mn', 'Mn3+'],
                      [5, 'MAGNETIC', 'Mn', 'Mn2+'], [6, 'MAGNETIC', 'Mn', 'Mn3+']]

        x = [0.001, 0.015,0.1648844464646, 0.00001313131313131563]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        sample_solution = copy.deepcopy(sample)

        idx = list(sample_solution.structure[3]['Mn'].polymorph).index('Mn2+')
        sample_solution.structure[3]['Mn'].mag_density[idx] = 0.001

        idx = list(sample_solution.structure[4]['Mn'].polymorph).index('Mn3+')
        sample_solution.structure[4]['Mn'].mag_density[idx] = 0.015

        idx = list(sample_solution.structure[5]['Mn'].polymorph).index('Mn2+')
        sample_solution.structure[5]['Mn'].mag_density[idx] = 0.1648844464646

        idx = list(sample_solution.structure[6]['Mn'].polymorph).index('Mn2+')
        sample_solution.structure[6]['Mn'].mag_density[idx] = 0.00001313131313131563

        for i, param in enumerate(parameters):
            idx_s = list(sample_solution.structure[param[0]][param[2]].polymorph).index(param[3])
            solution = sample_solution.structure[param[0]][param[2]].mag_density[idx_s]

            idx_t = list(sample_new.structure[param[0]][param[2]].polymorph).index(param[3])
            test_value = sample_new.structure[param[0]][param[2]].mag_density[idx_t]

            self.assertEqual(solution, test_value)


        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

    def test_ChangeSampleParams_eShift(self):
        
        # test scattering factors
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [['SCATTERING FACTOR', 'STRUCTURAL', 'Sr'], ['SCATTERING FACTOR', 'STRUCTURAL', 'Ti'],
                      ['SCATTERING FACTOR', 'STRUCTURAL', 'O'], ['SCATTERING FACTOR', 'STRUCTURAL', 'La'],
                      ['SCATTERING FACTOR', 'STRUCTURAL', 'Mn2'], ['SCATTERING FACTOR', 'STRUCTURAL', 'Mn3']]

        x = [-0.55,0.25,1.5,-2.5,1,-2]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)


        solutions = {'Sr': -0.55, 'Ti': 0.25, 'O':1.5, 'La':-2.5, 'Mn2':1,'Mn3':-2}

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

        for key in solutions.keys():
            self.assertEqual(sample_new.eShift[key], solutions[key])

    def test_ChangeSampleParams_eShift_mag(self):

        # test magnetic scattering factors
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [['SCATTERING FACTOR', 'MAGNETIC', 'Co'], ['SCATTERING FACTOR', 'MAGNETIC', 'Ni']]

        x = [0.1254545,-0.89797979]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        solutions = {'Co': 0.1254545, 'Ni': -0.89797979}

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

        for key in solutions.keys():
            self.assertEqual(sample_new.mag_eShift[key], solutions[key])

    def test_ChangeSampleParams_orbitals(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [['ORBITAL', 'DEXY', 'Mn2'], ['ORBITAL', 'DEXZYZ', 'Mn2'], ['ORBITAL', 'DEX2Y2', 'Mn2'], ['ORBITAL', 'DEZZ', 'Mn2']]

        x = [0.5,-0.25,1,-0.67]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        sample_solution = copy.deepcopy(sample)

        for i in range(0, 8):
            for key in list(sample_solution.structure[i].keys()):
                test_thickness = sample_new.structure[i][key].thickness
                sol_thickness = sample_solution.structure[i][key].thickness
                test_density = sample_new.structure[i][key].density
                sol_density = sample_solution.structure[i][key].density
                test_roughness = sample_new.structure[i][key].roughness
                sol_roughness = sample_solution.structure[i][key].roughness
                test_linked = sample_new.structure[i][key].linked_roughness
                sol_linked = sample_solution.structure[i][key].linked_roughness

                self.assertEqual(test_thickness, sol_thickness)
                self.assertEqual(test_density, sol_density)
                self.assertEqual(test_roughness, sol_roughness)
                self.assertEqual(test_linked, sol_linked)

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        orbitals_solution = {'Mn2': [0.5,-0.25,1,-0.67]}
        # test orbital
        self.assertListEqual(orbitals_solution['Mn2'], orbitals_new['Mn2'])

    def test_ChangeSampleParams_scaling(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [['SCALING FACTOR', '107_799.82_S'], ['SCALING FACTOR', '51-52_640.2_AL_Asymm'],
                      ['SCALING FACTOR', '110_E600.18_Th25.0_LC'], ['SCALING FACTOR', '114-115_E600.18_Th35.0_AC_Asymm'],
                      ['SCALING FACTOR', '17_E600.18_Th25.0_S']]

        x = [0.8, 1.2, 1, 1.1, 0.75]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        sample_solution = copy.deepcopy(sample)

        for i in range(0, 8):
            for key in list(sample_solution.structure[i].keys()):
                test_thickness = sample_new.structure[i][key].thickness
                sol_thickness = sample_solution.structure[i][key].thickness
                test_density = sample_new.structure[i][key].density
                sol_density = sample_solution.structure[i][key].density
                test_roughness = sample_new.structure[i][key].roughness
                sol_roughness = sample_solution.structure[i][key].roughness
                test_linked = sample_new.structure[i][key].linked_roughness
                sol_linked = sample_solution.structure[i][key].linked_roughness

                self.assertEqual(test_thickness, sol_thickness)
                self.assertEqual(test_density, sol_density)
                self.assertEqual(test_roughness, sol_roughness)
                self.assertEqual(test_linked, sol_linked)

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

        scaleF_solutions = copy.deepcopy(scaleF)
        scaleF_solutions['107_799.82_S'] = '0.8'
        scaleF_solutions['51-52_640.2_AL_Asymm'] = '1.2'
        scaleF_solutions['110_E600.18_Th25.0_LC'] = '1'
        scaleF_solutions['114-115_E600.18_Th35.0_AC_Asymm'] = '1.1'
        scaleF_solutions['17_E600.18_Th25.0_S'] = '0.75'

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF_solutions[skey], scaleF_new[skey])

    def test_ChangeSampleParams_all(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        print('self.my_path', self.my_path)
        # aux = input('Enter to continue')

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [['SCALING FACTOR', 'ALL SCANS'], ['BACKGROUND SHIFT', 'ALL SCANS']]

        x = [1.5, 1e-7]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        sample_solution = copy.deepcopy(sample)

        for i in range(0, 8):
            for key in list(sample_solution.structure[i].keys()):
                test_thickness = sample_new.structure[i][key].thickness
                sol_thickness = sample_solution.structure[i][key].thickness
                test_density = sample_new.structure[i][key].density
                sol_density = sample_solution.structure[i][key].density
                test_roughness = sample_new.structure[i][key].roughness
                sol_roughness = sample_solution.structure[i][key].roughness
                test_linked = sample_new.structure[i][key].linked_roughness
                sol_linked = sample_solution.structure[i][key].linked_roughness

                self.assertEqual(test_thickness, sol_thickness)
                self.assertEqual(test_density, sol_density)
                self.assertEqual(test_roughness, sol_roughness)
                self.assertEqual(test_linked, sol_linked)

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

        backS_solutions = copy.deepcopy(backS)
        for name in backS.keys():
            backS_solutions[name] = '1.000000e-07'

        scaleF_solutions = copy.deepcopy(scaleF)
        for name in scaleF.keys():
            scaleF_solutions[name] = '1.5'

        # test scaling factor
        for bkey in list(backS.keys()):
            self.assertEqual(backS_solutions[bkey], backS_new[bkey])
            # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF_solutions[skey], scaleF_new[skey])

    def test_ChangeSampleParams_shift(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)

        parameters = [['BACKGROUND SHIFT', '107_799.82_S'], ['BACKGROUND SHIFT', '51-52_640.2_AL_Asymm'],
                      ['BACKGROUND SHIFT', '110_E600.18_Th25.0_LC'],
                      ['BACKGROUND SHIFT', '114-115_E600.18_Th35.0_AC_Asymm'],
                      ['BACKGROUND SHIFT', '17_E600.18_Th25.0_S']]

        x = [5e-8, -5e-8, 1, 1e-9, -0.5e-9]  # parameter values

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals)

        sample_solution = copy.deepcopy(sample)

        for i in range(0, 8):
            for key in list(sample_solution.structure[i].keys()):
                test_thickness = sample_new.structure[i][key].thickness
                sol_thickness = sample_solution.structure[i][key].thickness
                test_density = sample_new.structure[i][key].density
                sol_density = sample_solution.structure[i][key].density
                test_roughness = sample_new.structure[i][key].roughness
                sol_roughness = sample_solution.structure[i][key].roughness
                test_linked = sample_new.structure[i][key].linked_roughness
                sol_linked = sample_solution.structure[i][key].linked_roughness

                self.assertEqual(test_thickness, sol_thickness)
                self.assertEqual(test_density, sol_density)
                self.assertEqual(test_roughness, sol_roughness)
                self.assertEqual(test_linked, sol_linked)

        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])

        backS_solutions = copy.deepcopy(backS)
        backS_solutions['107_799.82_S'] = '5.000000e-08'
        backS_solutions['51-52_640.2_AL_Asymm'] = '-5.000000e-08'
        backS_solutions['110_E600.18_Th25.0_LC'] = '1.000000e+00'
        backS_solutions['114-115_E600.18_Th35.0_AC_Asymm'] = '1.000000e-09'
        backS_solutions['17_E600.18_Th25.0_S'] = '-5.000000e-10'

        # test scaling factor
        for skey in list(backS.keys()):
            self.assertEqual(backS_solutions[skey], backS[skey])

    def test_ChangeSampleParams_script(self):

        # Tests to element fit
        sample = ds.ReadSampleHDF5(self.my_path)

        data, data_dict, sim_dict = ds.ReadDataHDF5(self.my_path)
        parameters = []
        x = []

        # create the background shift and scaling factor keys
        backS = dict()
        scaleF = dict()
        for name in list(data_dict.keys()):
            backS[name] = 0
            scaleF[name] = 1

        my_script, problem, my_error = checkscript(sample, fname=self.script_path, testing=True)  # load in the script
        orbitals = {'Mn2': [0, 0, 0, 0]}  # orbitals

        sample_new, backS_new, scaleF_new, orbitals_new = changeSampleParams(x, parameters, sample, backS, scaleF,
                                                                             my_script, orbitals, use_script=True)


        self.assertEqual(sample_new.structure[6]['Mn'].density, 0.028)
        self.assertEqual(sample.eShift['Mn2'], 0.1)
        self.assertEqual(sample.eShift['Mn3'], -1.1)
        self.assertEqual(sample.mag_eShift['Co'], 0.25)
        self.assertEqual(sample.mag_eShift['Ni'], -1.1)


        # test background shift
        for bkey in list(backS.keys()):
            self.assertEqual(backS[bkey], backS_new[bkey])

        # test scaling factor
        for skey in list(scaleF.keys()):
            self.assertEqual(scaleF[skey], scaleF_new[skey])

        # test orbital
        self.assertListEqual(orbitals['Mn2'], orbitals_new['Mn2'])


if __name__ == '__main__':
    unittest.main()