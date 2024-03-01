import os
import sys
# Add the parent directory of GUI_GO.py to the Python path
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import UTILS.material_model as mm
import numpy as np
import UTILS.material_structure as ms
import unittest

# This test script can be executed by inputting
#  ->  python -m unittest -v test_material_model.py
# into the terminal

class TestMaterialModel(unittest.TestCase):
    
    def test_form_factors_E(self):
        form_factors = ['La','Mn','O']
        energy = [0,200.56,300.48,456.25,788,1200]

        solution = [np.array([0, 0]), np.array([27.39474991,  2.19130417]), np.array([21.80144234,  7.77745056]), np.array([24.62693753,  9.74610276]), np.array([8.22143442, 6.59496319]), np.array([39.10212769, 27.12182917]), np.array([0, 0]), np.array([13.26296037,  5.51149566]), np.array([14.9012543 ,  4.34598672]), np.array([14.14087443,  2.78105032]), np.array([17.69976568, 14.30607907]), np.array([24.11308175,  8.02955793]), np.array([0, 0]), np.array([6.60493548, 1.09012143]), np.array([6.3101332 , 0.56313067]), np.array([5.42543942, 0.29771144]), np.array([7.76802606, 2.55593633]), np.array([8.34408952, 1.29813867])]
        test_case = []
        for ff in form_factors:
            for E in energy:
                test_case.append(mm.find_form_factor(ff,E,False))

        my_sum = 0
        for idx in range(len(test_case)):
            my_sum += sum(abs(solution[idx]-test_case[idx]))
        self.assertTrue(my_sum<1e-7)
        #self.assertFalse(are_nested_lists_equal(solution, test_case))

    def test_form_factors_Elist(self):

        solution = [np.array([0., 0.]), np.array([27.39474991,  2.19130417]), np.array([21.80144234,  7.77745056]), np.array([24.62693753,  9.74610276]), np.array([8.22143442, 6.59496319]), np.array([39.10212769, 27.12182917]), np.array([0., 0.]), np.array([13.26296037,  5.51149566]), np.array([14.9012543 ,  4.34598672]), np.array([14.14087443,  2.78105032]), np.array([17.69976568, 14.30607907]), np.array([24.11308175,  8.02955793]), np.array([0., 0.]), np.array([6.60493548, 1.09012143]), np.array([6.3101332 , 0.56313067]), np.array([5.42543942, 0.29771144]), np.array([7.76802606, 2.55593633]), np.array([8.34408952, 1.29813867])]

        form_factors = ['La', 'Mn', 'O']
        energy = [0, 200.56, 300.48, 456.25, 788, 1200]

        test_case = []
        for ff in form_factors:
            recast = [val for val in mm.find_form_factor(ff, energy, False)]
            test_case += recast

        my_sum = 0
        for idx in range(len(test_case)):
            my_sum += sum(abs(solution[idx] - test_case[idx]))

        self.assertTrue(my_sum < 1e-7)

    def test_form_factors_E_mag(self):
        form_factors = ['Co','Ni']
        energy = [400,600,625,641.51,648.25,800]

        solution = [np.array([0., 0.]), np.array([0., 0.]), np.array([-0.00681891,  0.0]), np.array([0.30686376, 2.47434212]), np.array([ 0.24675064, -0.34005255]), np.array([0., 0.]), np.array([0., 0.]), np.array([0., 0.]), np.array([0.1187884 , 0.00154338]), np.array([-4.77536493, -3.46825118]), np.array([-0.92991243, -0.23784557]), np.array([0., 0.])]
        test_case = []
        for ff in form_factors:
            for E in energy:
                test_case.append(mm.find_form_factor(ff,E,True))

        my_sum = 0
        for idx in range(len(test_case)):
            my_sum += sum(abs(solution[idx]-test_case[idx]))

        self.assertTrue(my_sum <1e-7)
        #self.assertFalse(are_nested_lists_equal(solution, test_case))

    def test_form_factors_Elist_mag(self):

        solution = [np.array([0., 0.]), np.array([0., 0.]), np.array([-0.00681891,  0.        ]), np.array([0.30686376, 2.47434212]), np.array([ 0.24675064, -0.34005255]), np.array([0., 0.]), np.array([0., 0.]), np.array([0., 0.]), np.array([0.1187884 , 0.00154338]), np.array([-4.77536493, -3.46825118]), np.array([-0.92991243, -0.23784557]), np.array([0., 0.])]

        form_factors = ['Co','Ni']
        energy = [400,600,625,641.51,648.25,800]

        test_case = []
        for ff in form_factors:
            recast = [val for val in mm.find_form_factor(ff, energy, True)]
            test_case += recast

        my_sum = 0
        for idx in range(len(test_case)):
            my_sum += sum(abs(solution[idx] - test_case[idx]))

        self.assertTrue(my_sum <1e-7)

    def test_MOC(self):
        filename = 'optical_energy.txt'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        solution = np.loadtxt(my_path)

        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

        sample.magnetization(1, ['Mn2', 'Mn3'], [0.015, 0.01], ['Co', 'Ni'])

        sample.energy_shift()

        energy = np.array([641, 642.1, 650.5, 653])

        thickness, density, density_magnetic = sample.density_profile(step=0.1)  # Computes the density profile
        # Magnetic Scattering Factor
        sfm = dict()
        sf = dict()

        for e in sample.find_sf[0].keys():
            sf[e] = mm.find_form_factor(sample.find_sf[0][e], energy, False)
        # Magnetic Scattering Factor
        for em in sample.find_sf[1].keys():
            sfm[em] = mm.find_form_factor(sample.find_sf[1][em], energy, True)

        d_len = len(thickness)

        delta_m, beta_m = mm.MOC(density_magnetic, sfm, energy,
                                 d_len)  # absorptive and dispersive components for magnetic components

        total = 0
        total += sum(abs(delta_m[0, :] - solution[8, :])) + sum(abs(delta_m[1, :] - solution[9, :])) + sum(
            abs(delta_m[2, :] - solution[10, :])) + sum(abs(delta_m[3, :] - solution[11, :]))
        total += sum(abs(beta_m[0, :] - solution[12, :])) + sum(abs(beta_m[1, :] - solution[13, :])) + sum(
            abs(beta_m[2, :] - solution[14, :])) + sum(abs(beta_m[3, :] - solution[15, :]))

        self.assertTrue(total < 1e-7)

    def test_magnetic_optical_constant(self):

        filename = 'optical_theta.txt'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        solution = np.loadtxt(my_path)

        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

        sample.magnetization(1, ['Mn2', 'Mn3'], [0.015, 0.01], ['Co', 'Ni'])

        sample.energy_shift()

        E = 641

        thickness, density, density_magnetic = sample.density_profile()

        sfm = dict()  # scattering factors of magnetic components
        sfm_dict = {}


        # Magnetic Scattering Factor
        for em in sample.find_sf[1].keys():
            sfm[em] = mm.find_form_factor(sample.find_sf[1][em], E, True)

        delta_m, beta_m = mm.magnetic_optical_constant(density_magnetic, sfm,E)  # calculates depth-dependent magnetic components

        total_delta = sum(abs(solution[2] - delta_m))
        total_beta = sum(abs(solution[3] - beta_m))

        self.assertTrue(total_delta<1e-7)
        self.assertTrue(total_beta<1e-7)

    def test_IoR(self):
        filename = 'optical_energy.txt'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        solution = np.loadtxt(my_path)

        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

        sample.magnetization(1, ['Mn2', 'Mn3'], [0.015, 0.01], ['Co', 'Ni'])

        sample.energy_shift()

        energy = np.array([641, 642.1, 650.5, 653])

        thickness, density, density_magnetic = sample.density_profile(step=0.1)  # Computes the density profile
        # Magnetic Scattering Factor
        sfm = dict()
        sf = dict()

        for e in sample.find_sf[0].keys():
            sf[e] = mm.find_form_factor(sample.find_sf[0][e], energy, False)

        delta, beta = mm.IoR(density, sf, energy)  # gets absorptive and dispersive components of refractive index


        total = 0
        total += sum(abs(delta[0,:]-solution[0,:])) + sum(abs(delta[1,:]-solution[1,:])) + sum(abs(delta[2,:]-solution[2,:])) + sum(abs(delta[3,:]-solution[3,:]))
        total += sum(abs(beta[0,:]-solution[4,:])) + sum(abs(beta[1,:]-solution[5,:])) + sum(abs(beta[2,:]-solution[6,:])) + sum(abs(beta[3,:]-solution[7,:]))

        print(total)
        self.assertTrue(total < 1e-7)

    def test_index_of_refraction(self):
        filename = 'optical_theta.txt'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
        else:
            my_path = os.getcwd() + '/test_data/' + filename

        solution = np.loadtxt(my_path)

        sample = ms.slab(3)
        sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
        sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                        linked_roughness=[False, 0.25, False])
        sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                        linked_roughness=[4, 4.5, False])

        sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

        sample.magnetization(1, ['Mn2', 'Mn3'], [0.015, 0.01], ['Co', 'Ni'])

        sample.energy_shift()

        E = 641

        thickness, density, density_magnetic = sample.density_profile()

        sf = dict()  # scattering factors of non-magnetic components
        sfm = dict()  # scattering factors of magnetic components
        sf_dict = {}
        sfm_dict = {}

        # Non-Magnetic Scattering Factor - no need to access original
        for e in sample.find_sf[0].keys():
            sf[e] = mm.find_form_factor(sample.find_sf[0][e], E, False)  # find the scattering factor at energy E + dE


        delta, beta = mm.index_of_refraction(density, sf, E)  # calculates depth-dependent refractive index components

        total_delta = sum(abs(delta-solution[0]))
        total_beta = sum(abs(beta - solution[1]))

        self.assertTrue(total_delta < 1e-7)
        self.assertTrue(total_beta < 1e-7)

if __name__ == '__main__':
    unittest.main()