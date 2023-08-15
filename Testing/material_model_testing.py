import material_model as mm
import numpy as np
import unittest


def are_nested_lists_equal(list1, list2):
    if len(list1) != len(list2):
        return False

    for sublist1, sublist2 in zip(list1, list2):
        if sublist1 != sublist2:
            return False

    return True

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
            my_sum += sum(solution[idx]-test_case[idx])

        self.assertFalse(my_sum==0)
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
            my_sum += sum(solution[idx] - test_case[idx])

        self.assertFalse(my_sum == 0)

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
            my_sum += sum(solution[idx]-test_case[idx])

        self.assertFalse(my_sum==0)
        #self.assertFalse(are_nested_lists_equal(solution, test_case))

    def test_form_factors_Elist_mag(self):

        solution = solution = [np.array([0., 0.]), np.array([0., 0.]), np.array([-0.00681891,  0.        ]), np.array([0.30686376, 2.47434212]), np.array([ 0.24675064, -0.34005255]), np.array([0., 0.]), np.array([0., 0.]), np.array([0., 0.]), np.array([0.1187884 , 0.00154338]), np.array([-4.77536493, -3.46825118]), np.array([-0.92991243, -0.23784557]), np.array([0., 0.])]

        form_factors = ['Co','Ni']
        energy = [400,600,625,641.51,648.25,800]

        test_case = []
        for ff in form_factors:
            recast = [val for val in mm.find_form_factor(ff, energy, True)]
            test_case += recast

        my_sum = 0
        for idx in range(len(test_case)):
            my_sum += sum(solution[idx] - test_case[idx])

        self.assertFalse(my_sum == 0)

    def test_MOC(self):
        pass
    def test_magnetic_optical_constant(self):
        pass
    def test_IoR(self):
        pass
    def test_index_of_refraction(self):
        pass

if __name__ == '__main__':
    unittest.main()