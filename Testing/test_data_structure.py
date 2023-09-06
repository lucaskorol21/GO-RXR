import unittest
from data_structure import *

class TestMaterialStructure(unittest.TestCase):
    def test_getTitleInfo(self):
        # testing function that retrieves information from ReMagX data names
        tests = ['1_A_800.0_S', '57_A_515.2_P','42-43_A_55.5_S-P_Asymm', '75_E499.93_Th20.0_S',
                 '76_E592.1_Th25.2_P', '75-76_E499.93_Th5.0_S-P_Asymm', '88_A_640.2_R',
                 '99_A_642.2_L', '100-99_A_642.2_R-L_Asymm', '103_E640.2_Th25.0_L',
                 '104_E620.0_Th22.2_R', '103-104_E630.0_Th55.0_L-R_Asymm']
        scanNumber = ['1','57','42-43','75','76','75-76','88','99','100-99', '103', '104','103-104']
        scanType = ['Reflectivity','Reflectivity', 'Reflectivity','Energy','Energy','Energy','Reflectivity',
                    'Reflectivity', 'Reflectivity', 'Energy', 'Energy', 'Energy']
        energy = ['800.0','515.2', '55.5','E499.93','E592.1','E499.93','640.2', '642.2','642.2','E640.2', 'E620.0', 'E630.0']
        polarization = ['S','P', 'AL','S','P','AL', 'RC', 'LC', 'AC', 'LC','RC', 'AC']
        angle = [None, None, None,'20.0','25.2','5.0', None, None, None, '25.0','22.2', '55.0']

        for i, test in enumerate(tests):
            num, ty, E, pol, A = getTitleInfo(test)
            self.assertEqual(num, scanNumber[i])
            self.assertEqual(ty, scanType[i])
            self.assertEqual(E, energy[i])
            self.assertEqual(pol,polarization[i])
            self.assertEqual(A, angle[i])

    def test_ReadReMagX(self):
        pass

    def test_evaluateList(self):
        # tests to make sure that a simple list can be evaluated
        tests = ['[0,1]', '[0.01, 0.52, 20, 800]', '[0.01   , 0.21 , 0.55, 2]', '[0 2 6 1313 22]', '[0    1  2 0.02]']
        solutions = [[0.0, 1.0], [0.01, 0.52, 20.0, 800.0], [0.01, 0.21, 0.55, 2.0], [0.0, 2.0, 6.0, 1313.0, 22.0],
                     [0.0, 1.0, 2.0, 0.02]]

        for i, test in enumerate(tests):
            value = evaluate_list(test)
            self.assertListEqual(value, solutions[i])

    def test_find_parameter_bound(self):
        # seperate the parameter boundaries
        tests = ['[1,[0,1]]', '[0.5, [ 0.25, 1]], [10, [5 ,12]]', '[0.028, [ 0.01, 0.5]], [4, [2 ,11.5]], [6, [5 ,7]]',
                 '[0.5, [ 0.25, 1]], [10, [5 ,12]], [0.01, [ 0.25, 1]], [4, [5 ,12]]',
                 '[6, [2, 100]], [55, [2.5 ,55.5]]']

        solution = [['[1,[0,1]]'], ['[0.5, [ 0.25, 1]]','[10, [5 ,12]]'], ['[0.028, [ 0.01, 0.5]]','[4, [2 ,11.5]]','[6, [5 ,7]]'],
                  ['[0.5, [ 0.25, 1]]','[10, [5 ,12]]','[0.01, [ 0.25, 1]]','[4, [5 ,12]]'],
                  ['[6, [2, 100]]','[55, [2.5 ,55.5]]']]

        for i,test in enumerate(tests):
            value = find_parameter_bound(test)
            self.assertListEqual(value, solution[i])


    def test_find_each_bound(self):
        tests = ['[[0.01,0.5],[0.5,0.9]],[[0.01,0.9]]', '[[200,255],[255,600],[650,800]]',
                 '[[255.5 ,275]],[[ 800,895 ],[900, 921.5]]',
                 '[[0.01,0.52],[0.52,0.65]], [[0.55,0.5678778]], [[125.5,555],[587,598.9],[698,755.5]]']
        solution = [['[[0.01,0.5],[0.5,0.9]]', '[[0.01,0.9]]'], ['[[200,255],[255,600],[650,800]]'],
                    ['[[255.5 ,275]]','[[ 800,895 ],[900, 921.5]]'],
                    ['[[0.01,0.52],[0.52,0.65]]', '[[0.55,0.5678778]]','[[125.5,555],[587,598.9],[698,755.5]]']]

        for i,test in enumerate(tests):
            value = find_each_bound(test)
            self.assertListEqual(value, solution[i])

        pass

    def test_find_the_bound(self):
        tests = ['[[0.01,0.5],[0.5,0.9]]', '[[0.01,0.9]]', '[[200,255],[255,600],[650,800]]',
                 '[[255.5 ,275]]','[[ 800,895 ],[900, 921.5]]',
                 '[[0.01,0.52],[0.52,0.65]]', '[[0.55,0.5678778]]','[[125.5,555],[587,598.9],[698,755.5]]']

        solutions = [[[0.01,0.5],[0.5,0.9]],[[0.01,0.9]],[[200.0,255.0],[255.0,600.0],[650.0,800.0]],
                     [[255.5, 275.0]], [[ 800.0,895.0],[900.0, 921.5]],
                     [[0.01,0.52],[0.52,0.65]], [[0.55,0.5678778]], [[125.5,555],[587.0,598.9],[698.0,755.5]]]

        for i,test in enumerate(tests):
            value = find_the_bound(test)
            self.assertListEqual(value,solutions[i])

    def test_evaluated_bound(self):
        pass

    def test_evaluate_weights(self):
        pass

    def test_find_parameter_values(self):
        # retrieve the parameter boundary values
        tests = ['[1,[0,1]]','[0.5, [ 0.25, 1]]', '[10, [5 ,12]]', '[0.028, [ 0.01, 0.5]]', '[4, [2 ,11.5]]',
                 '[6, [5 ,7]]', '[0.5, [ 0.25, 1]]', '[10, [5 ,12]]', '[0.01, [ 0.25, 1]]', '[4, [5 ,12]]',
                 '[6, [2, 100]]', '[55, [2.5 ,55.5]]']

        solutions = [[1.0, [0.0, 1.0]],[0.5, [0.25, 1.0]], [10.0, [5.0, 12.0]],[0.028, [0.01, 0.5]],[4.0, [2.0, 11.5]],
                     [6.0, [5.0, 7.0]],[0.5,  [0.25, 1.0]], [10.0, [5.0, 12.0]], [0.01, [0.25, 1.0]], [4.0, [5.0, 12.0]],
                     [6.0, [2.0, 100.0]], [55.0, [2.5 ,55.5]]]

        for i, test in enumerate(tests):
            value = find_parameter_values(test)
            self.assertListEqual(value, solutions[i])


    def test_evaluate_parameters(self):
        # retrieve parameter bounds from string
        tests = ['[[1,[0,1]]]', '[[0.5, [ 0.25, 1]], [10, [5 ,12]]]', '[[0.028, [ 0.01, 0.5]], [4, [2 ,11.5]], [6, [5 ,7]]]',
                 '[[0.5, [ 0.25, 1]], [10, [5 ,12]], [0.01, [ 0.25, 1]], [4, [5 ,12]]]',
                 '[[6, [2, 100]], [55, [2.5 ,55.5]]]']

        solution = [[[1.0,[0.0,1.0]]], [[0.5, [ 0.25, 1.0]], [10.0, [5.0, 12.0]]],
                    [[0.028, [ 0.01, 0.5]], [4.0, [2.0 ,11.5]], [6.0, [5.0 ,7.0]]],
                    [[0.5, [ 0.25, 1.0]], [10.0, [5.0 ,12.0]], [0.01, [ 0.25, 1.0]], [4.0, [5.0 ,12.0]]],
                    [[6, [2, 100.0]], [55.0, [2.5 ,55.5]]]]

        for i, test in enumerate(tests):
            value = evaluate_parameters(test)
            self.assertListEqual(value, solution[i])