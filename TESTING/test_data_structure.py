import unittest
from data_structure import *

# This test script can be executed by inputting
#  ->  python -m unittest -v Testing/test_data_structure.py
# into the terminal

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
        filename = 'Pim7uc.all'
        filename2 = 'Pim7uc.pkl'
        if os.getcwd().split('\\')[-1] == 'Testing':
            my_path = os.getcwd() + '/test_data/' + filename
            my_path_2 = os.getcwd() + '/test_data/' + filename2
        else:
            my_path = os.getcwd() + '/Testing/test_data/' + filename
            my_path_2 = os.getcwd() + '/Testing/test_data/' + filename2

        data_info, data_dict = Read_ReMagX(my_path)

        data_info_solution = [[2, 'Reflectivity', '8_399.39_S'], [3, 'Reflectivity', '13_499.93_S'],
                              [4, 'Reflectivity', '18_600.18_S'], [5, 'Reflectivity', '23_700.14_S'],
                              [6, 'Reflectivity', '28_799.82_S'], [7, 'Reflectivity', '33_899.22_S'],
                              [8, 'Reflectivity', '99_799.82_RC'], [9, 'Reflectivity', '38_455.73_S'],
                              [10, 'Reflectivity', '39_455.73_P'], [11, 'Reflectivity', '38-39_455.73_AL_Asymm'],
                              [12, 'Reflectivity', '42_459.75_S'], [13, 'Reflectivity', '43_459.75_P'],
                              [14, 'Reflectivity', '42-43_459.75_AL_Asymm'], [15, 'Reflectivity', '48_640.2_S'],
                              [16, 'Reflectivity', '49_640.2_P'], [17, 'Reflectivity', '48-49_640.2_AL_Asymm'],
                              [18, 'Reflectivity', '52_642.2_S'], [19, 'Reflectivity', '53_642.2_P'],
                              [20, 'Reflectivity', '52-53_642.2_AL_Asymm'], [21, 'Reflectivity', '58_833.85_S'],
                              [22, 'Reflectivity', '59_836.04_S'], [23, 'Reflectivity', '72_640.2_LC'],
                              [24, 'Reflectivity', '73_640.2_RC'], [25, 'Reflectivity', '72-73_640.2_AC_Asymm'],
                              [26, 'Reflectivity', '78_642.2_LC'], [27, 'Reflectivity', '79_642.2_RC'],
                              [28, 'Reflectivity', '78-79_642.2_AC_Asymm'], [29, 'Reflectivity', '104_640.2_LC'],
                              [30, 'Reflectivity', '105_640.2_RC'], [31, 'Reflectivity', '104-105_640.2_AC_Asymm'],
                              [32, 'Reflectivity', '110_642.2_LC'], [33, 'Reflectivity', '111_642.2_RC'],
                              [34, 'Reflectivity', '110-111_642.2_AC_Asymm'], [35, 'Energy', '54_E429.58_Th5.0_S'],
                              [36, 'Energy', '55_E429.58_Th5.0_P'], [37, 'Energy', '54-55_E429.58_Th5.0_AL_Asymm'],
                              [38, 'Energy', '44_E429.58_Th10.0_S'], [39, 'Energy', '45_E429.58_Th10.0_P'],
                              [40, 'Energy', '44-45_E429.58_Th10.0_AL_Asymm'], [41, 'Energy', '34_E429.58_Th15.0_S'],
                              [42, 'Energy', '35_E429.58_Th15.0_P'], [43, 'Energy', '34-35_E429.58_Th15.0_AL_Asymm'],
                              [44, 'Energy', '24_E429.58_Th20.0_S'], [45, 'Energy', '25_E429.58_Th20.0_P'],
                              [46, 'Energy', '24-25_E429.58_Th20.0_AL_Asymm'], [47, 'Energy', '14_E429.58_Th25.0_S'],
                              [48, 'Energy', '15_E429.58_Th25.0_P'], [49, 'Energy', '14-15_E429.58_Th25.0_AL_Asymm'],
                              [50, 'Energy', '9_E429.58_Th30.0_S'], [51, 'Energy', '10_E429.58_Th30.0_P'],
                              [52, 'Energy', '9-10_E429.58_Th30.0_AL_Asymm'], [53, 'Energy', '19_E429.58_Th35.0_S'],
                              [54, 'Energy', '20_E429.58_Th35.0_P'], [55, 'Energy', '19-20_E429.58_Th35.0_AL_Asymm'],
                              [56, 'Energy', '29_E429.58_Th40.0_S'], [57, 'Energy', '30_E429.58_Th40.0_P'],
                              [58, 'Energy', '29-30_E429.58_Th40.0_AL_Asymm'], [59, 'Energy', '40_E429.58_Th45.0_S'],
                              [60, 'Energy', '50_E429.58_Th50.0_S'], [61, 'Energy', '60_E429.58_Th55.0_S'],
                              [62, 'Energy', '56_E600.18_Th5.0_S'], [63, 'Energy', '57_E600.18_Th5.0_P'],
                              [64, 'Energy', '56-57_E600.18_Th5.0_AL_Asymm'], [65, 'Energy', '46_E600.18_Th10.0_S'],
                              [66, 'Energy', '47_E600.18_Th10.0_P'], [67, 'Energy', '46-47_E600.18_Th10.0_AL_Asymm'],
                              [68, 'Energy', '36_E600.18_Th15.0_S'], [69, 'Energy', '37_E600.18_Th15.0_P'],
                              [70, 'Energy', '36-37_E600.18_Th15.0_AL_Asymm'], [71, 'Energy', '26_E600.18_Th20.0_S'],
                              [72, 'Energy', '27_E600.18_Th20.0_P'], [73, 'Energy', '26-27_E600.18_Th20.0_AL_Asymm'],
                              [74, 'Energy', '16_E600.18_Th25.0_S'], [75, 'Energy', '17_E600.18_Th25.0_P'],
                              [76, 'Energy', '16-17_E600.18_Th25.0_AL_Asymm'], [77, 'Energy', '11_E600.18_Th30.0_S'],
                              [78, 'Energy', '12_E600.18_Th30.0_P'], [79, 'Energy', '11-12_E600.18_Th30.0_AL_Asymm'],
                              [80, 'Energy', '21_E600.18_Th35.0_S'], [81, 'Energy', '22_E600.18_Th35.0_P'],
                              [82, 'Energy', '21-22_E600.18_Th35.0_AL_Asymm'], [83, 'Energy', '31_E600.18_Th40.0_S'],
                              [84, 'Energy', '32_E600.18_Th40.0_P'], [85, 'Energy', '31-32_E600.18_Th40.0_AL_Asymm'],
                              [86, 'Energy', '41_E600.18_Th45.0_S'], [87, 'Energy', '51_E600.18_Th50.0_S'],
                              [88, 'Energy', '61_E600.18_Th55.0_S'], [89, 'Energy', '88_E600.18_Th5.0_LC'],
                              [90, 'Energy', '89_E600.18_Th5.0_RC'], [91, 'Energy', '88-89_E600.18_Th5.0_AC_Asymm'],
                              [92, 'Energy', '86_E600.18_Th10.0_LC'], [93, 'Energy', '87_E600.18_Th10.0_RC'],
                              [94, 'Energy', '86-87_E600.18_Th10.0_AC_Asymm'], [95, 'Energy', '82_E600.18_Th15.0_LC'],
                              [96, 'Energy', '83_E600.18_Th15.0_RC'], [97, 'Energy', '82-83_E600.18_Th15.0_AC_Asymm'],
                              [98, 'Energy', '76_E600.18_Th20.0_LC'], [99, 'Energy', '77_E600.18_Th20.0_RC'],
                              [100, 'Energy', '76-77_E600.18_Th20.0_AC_Asymm'], [101, 'Energy', '70_E600.18_Th25.0_LC'],
                              [102, 'Energy', '71_E600.18_Th25.0_RC'], [103, 'Energy', '70-71_E600.18_Th25.0_AC_Asymm'],
                              [104, 'Energy', '68_E600.18_Th30.0_LC'], [105, 'Energy', '69_E600.18_Th30.0_RC'],
                              [106, 'Energy', '68-69_E600.18_Th30.0_AC_Asymm'], [107, 'Energy', '74_E600.18_Th35.0_LC'],
                              [108, 'Energy', '75_E600.18_Th35.0_RC'], [109, 'Energy', '74-75_E600.18_Th35.0_AC_Asymm'],
                              [110, 'Energy', '80_E600.18_Th40.0_LC'], [111, 'Energy', '81_E600.18_Th40.0_RC'],
                              [112, 'Energy', '80-81_E600.18_Th40.0_AC_Asymm'], [113, 'Energy', '84_E600.18_Th45.0_LC'],
                              [114, 'Energy', '85_E600.18_Th45.0_RC'], [115, 'Energy', '84-85_E600.18_Th45.0_AC_Asymm'],
                              [116, 'Energy', '120_E600.18_Th5.0_LC'], [117, 'Energy', '121_E600.18_Th5.0_RC'],
                              [118, 'Energy', '120-121_E600.18_Th5.0_AC_Asymm'], [119, 'Energy', '118_E600.18_Th10.0_LC'],
                              [120, 'Energy', '119_E600.18_Th10.0_RC'], [121, 'Energy', '118-119_E600.18_Th10.0_AC_Asymm'],
                              [122, 'Energy', '114_E600.18_Th15.0_LC'], [123, 'Energy', '115_E600.18_Th15.0_RC'],
                              [124, 'Energy', '114-115_E600.18_Th15.0_AC_Asymm'], [125, 'Energy', '108_E600.18_Th20.0_LC'],
                              [126, 'Energy', '109_E600.18_Th20.0_RC'], [127, 'Energy', '108-109_E600.18_Th20.0_AC_Asymm'],
                              [128, 'Energy', '102_E600.18_Th25.0_LC'], [129, 'Energy', '103_E600.18_Th25.0_RC'],
                              [130, 'Energy', '102-103_E600.18_Th25.0_AC_Asymm'], [131, 'Energy', '100_E600.18_Th30.0_LC'],
                              [132, 'Energy', '101_E600.18_Th30.0_RC'], [133, 'Energy', '100-101_E600.18_Th30.0_AC_Asymm'],
                              [134, 'Energy', '106_E600.18_Th35.0_LC'], [135, 'Energy', '107_E600.18_Th35.0_RC'],
                              [136, 'Energy', '106-107_E600.18_Th35.0_AC_Asymm'], [137, 'Energy', '112_E600.18_Th40.0_LC'],
                              [138, 'Energy', '113_E600.18_Th40.0_RC'], [139, 'Energy', '112-113_E600.18_Th40.0_AC_Asymm'],
                              [140, 'Energy', '116_E600.18_Th45.0_LC'], [141, 'Energy', '117_E600.18_Th45.0_RC'],
                              [142, 'Energy', '116-117_E600.18_Th45.0_AC_Asymm'], [143, 'Energy', '62_E814.75_Th10.0_S'],
                              [144, 'Energy', '63_E814.75_Th25.0_S'], [145, 'Energy', '64_E499.93_Th5.0_S'],
                              [146, 'Energy', '65_E499.93_Th5.0_P'], [147, 'Energy', '64-65_E499.93_Th5.0_AL_Asymm'],
                              [148, 'Energy', '66_E499.93_Th20.0_S'], [149, 'Energy', '67_E499.93_Th20.0_P'],
                              [149, 'Energy', '66-67_E499.93_Th20.0_AL_Asymm']]

        self.assertListEqual(data_info, data_info_solution)

        with open(my_path_2, 'rb') as file:
            data_dict_solution = pickle.load(file)
        for key in data_dict.keys():
            for param in data_dict[key].keys():
                if param == 'Polarization':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])
                elif param == 'DatasetNumber':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])
                elif param == 'Background Shift':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])
                elif param == 'Scaling Factor':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])
                elif param == 'Energy':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])
                elif param == 'DataPoints':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])
                elif param == 'Data':
                    test_list = list(data_dict[key][param])
                    solution_list = list(data_dict_solution[key][param])
                    for i in range(len(test_list)):
                        self.assertListEqual(list(test_list[i]), list(solution_list[i]))
                elif param == 'Angle':
                    self.assertEqual(data_dict[key][param], data_dict_solution[key][param])



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
        tests = ['[[[0.01,0.5],[0.5,0.9]],[[0.01,0.9]]]', '[[[200,255],[255,600],[650,800]]]',
                 '[[[255.5 ,275]],[[ 800,895 ],[900, 921.5]]]',
                 '[[[0.01,0.52],[0.52,0.65]], [[0.55,0.5678778]], [[125.5,555],[587,598.9],[698,755.5]]]']
        solutions = [[[[0.01,0.5],[0.5,0.9]],[[0.01,0.9]]], [[[200,255],[255,600],[650,800]]],
                     [[[255.5 ,275]],[[ 800,895 ],[900, 921.5]]],
                     [[[0.01,0.52],[0.52,0.65]], [[0.55,0.5678778]], [[125.5,555],[587,598.9],[698,755.5]]]]

        for i, test in enumerate(tests):
            value = evaluate_bounds(test)
            self.assertListEqual(value,solutions[i])

    def test_evaluate_weights(self):
        tests = ['[[1,2,3,4,5]]', '[[1.2,   3, 6.8, 10]]', '[[23,1.2,32.2],  [1, 3, 9.0]]',
                '[[1,2,3,4,5],[1,1,3]]', '[[3.3, 4.4], [1], [6.5,7.7]]',
                '[[1], [ 2 ], [3,     4]]']

        solutions = [[[1.0,2.0,3.0,4.0,5.0]], [[1.2, 3.0, 6.8, 10.0]], [[23.0,1.2,32.2],  [1.0, 3.0, 9.0]],
                [[1.0, 2.0, 3.0, 4.0, 5.0], [1.0, 1.0, 3.0]], [[3.3, 4.4], [1.0], [6.5,7.7]],
                [[1.0], [2.0], [3.0, 4.0]]]

        for i, test in enumerate(tests):
            value = evaluate_weights(test)
            self.assertListEqual(value, solutions[i])


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