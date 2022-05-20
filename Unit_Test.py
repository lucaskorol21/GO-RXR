
from material_structure import *

import matplotlib
import numpy as np

if __name__ == "__main__":
    # ------------------------------------------TEST 0---------------------------------------------------------------- #
    # ---------------------------------------Class: element----------------------------------------------------------- #
    Success = True

    test_00 = [{'input':['Fe', 1], 'output': ['Fe', 0, 0, 0, 0, 1, 1, [], 0, 0, None, 'Fe', None, None], 'reason': 0},
               {'input': ['Sr', 7], 'output': ['Sr', 0, 0, 0, 0, 7, 1, [], 0, 0, None, 'Sr', None, None], 'reason': 0},
               {'input': ['K', 1], 'output': ['K', 0, 0, 0, 0, 1, 1, [], 0, 0, None, 'K', None, None], 'reason': 0},
               {'input': ['C', 2], 'output': ['C', 0, 0, 0, 0, 2, 1, [], 0, 0, None, 'C', None, None], 'reason': 0},
               {'input':['Ga', 25], 'output': ['Ga', 0, 0, 0, 0, 25, 1, [], 0, 0, None, 'Ga', None, None], 'reason': 0},
               {'input':['V', 1], 'output': ['V', 0, 0, 0, 0, 1, 1, [], 0, 0, None, 'V', None, None], 'reason': 0},
               {'input':['N', 1], 'output': ['N', 0, 0, 0, 0, 1, 1, [], 0, 0, None, 'N', None, None], 'reason': 0},
               {'input':['Hg', 184], 'output': ['Hg', 0, 0, 0, 0, 184, 1, [], 0, 0, None, 'Hg', None, None], 'reason': 0}]

    for test in test_00:
        test_element = element(test['input'][0], test['input'][1])

        if test_element.name != test['output'][0]:
            print('Element name not initialized properly')
            Success = False
        elif test_element.molar_mass != test['output'][1]:
            print('Element molar mass not initialized properly')
            Success = False
        elif test_element.density != test['output'][2]:
            print('Element density not initialized properly')
            Success = False
        elif test_element.thickness != test['output'][3]:
            print('Element thickness not initialized properly')
            Success = False
        elif test_element.roughness != test['output'][4]:
            print('Element roughness not initialized properly')
            Success = False
        elif test_element.stoichiometry != test['output'][5]:
            print('Element stoichiometry not initialized properly')
            Success = False
        elif test_element.poly_ratio != test['output'][6]:
            print('Element polymorphous ratio not initialized properly')
            Success = False
        elif test_element.polymorph != test['output'][7]:
            print('Element polymorph not initialized properly')
            Success = False
        elif test_element.gamma != test['output'][8]:
            print('Element magnetization direction theta not initialized properly')
            Success = False
        elif test_element.phi != test['output'][9]:
            print('Element magnetization direction phi not initialized properly')
            Success = False
        elif test_element.mag_density != test['output'][10]:
            print('Element magnetization density not initialized properly')
            Success = False
        elif test_element.scattering_factor != test['output'][11]:
            print('Element scattering factor not initialized properly')
            Success = False
        elif test_element.mag_scattering_factor != test['output'][12]:
            print('Element magnetization scattering factor not initialized properly')
            Success = False
        elif test_element.position != test['output'][13]:
            print('Element position not initialized properly')
            Success = False

    print('0. element class testing complete')
    # ------------------------------------------TEST 1---------------------------------------------------------------- #
    # -----------------------------------------get_number------------------------------------------------------------- #

    test_01 = [{'input': '0', 'output': ['',0], 'reason': 0},
               {'input': '99', 'output': ['',99], 'reason': 0},
               {'input': '2O', 'output': ['O',2], 'reason': 0},
               {'input': '546213La14MnO', 'output': ['La14MnO', 546213], 'reason': 0},
               {'input': '74Cufhgbfdjvbkjhf', 'output': ['Cufhgbfdjvbkjhf', 74], 'reason': 0}]

    for test in test_01:
        result1, result2 = get_number(test['input'])

        if result1 != test['output'][0]:
            print('String ' + test['input'] + ' incorrect ')
            Success = False
        elif result2 != test['output'][1]:
            print('Number ' + test['input']+ 'incorrect')
            Success = False

    print('1. get_number testing complete')

    # ------------------------------------------TEST 02--------------------------------------------------------------- #
    # ----------------------------------------checkstring------------------------------------------------------------- #

    test_02 = [{'input':'LaMn3', 'output': ['Mn3', ['La',1]], 'reason':0},
               {'input':'La2MnO3', 'output': ['MnO3', ['La',2]], 'reason':0},
               {'input':'SO3', 'output': ['O3', ['S',1]], 'reason':0},
               {'input':'As77B0N3', 'output': ['B0N3', ['As',77]], 'reason':0},
               {'input':'P10KO7N3', 'output': ['KO7N3', ['P',10]], 'reason':0},
               {'input':'SrTi7Ni', 'output': ['Ti7Ni', ['Sr',1]], 'reason':0},
               {'input':'Hg4KZn2', 'output': ['KZn2', ['Hg',4]], 'reason':0}]

    for test in test_02:
        result1, result2 = checkstring(test['input'])

        if result1 != test['output'][0]:
            print('The formula ' + test['input'] + ' is reduced to ' + result1 + ' and not ' + test['output'][0])
            Success = False

        if result2[0] != test['output'][1][0]:
            print('Element not removed properly.')
            Success = False

        if result2[1] != test['output'][1][1]:
            print('Stochicometry not correct.')
            Success = False

    print('2. checkstring testing complete')

    # ------------------------------------------TEST 03--------------------------------------------------------------- #
    # ----------------------------------------atomic_mass------------------------------------------------------------- #

    test_03 = [{'input':'H', 'output':1.00797, 'reason': 0},
               {'input':'Be', 'output':9.01218, 'reason': 0},
               {'input':'O', 'output':15.9994, 'reason': 0},
               {'input':'La', 'output':138.9055, 'reason': 0},
               {'input':'Mn', 'output':54.938, 'reason': 0},
               {'input':'Sr', 'output':87.62, 'reason': 0},
               {'input':'Ti', 'output':47.9, 'reason': 0},
               {'input':'Zr', 'output':91.22, 'reason': 0},
               {'input':'K', 'output':39.0983, 'reason': 0},
               {'input':'Hg', 'output':200.59, 'reason': 0},
               {'input':'Hf', 'output':178.49, 'reason': 0}]


    for test in test_03:
        result1 = atomic_mass(test['input'])
        if result1 != test['output']:
            print('Expecting molar mass of ' + str(test['output'])+ 'g/mol for '+ test['input']+ ', but got the value of ' + str(result1) + 'g/mol.')

    print('3. atomic mass testing complete')
    # ------------------------------------------TEST 04--------------------------------------------------------------- #
    # -------------------------------------perovskite_density--------------------------------------------------------- #

    test_04 = [{'input': 'SrTiO3', 'output': 5.12, 'reason': 0},
               {'input': 'CaTiO3', 'output': 4.10, 'reason': 0},
               {'input': 'PbTiO3', 'output': 7.52, 'reason': 0},
               {'input': 'BiFeO3', 'output': 8.22, 'reason': 0},
               {'input': 'LaYbO3', 'output': 8.08, 'reason': 0},
               {'input': 'LaMnO3', 'output': 6.52, 'reason': 0},
               {'input': 'AlYO3', 'output': 5.35, 'reason': 0},
               {'input': 'LaAlO3', 'output': 6.52, 'reason': 0}]

    for test in test_04:
        result1 = perovskite_density(test['input'])

        if result1 != test['output']:
            print('Density of ' + str(test['output']) + 'g/cm^3  for ' + test['input'] + '  not ' + str(result1) + 'g/cm^3.')
            Success = False

    print('4. pervoskite_density testing complete')



    # ------------------------------------------TEST 05--------------------------------------------------------------- #
    # -------------------------------------find_stoichiometry--------------------------------------------------------- #
    test_05 = [{'input': 'FeNeO3', 'output': {'Fe': element('Fe', 1), 'Ne': element('Ne',1), 'O': element('O', 3)}, 'reason': 'FeNeO3 Incorrect'},
               {'input': 'La2222Mn7784848O', 'output': {'La': element('La',2222), 'Mn': element('Mn', 7784848), 'O': element('O',1)}, 'reason': 'La2222Mn7784848O incorrect'},
               {'input': 'FeMn2', 'output': {'Fe':element('Fe',1), 'Mn': element('Mn',2)}, 'reason': 'FeMn2 incorrect'},
               {'input': 'CoZn', 'output': {'Co':element('Co',1), 'Zn': element('Zn',1)}, 'reason': 'CoZn incorrect'},
               {'input': 'BaLa6K2', 'output': {'Ba':element('Ba',1), 'La': element('La',6), 'K': element('K',2)}, 'reason': 'BLa6K2 incorrect'},
               {'input': 'CuK7MnO', 'output': {'Cu':element('Cu',1), 'K': element('K',7), 'Mn': element('Mn',1), 'O': element('O',1)}, 'reason': 'CuK7MnO incorrect'},
               {'input': 'FOs2', 'output': {'F': element('F', 1), 'Os': element('Os', 2)}, 'reason': 'FOs2 incorrect'},
               {'input': 'CZn', 'output': {'C': element('C', 1), 'Zn': element('Zn', 1)}, 'reason': 'CZn incorrect'},
               {'input': 'BLa6N2', 'output': {'B': element('B', 1), 'La': element('La', 6), 'N': element('N', 2)}, 'reason': 'BLa6N2 incorrect'},
               {'input': 'VK7MnS', 'output': {'V': element('V', 1), 'K': element('K', 7), 'Mn': element('Mn', 1), 'S': element('S', 1)}, 'reason': 'VK7MnS incorrect'}]


    for test in test_05:
        result1,result2 = find_stoichiometry(test['input'])

        # Check if there is the same number of dictionary elements
        if len(result1) != len(test['output']):
            print('Incorrect number of dictionary elements for finding formula stoichiometry.')
            print(test['reason'])
            Success = False

        # Testing keys
        if set(list(result1.keys())) != set(list(test['output'].keys())):
            print('Keys for find_stoichiometry not as expected.')
            print(test['reason'])
            Success = False

        # Test key order
        if list(result1.keys()) != list(test['output'].keys()):
            print('Keys and values not in correct order.')
            print(test['reason'])
            Success = False

        for check in list(result1.keys()):
            answer = test['output'][check]
            test_function = result1[check]

            # Checking that name input is correct in element class
            if answer.name != test_function.name:
                print('Name incorrect in element class.')
                print(test['reason'])
                Success = False

            # Checking that stoichiometry input is correct in element class
            if answer.stoichiometry != test_function.stoichiometry:
                print('Stoichiometry of element incorrect.')
                print(test['reason'])
                Success = False


    print('5. find_stochiometry testing complete')

    # ---------------------------------------------TEST--------------------------------------------------------------- #
    # -----------------------------------------Class: slab------------------------------------------------------------ #

    # -------------------------------------------TEST 06-------------------------------------------------------------- #
    # ---------------------------------------Slab initialization------------------------------------------------------ #

    test_06 = [{'input': 1, 'output': 1, 'reason': 'single digit'},
               {'input': 4, 'output': 4, 'reason': 'single digit' },
               {'input': 9, 'output': 9, 'reason': 'single digit' },
               {'input': 10, 'output': 10, 'reason': 'double digit' },
               {'input': 27, 'output': 27, 'reason': 'double digit'},
               {'input': 77, 'output': 77, 'reason': 'double digit'},
               {'input': 98, 'output': 98, 'reason': 'double digit'},
               {'input': 156, 'output': 156, 'reason': 'triple digit'},
               {'input': 377, 'output': 377, 'reason': 'triple digit'},
               {'input': 1025, 'output': 1025, 'reason': 'quadruple digit'},
               {'input': -5, 'output': 5, 'reason': 'negative digit'},
               {'input': -155, 'output': 155, 'reason': 'negative digit'},
               {'input': -2098, 'output': 2098, 'reason': 'negative digit'}]

    for test in test_06:
        sample_test = slab(test['input'])

        if len(sample_test.structure) != test['output']:
            print('Incorrect number of layers initialized for ' + test['reason'])

    print('6. slab initialization testing complete')

    # -------------------------------------------TEST 07-------------------------------------------------------------- #
    # -------------------------------------------addlayer------------------------------------------------------------- #

    test_07 = [{'input': [0, 'LaMnO3', 50, 2.5,2,[True,True, False]], 'output': ['LaMnO', [1,1,3],50, 2.5,2,[True,True, False]], 'reason': 1},
               {'input': [1, 'LaAlO3', 2, 0.5,2,[True,False, False]], 'output': ['LaAlO', [1,1,3], 2, 0.5,2,[True,False, False]], 'reason': 1},
               {'input': [2, 'SrTiO3', 17.5, 1,0.1,[True,True, True]], 'output': ['SrTiO', [1,1,3], 17.5, 1,0.1,[True,True, True]], 'reason': 1},
               {'input': [3, 'LaMnO3', 5.55, 2.5,3,[False,False, False]], 'output': ['LaMnO', [1,1,3], 5.55, 2.5,3,[False,False, False]], 'reason': 1},
               {'input': [4, 'LaMnO3', 101, 7.7,2.2,[True,True, False]], 'output': ['LaMnO', [1,1,3], 101, 7.7,2.2,[True,True, False]], 'reason': 1},
               {'input': [5, 'SrTiO3', 111.111, 2.5,2,[True,False, False]], 'output': ['SrTiO', [1,1,3], 111.111, 2.5,2,[True,False, False]], 'reason': 1},
               {'input': [6, 'LaAlO3', 1, 2.5,0.66,[False,True, True]], 'output': ['LaAlO', [1,1,3], 1, 2.5,0.66,[False,True, True]], 'reason': 1},
               {'input': [7, 'LaYbO3', 6.67, 9,2,[False,True, False]], 'output': ['LaYbO', [1,1,3], 6.67, 9,2,[False,True, False]], 'reason': 1},
               {'input': [8, 'AlYO3', 42, 2.5,1.89,[False,False, True]], 'output': ['AlYO', [1,1,3], 42, 2.5,1.89,[False,False, True]], 'reason': 1},
               {'input': [9, 'AlYO3', 42.1, 7.5,1.1,[True,False, True]], 'output': ['AlYO', [1,1,3], 42.1, 7.5,1.1,[True,False, True]], 'reason': 1}]

    myelements = ['La', 'Sr','Mn', 'Al', 'Ti',' Yb', 'Y', 'O']
    sample_test = slab(10)

    for layer in test_07:
        sample_test.addlayer(layer['input'][0], layer['input'][1], layer['input'][2], density=layer['input'][3], roughness=layer['input'][4], link=layer['input'][5])

    for lay in range(len(sample_test.structure)):
        test_layer = test_07[lay]
        layer = sample_test.structure[lay]

        element_name = ''.join(list(layer.keys()))

        if sample_test.link[lay] != test_layer['output'][5]:
            print('Element link incorrect.')
            Success = False

        if element_name != test_layer['output'][0]:
            print('Element names do not match original chemical formula')
            Success = False

        idx = 0
        for ele in list(layer.keys()):
            if layer[ele].stoichiometry != test_layer['output'][1][idx]:
                print('Stocihiometry in layer ' + str(lay) + ' is incorrect for' + ele)
                Success = False

            if layer[ele].thickness != test_layer['output'][2]:
                print('Thickness in layer ' + str(lay) + ' is incorrect for' + ele)
                Success = False

            if layer[ele].density != test_layer['output'][3]:
                print('Density in layer ' + str(lay) + ' is incorrect for' + ele)
                Success = False

            if layer[ele].roughness != test_layer['output'][4]:
                print('Roughness in layer ' + str(lay) + ' is incorrect for' + ele)
                Success = False

            idx = idx + 1

    print('7. addlayer testing complete')

    # -------------------------------------------TEST 08-------------------------------------------------------------- #
    # -------------------------------------------polymorph------------------------------------------------------------ #

    test_08 = [{'input': [0,'Mn', ['Mn2+', 'Mn3+'],[0.5,0.5], ['Mn','Fe']], 'output':['Mn', ['Mn2+', 'Mn3+'],[0.5,0.5], ['Mn','Fe']], 'reason': 1},
               {'input': [2,'Ti', ['Ti3+', 'Ti4+'],[0.1,0.9], ['K','Ge']], 'output':['Ti', ['Ti3+', 'Ti4+'],[0.1,0.9], ['K','Ge']], 'reason': 1},
               {'input': [3, 'Mn', ['Mn2+', 'Mn3+'], [0.66, 0.34], ['Mn', 'Fe']],'output': ['Mn', ['Mn2+', 'Mn3+'], [0.66, 0.34], ['Mn', 'Fe']], 'reason': 1},
               {'input': [4,'Mn', ['Mn2+', 'Mn3+'],[0.85,0.15], ['Mn','Fe']], 'output':['Mn', ['Mn2+', 'Mn3+'],[0.85,0.15], ['Mn','Fe']], 'reason': 1},
               {'input': [5,'Ti', ['Ti3+', 'Ti4+'],[0.74,0.26], ['K','Ge']], 'output':['Ti', ['Ti3+', 'Ti4+'],[0.74,0.26], ['K','Ge']], 'reason': 1}]

    for test in test_08:
        sample_test.polymorphous(test['input'][0], test['input'][1], test['input'][2], test['input'][3], sf=test['input'][4])

    idx = 0
    for lay in range(len(sample_test.structure)):
        layer = sample_test.structure[lay]
        if idx <=len(test_08)-1:
            if test_08[idx]['input'][0] == lay:
                test = test_08[idx]
                idx = idx + 1
        for ele in list(layer.keys()):
            if ele in list(sample_test.poly_elements.keys()):  # polymorphous element
                if layer[ele].polymorph != test['output'][1]:
                    print('Polymorphous element ' + ele + ' incorrectly initialized')
                    Success = False

                for rat in range(len(layer[ele].poly_ratio)):
                    if layer[ele].poly_ratio[rat] != test['output'][2][rat]:
                        print('Polymorphous element ratio ' + ele + ' incorrect initialized')
                        Success = False
                if layer[ele].scattering_factor != test['output'][3]:
                    print('Polymorphous element scattering factor ' + ele + ' incorrect initialized')
                    Success = False
            else:
                if layer[ele].polymorph != []:
                    print('Non-polymorphous element '+ ele + ' incorrectly initialized' )
                    Success = False
                if layer[ele].poly_ratio != 1:
                    print('Non-polymorphous element ratio ' + ele + ' incorrect initialized')
                    Success = False
                if layer[ele].scattering_factor != ele:
                    print('Non-polymorphous element scattering factor ' + ele + ' incorrect initialized')
                    Success = False


    print('8. polymorph testing complete')

    # -------------------------------------------TEST 09-------------------------------------------------------------- #
    # ----------------------------------------magnetization----------------------------------------------------------- #

    test_09 = [{'input': [0, ['Mn2+', 'Mn3+'], [5,6],['Co', 'Fe'], 0, 0], 'output': [ [5,6],['Co', 'Fe'], 0, 0], 'reason': 1},
               {'input': [3, ['Mn3+', 'Mn2+'], [0,3],['Fe', 'Co'], 45, 60], 'output': [[3,0],['Co', 'Fe'], 45, 60], 'reason': 1},
               {'input': [4, ['Mn2+', 'Mn3+'], [4.4,7.9],['Co', 'Fe'], 16, 181], 'output': [[4.4,7.9],['Co', 'Fe'], 16, 181], 'reason': 1}]

    for test in test_09:
        sample_test.magnetization(test['input'][0], test['input'][1], test['input'][2], test['input'][3], gamma=test['input'][4], phi=test['input'][5])

    idx = 0
    for lay in range(len(sample_test.structure)):
        layer = sample_test.structure[lay]
        if idx <=len(test_09)-1:
            if test_09[idx]['input'][0] == lay:
                test = test_09[idx]
                idx = idx + 1

        for ele in list(layer.keys()):
            if ele in list(sample_test.mag_elements.keys()):  # polymorphous element
                for md in range(len(layer[ele].mag_density)):
                    if layer[ele].mag_density[md] != test['output'][0][md]:
                        print('Magnetic density of element ' + ele + ' incorrectly initialized')
                        Success = False
                if layer[ele].mag_scattering_factor != test['output'][1]:
                    print('Magnetic scattering factor of element ' + ele + ' incorrect initialized')
                    Success = False
                if layer[ele].gamma != test['output'][2]:
                    print('Magnetic phi of element ' + ele + ' incorrect initialized')
                    Success = False
                if layer[ele].phi != test['output'][3]:
                    print('Magnetic theta of element scattering factor ' + ele + ' incorrect initialized')
                    Success = False
            else:
                if layer[ele].mag_density != None:
                    print('Non-magnetic element ' + ele + ' incorrectly initialized')
                    Success = False
                if layer[ele].mag_scattering_factor != None:
                    print('Non-magnetic scattering factor of element ' + ele + ' incorrect initialized')
                    Success = False
                if layer[ele].gamma != 0:
                    print('Non-magnetic phi of element ' + ele + ' incorrect initialized')
                    Success = False
                if layer[ele].phi != 0:
                    print('Non-magnetic theta of element ' + ele + ' incorrect initialized')
                    Success = False

    print('9. magnetization testing complete')
    if Success:
        print()
        print()
        print('UNIT TESTING SUCCESSFUL :)')