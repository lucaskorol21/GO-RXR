
from material_structure import *
import matplotlib
import numpy as np

if __name__ == "__main__":
    # ------------------------------------------TEST 0---------------------------------------------------------------- #
    # ---------------------------------------Class: element----------------------------------------------------------- #

    test_00 = [{'input':['Fe', 1], 'output': ['Fe', 0, 0, 0, 0, 1, 1, [], None, None, 'Fe', None], 'reason': 0},
               {'input': ['Sr', 7], 'output': ['Sr', 0, 0, 0, 0, 7, 1, [], None, None, 'Sr', None], 'reason': 0},
               {'input': ['K', 1], 'output': ['K', 0, 0, 0, 0, 1, 1, [], None, None, 'K', None], 'reason': 0},
               {'input': ['C', 2], 'output': ['C', 0, 0, 0, 0, 2, 1, [], None, None, 'C', None], 'reason': 0},
               {'input':['Ga', 25], 'output': ['Ga', 0, 0, 0, 0, 25, 1, [], None, None, 'Ga', None], 'reason': 0},
               {'input':['V', 1], 'output': ['V', 0, 0, 0, 0, 1, 1, [], None, None, 'V', None], 'reason': 0},
               {'input':['N', 1], 'output': ['N', 0, 0, 0, 0, 1, 1, [], None, None, 'N', None], 'reason': 0},
               {'input':['Hg', 184], 'output': ['Hg', 0, 0, 0, 0, 184, 1, [], None, None, 'Hg', None], 'reason': 0}]

    for test in test_00:
        test_element = element(test['input'][0], test['input'][1])

        if test_element.name != test['output'][0]:
            print('Element name not initialized properly')
        elif test_element.molar_mass != test['output'][1]:
            print('Element molar mass not initialized properly')
        elif test_element.density != test['output'][2]:
            print('Element density not initialized properly')
        elif test_element.thickness != test['output'][3]:
            print('Element thickness not initialized properly')
        elif test_element.roughness != test['output'][4]:
            print('Element roughness not initialized properly')
        elif test_element.stoichiometry != test['output'][5]:
            print('Element stoichiometry not initialized properly')
        elif test_element.poly_ratio != test['output'][6]:
            print('Element polymorphous ratio not initialized properly')
        elif test_element.polymorph != test['output'][7]:
            print('Element polymorph not initialized properly')
        elif test_element.mag_dir != test['output'][8]:
            print('Element magnetization direction not initialized properly')
        elif test_element.mag_density != test['output'][9]:
            print('Element magnetization density not initialized properly')
        elif test_element.scattering_factor != test['output'][10]:
            print('Element scattering factor not initialized properly')
        elif test_element.mag_scattering_factor != test['output'][11]:
            print('Element magnetization scattering factor not initialized properly')

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
        elif result2 != test['output'][1]:
            print('Number ' + test['input']+ 'incorrect')

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

        if result2[0] != test['output'][1][0]:
            print('Element not removed properly.')

        if result2[1] != test['output'][1][1]:
            print('Stochicometry not correct.')

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

    test_04 = [{'input': 'SrTiO3', 'output': 4.81, 'reason': 0},
               {'input': 'CaTiO3', 'output': 3.98, 'reason': 0},
               {'input': 'PbTiO3', 'output': 7.52, 'reason': 0},
               {'input': 'BiFeO3', 'output': 8.22, 'reason': 0},
               {'input': 'LaYbO3', 'output': 8.08, 'reason': 0},
               {'input': 'LaMnO3', 'output': 6.50, 'reason': 0},
               {'input': 'AlYO3', 'output': 5.35, 'reason': 0},
               {'input': 'LaAlO3', 'output': 6.42, 'reason': 0}]

    for test in test_04:
        result1 = perovskite_density(test['input'])

        if result1 != test['output']:
            print('Density of ' + str(test['output']) + 'g/cm^3  for ' + test['input'] + '  not ' + str(result1) + 'g/cm^3.')

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

        # Testing keys
        if set(list(result1.keys())) != set(list(test['output'].keys())):
            print('Keys for find_stoichiometry not as expected.')
            print(test['reason'])

        # Test key order
        if list(result1.keys()) != list(test['output'].keys()):
            print('Keys and values not in correct order.')
            print(test['reason'])

        for check in list(result1.keys()):
            answer = test['output'][check]
            test_function = result1[check]

            # Checking that name input is correct in element class
            if answer.name != test_function.name:
                print('Name incorrect in element class.')
                print(test['reason'])

            # Checking that stoichiometry input is correct in element class
            if answer.stoichiometry != test_function.stoichiometry:
                print('Stoichiometry of element incorrect.')
                print(test['reason'])


    print('5. find_stociometry testing complete')

    # -------------------------------------------TEST 06-------------------------------------------------------------- #
    # -----------------------------------------Class: slab------------------------------------------------------------ #

