
from material_structure import *
import matplotlib
import numpy as np

if __name__ == "__main__":


    # ------------------------------------------TEST 01--------------------------------------------------------------- #
    # -------------------------------------find_stoichiometry--------------------------------------------------------- #
    test_01 = [{'input': 'FeNeO3', 'output': {'Fe': element('Fe', 1), 'Ne': element('Ne',1), 'O': element('O', 3)}, 'reason': 'FeNeO3 Incorrect'},
               {'input': 'La2222Mn7784848O', 'output': {'La': element('La',2222), 'Mn': element('Mn', 7784848), 'O': element('O',1)}, 'reason': 'La2222Mn7784848O incorrect'},
               {'input': 'FeMn2', 'output': {'Fe':element('Fe',1), 'Mn': element('Mn',2)}, 'reason': 'FeMn2 incorrect'},
               {'input': 'CoZn', 'output': {'Co':element('Co',1), 'Zn': element('Zn',1)}, 'reason': 'CoZn incorrect'},
               {'input': 'BaLa6K2', 'output': {'Ba':element('Ba',1), 'La': element('La',6), 'K': element('K',2)}, 'reason': 'BLa6K2 incorrect'},
               {'input': 'CuK7MnO', 'output': {'Cu':element('Cu',1), 'K': element('K',7), 'Mn': element('Mn',1), 'O': element('O',1)}, 'reason': 'CuK7MnO incorrect'},
               {'input': 'FOs2', 'output': {'F': element('F', 1), 'Os': element('Os', 2)}, 'reason': 'FOs2 incorrect'},
               {'input': 'CZn', 'output': {'C': element('C', 1), 'Zn': element('Zn', 1)}, 'reason': 'CZn incorrect'},
               {'input': 'BLa6N2', 'output': {'B': element('B', 1), 'La': element('La', 6), 'N': element('N', 2)}, 'reason': 'BLa6N2 incorrect'},
               {'input': 'VK7MnS', 'output': {'V': element('V', 1), 'K': element('K', 7), 'Mn': element('Mn', 1), 'S': element('S', 1)}, 'reason': 'VK7MnS incorrect'}]


    for test in test_01:
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


    print('1. find_stociometry test succesful')