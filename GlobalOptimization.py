from scipy import optimize
from material_structure import *
import numpy as np
from data_structure import *


# Begin with creating a function that I can use global optimization on!
def global_optimization():
    # The idea for this would be to allow the user to select the scans they would like to
    # use in the global optimization algorithm

    # We will then ask the user which parameters they would like to vary, along with the range

    # We may need a very complex function to perform the optimization...


    print('good try')

def brute_force_optimization():
    # This function will use the brute optimization technique on our system
    print('hello')

def selectOptimize(fname):
    sample = ReadSampleHDF5(fname)  # import the sample information
    data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    # Sets up and prints scan table
    header = ['#', 'Scan Type', 'Scan']
    tab = PrettyTable(header)
    tab.add_rows(data_info)
    print(tab)

    scans = list()
    sample_params = list()
    sample_values = list()
    upperbound = list()
    lowerbound = list()

    # ------------------------------------------ Scan Selection ------------------------------------------------------ #
    print('SCAN SELECTION')
    print()

    scan = input("Select a scan for global optimization ("+data_info[0,0]+"-"+data_info[-1,0] + "): ")
    scans.append(int(scan))
    print()

    # Requires user to select a new scan
    while scan not in data_info[:,0]:
        scan = input("Please select a scan from " + data_info[0,0]+" to "+data_info[-1,0])

    print()
    cont = input('Would you like to select another scan (y/n)?')

    while(cont.upper() == 'Y'):
        scan = input('Select another scan: ')

        # Checks to make sure there are no repeats
        while scan in scans:
            scan = input("Scan "+ scan+" has already been selected. Please select a different scan or typ \'n\' if you no longer want to enter a new scan: ")
            if scan.upper() == 'N':
                break
        if scan.upper() == 'N':
            break

        # Checks to make sure scan exists
        while scan not in data_info[:,0]:
            scan = input("Please select a scan from " + data_info[0,0]+" to "+data_info[-1,0]+ ". If you no longer want to select another scan please enter \'n\': ")

            if scan.upper() == 'N':
                break
        if scan.upper() == 'N':
            break
        else:
            scans.append(scan)
            print()
            cont = input('Would you like to select another scan (y/n)?')

    # ----------------------------------------- Parameter Selection -------------------------------------------------- #
    print()
    print('PARAMETER SELECTION')
    print()
    number_layers = len(sample.structure)
    layer = input('Select layer you would like to optimize (0-'+str(number_layers-1)+"): ")
    while int(layer) < 0 or int(layer)>number_layers-1:
        layer = input('Please select a layer between 0 and ' + str(number_layers-1) + " : ")

    elements = list(sample.structure[int(layer)].keys())

    element = input('Select element you would like to optimize ' + str(elements)+ " : ")
    while element not in elements:
        element = input('Please select an element from the list ' + str(elements) + ' : ')


    s = 'structure'
    ele = sample.structure[int(layer)][element]

    if len(ele.polymorph) != 0:
        s = s + '/polymorph'

    if len(ele.mag_density) != 0:
        s = s + '/magnetic'

    whatKind = input('Select property you would like to optimize ('+s+') : ')
    while whatKind.upper() != 'STRUCTURE' and whatKind.upper() != 'POLYMORPH' and whatKind.upper() != 'MAGNETIC':
        whatKind = input('Please enter property as ' + s + ' : ')

    if whatKind.upper() == 'STRUCTURE':
        print('hello')
    elif whatKind.upper() == 'POLYMORPH':
        print(1)
    elif whatKind.upper() == 'MAGNETIC':
        print()




    my_params = list()
    my_params.append(int(input('Please select layer you would like to optimize: ')))
    my_params.append(input('Select element you want to optimize: '))
    my_params.append(input('Which parameters would you like to optimize (struct/poly/mag)?'))
    my_params.append(input)


if __name__ == "__main__":
    fname = 'Pim10uc.h5'
    selectOptimize(fname)