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


    sample_params = list()
    upperbound = list()
    lowerbound = list()

    # ------------------------------------------ Scan Selection ------------------------------------------------------ #
    print('SCAN SELECTION')
    print()
    scans = input('Input scans from '+str(data_info[0,0])+' to '+ str(data_info[-1,0]) + ' for global optimization. Separate each scan using a space: ')
    if scans.upper() == 'EXIT':
        quit()
    scans = scans.split()
    scans = list(dict.fromkeys(scans))  # removes any duplicates
    for scan in scans:
        if scan not in data_info[:,0]:
            val = input('This scan number ' + scan + ' does not exist. Would you like to select another scan (y/n)?')
            if val.upper() == 'EXIT':
                quit()
            scans.remove(scan)
            if val.upper() == 'Y':
                scan = input('Input the new scan number: ')
                if scan.upper() == 'EXIT':
                    quit()
                while scan not in data_info[:,0]:
                    scan = input('Please input a scan that lies on the interval '+str(data_info[0,0])+' to '+ str(data_info[-1,0])+' : ')
                    if scan.upper() == 'EXIT':
                        quit()
                scans.append(scan)
    scans = [int(scan) for scan in scans]
    print(scans)



    # ----------------------------------------- Parameter Selection -------------------------------------------------- #
    print()
    print('PARAMETER SELECTION')
    print()

    number_layers = len(sample.structure)

    cont = 'y'
    while cont.upper() == 'Y':
        layer = input('Select layer you would like to optimize (0-'+str(number_layers-1)+"): ")
        if layer.upper() == 'EXIT':
            quit()

        poly = list()
        mag = list()

        structural = True
        polymorphous = False
        magnetic = False

        layer = int(layer)
        elements = list(sample.structure[layer].keys())
        for ele in elements:
            if len(sample.structure[layer][ele].polymorph) > 0:
                poly.append(ele)
                polymorphous = True
            if len(sample.structure[layer][ele].mag_density) > 0:
                mag.append(ele)
                magnetic = True

        s = 'structural'
        if polymorphous:
            s = s + "/polymorphous"
        if magnetic:
            s = s + "/magnetic"

        prop = input('Select the property you would like to vary ('+ s + ') : ')
        if prop.upper() == 'EXIT':
            quit()
        while prop.upper() != 'STRUCTURAL' and prop.upper() != 'POLYMORPHOUS' and prop.upper() != 'MAGNETIC':
            prop = input('Please select one of the properties (' + s + ') : ')
            if prop.upper() == 'EXIT':
                quit()

        if prop.upper() == 'STRUCTURAL':
            mode = input('Select mode (element/compound): ')
            if mode.upper() == 'EXIT':
                quit()
            while mode.upper() != 'ELEMENT' and mode.upper() != 'COMPOUND':
                mode = input('Please select (element/compound): ')
                if mode.upper() == 'EXIT':
                    quit()

            if mode.upper() == 'COMPOUND':
                char_list = ['THICKNESS','DENSITY', 'ROUGHNESS', 'LINKED ROUGHNESS']
                val = 'y'
                while val.upper() == 'Y' and len(char_list) > 0:
                    characteristic = input('Select characteristic (' + '/'.join([char.lower() for char in char_list]) + '): ')
                    if characteristic.upper() == 'EXIT':
                        quit()
                    while characteristic.upper() != 'THICKNESS' and characteristic.upper() != 'DENSITY' and characteristic.upper() != 'ROUGHNESS' and characteristic.upper() != 'LINKED ROUGHNESS':
                        characteristic = input('Select characteristic (' + '/'.join([char.lower() for char in char_list]) + '): ')
                        if characteristic.upper() == 'EXIT':
                            quit()

                    if characteristic.upper() == 'THICKNESS' and 'THICKNESS' in char_list:
                        char_list.remove('THICKNESS')
                        lw = 1
                        up = 0
                        while float(lw) > float(up):
                            lw = input('Select lower bound in units of Angstrom: ')
                            if lw.upper() == 'EXIT':
                                quit()
                            up = input('Select upper bound in units of Angstrom: ')
                            if up.upper() == 'EXIT':
                                quit()

                        lowerbound.append(lw)
                        upperbound.append(up)
                    if characteristic.upper() == 'DENSITY' and 'DENSITY' in char_list:
                        char_list.remove('DENSITY')
                        lw = 1
                        up = 0
                        while float(lw) > float(up):
                            lw = input('Select lower bound in units of g/cm^3: ')
                            if lw.upper() == 'EXIT':
                                quit()
                            up = input('Select upper bound in units of g/cm^3: ')
                            if up.upper() == 'EXIT':
                                quit()
                        lowerbound.append(lw)
                        upperbound.append(up)
                    if characteristic.upper() == 'ROUGHNESS' and 'ROUGHNESS' in char_list:
                        char_list.remove('ROUGHNESS')
                        lw = 1
                        up = 0
                        while float(lw) > float(up):
                            lw = input('Select lower bound in units of Angstrom: ')
                            if lw.upper() == 'EXIT':
                                quit()
                            up = input('Select upper bound in units of Angstrom: ')
                            if up.upper() == 'EXIT':
                                quit()
                        lowerbound.append(lw)
                        upperbound.append(up)
                    if characteristic.upper() == 'LINKED ROUGHNESS' and 'LINKED ROUGHNESS' in char_list:
                        char_list.remove('LINKED ROUGHNESS')
                        lw = 1
                        up = 0
                        while float(lw) > float(up):
                            lw = input('Select lower bound in units of Angstrom: ')
                            if lw.upper() == 'EXIT':
                                quit()
                            up = input('Select upper bound in units of Angstrom: ')
                            if up.upper() == 'EXIT':
                                quit()
                        lowerbound.append(lw)
                        upperbound.append(up)

                    sample_params.append([layer,prop.upper(),mode.upper(),characteristic.upper()])

                    if len(char_list) > 0:
                        val = input('Would you like to select another characteristic (y/n)?')
                    else:
                        val = 'N'


            elif mode.upper() == 'ELEMENT':
                char_list = ['THICKNESS','DENSITY', 'ROUGHNESS', 'LINKED ROUGHNESS']
                val = 'y'
                while val.upper() == 'Y':
                    element = input('Select element (' + str(elements) + ") : ")
                    while element not in elements:
                        element = input('Please select element in list (' + str(elements) + ") : ")

                    again = 'y'
                    while again.upper() == 'Y' and len(char_list)>0:
                        characteristic = input('Select characteristic (' + '/'.join([char.lower() for char in char_list]) + '): ')
                        while characteristic.upper() != 'THICKNESS' and characteristic.upper() != 'DENSITY' and characteristic.upper() != 'ROUGHNESS' and characteristic.upper() != 'LINKED ROUGHNESS':
                            characteristic = input('Select characteristic (' + '/'.join([char.lower() for char in char_list]) + '): ')

                        if characteristic.upper() == 'THICKNESS' and 'THICKNESS' in char_list:
                            char_list.remove('THICKNESS')
                            lw = 1
                            up = 0
                            while float(lw) > float(up):
                                lw = input('Select lower bound in units of Angstrom: ')
                                if lw.upper() == 'EXIT':
                                    quit()
                                up = input('Select upper bound in units of Angstrom: ')
                                if up.upper() == 'EXIT':
                                    quit()
                            lowerbound.append(lw)
                            upperbound.append(up)
                        if characteristic.upper() == 'DENSITY' and 'DENSITY' in char_list:
                            char_list.remove('DENSITY')
                            lw = 1
                            up = 0
                            while float(lw) > float(up):
                                lw = input('Select lower bound in units of mol/cm^3: ')
                                if lw.upper() == 'EXIT':
                                    quit()
                                up = input('Select upper bound in units of mol/cm^3: ')
                                if up.upper() == 'EXIT':
                                    quit()
                            lowerbound.append(lw)
                            upperbound.append(up)
                        if characteristic.upper() == 'ROUGHNESS' and 'ROUGHNESS' in char_list:
                            char_list.remove('ROUGHNESS')
                            lw = 1
                            up = 0
                            while float(lw) > float(up):
                                lw = input('Select lower bound in units of Angstrom: ')
                                if lw.upper() == 'EXIT':
                                    quit()
                                up = input('Select upper bound in units of Angstrom: ')
                                if up.upper() == 'EXIT':
                                    quit()
                            lowerbound.append(lw)
                            upperbound.append(up)
                        if characteristic.upper() == 'LINKED ROUGHNESS' and 'LINKED ROUGHNESS' in char_list:
                            char_list.remove('LINKED ROUGHNESS')
                            lw = 1
                            up = 0
                            while float(lw) > float(up):
                                lw = input('Select lower bound in units of Angstrom: ')
                                if lw.upper() == 'EXIT':
                                    quit()
                                up = input('Select upper bound in units of Angstrom: ')
                                if up.upper() == 'EXIT':
                                    quit()
                            lowerbound.append(lw)
                            upperbound.append(up)

                        sample_params.append([layer,prop.upper(),mode.upper(),element, characteristic.upper()])
                        if len(char_list) > 0:
                            again = input('Would you liked to select another characteristic for '+ element+" (y/n): ")
                        else:
                            again = 'N'
                    val = input('Would you like to select another element (y/n)?')

        cont = input('Would you liked to select another layer (y/n): ')






    return


if __name__ == "__main__":
    fname = 'Pim10uc.h5'
    selectOptimize(fname)