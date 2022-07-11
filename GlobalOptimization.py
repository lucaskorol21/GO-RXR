from scipy import optimize
from material_structure import *
import numpy as np
from data_structure import *
import matplotlib.pyplot as plt
import time
from tkinter import *
from material_model import *


def changeSampleParams(x, parameters, sample):
    """
    Purpose: Change the parameters for the input sample for x-values
    :param x: A list that contains the new parameters values
    :param parameters: A list that defines which parameter to change in the sample definition
    :param sample: A slab object class
    :return: The new sample
    """
    # Loop through each sample parameter
    for p in range(len(parameters)):
        params = parameters[p]

        if params[0] == "SCALING FACTOR":
            sample.scaling_factor = x[p]
        elif params[0] == "BACKGROUND SHIFT":
            sample.background_shift = x[p]
        elif params[0] == 'SCATTERING FACTOR':
            mode =params[1]
            element = params[2]
            dE = x[p]

            if mode == 'STRUCTURAL':
                FfEnergyShift(element, dE, opt=True)
            elif mode == 'MAGNETIC':
                FfmEnergyShift(element, dE, opt=True)
        else:
            layer = params[0]  # Determine which layer to change
            property = params[1]  # structural/polymorphous/magnetic
            if property == 'STRUCTURAL':
                mode = params[2]
                # Determine if user wants to use compound or element mode
                if mode == 'COMPOUND':
                    characteristic = params[3]  # thickness/density/roughness/linked roughness
                    for ele in list(sample.structure[layer].keys()):
                        if characteristic == 'THICKNESS':
                            sample.structure[layer][ele].thickness = x[p]
                        elif characteristic == 'DENSITY':
                            stoich = sample.structure[layer][ele].stoichiometry  # stoichiometry
                            molar_mass = sample.structure[layer][ele].molar_mass  # molar mass
                            sample.structure[layer][ele].density = x[p]*stoich/molar_mass  # set density
                        elif characteristic == 'ROUGHNESS':
                            sample.structure[layer][ele].roughness = x[p]
                        elif characteristic == 'LINKED ROUGHNESS':
                            sample.structure[layer][ele].linked_roughness = x[p]

                elif mode == 'ELEMENT':
                    element = params[3]  # determines the element to change
                    characteristic = params[4]  # thickness/density/roughness/linked roughness
                    if characteristic == 'THICKNESS':
                        sample.structure[layer][element].thickness = x[p]
                    elif characteristic == 'DENSITY':
                        sample.structure[layer][element].density = x[p]
                    elif characteristic == 'ROUGHNESS':
                        sample.structure[layer][element].roughness = x[p]
                    elif characteristic == 'LINKED ROUGHNESS':
                        sample.structure[layer][element].linked_roughness = x[p]


            elif property == 'POLYMORPHOUS':
                element = params[3]  # determines the element that contains the polymorph
                polymorph = params[4]  # determines the polymorph to change

                ratio = 1 - x[p]  # Assumes only two possible polymorphs for now and computes other polymorph ratio

                poly = np.where(sample.structure[layer][element].polymorph == polymorph)  # determines location of polymorph

                # sets poly_ratio value making sure sum equals to 1
                if poly == 0:
                    sample.structure[layer][element].poly_ratio[0] = x[p]
                    sample.structure[layer][element].poly_ratio[1] = ratio
                elif poly == 1:
                    sample.structure[layer][element].poly_ratio[1] = x[p]
                    sample.structure[layer][element].poly_ratio[0] = ratio

            elif property == 'MAGNETIC':
                element = params[3]  # determines the magnetic element to use

                # determines if magnetic element is polymorphous
                if len(params) == 3:
                    sample.structure[layer][element].mag_density[0] = x[p]  # non-polymorphous case
                else:
                    polymorph = params[4]  # polymorphous case
                    poly = np.where(sample.structure[layer][element].polymorph == polymorph)
                    sample.structure[layer][element].mag_density[poly] = x[p]

    return sample

def scanCompute(x, *args):

    chi2 = 0  # what we will use to determine some values

    sample = args[0]
    scans = args[1]
    data = args[2]
    sims = args[3]
    parameters = args[4]

    sample = changeSampleParams(x, parameters, sample)



    for scan in scans:
        scanType = scan[1]
        name = scan[2]
        if scanType == 'Reflectivity':
            myDataScan = data[name]
            myData = list(myDataScan)
            E = myDataScan.attrs['Energy']
            pol = myDataScan.attrs['Polarization']
            qz = np.array(myData[0])

            Rdat = np.log(np.array(myData[2]))
            qz, Rsim = sample.reflectivity(E, qz)
            Rsim = np.log10(Rsim[pol])*sample.scaling_factor + np.ones(len(Rsim[pol]))*sample.background_shift
            chi2 = chi2 + sum((Rdat-Rsim)**2/abs(Rsim))




        elif scanType == 'Energy':
            myDataScan = data[name]
            myData = list(myDataScan)
            Theta = myDataScan.attrs['Angle']
            E = np.array(myData[3])
            pol = myDataScan.attrs['Polarization']

            Rdat = np.log(np.array(myData[2]))
            Rsim = sample.energy_scan(Theta, E)
            Rsim = np.log(Rsim[pol])

            chi2 = chi2 + sum((Rdat-Rsim)**2/abs(Rsim))


    return chi2


def global_optimization(fname, scan, parameters, bounds, algorithm = 'differential_evolution'):
    sample = ReadSampleHDF5(fname)  # import the sample information

    f, data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data


    # makes sure that scan is a list
    if type(scan) != list or type(scan) != np.ndarray:
        scan = [scan]

    scans = data_info[tuple(scan)]

    params = [sample, scans, data, sims, parameters]  # required format for function scanCompute

    if algorithm == 'differential_evolution':
        # This line will be used to select and use different global optimization algorithms
        ret = optimize.differential_evolution(scanCompute, bounds, args=params, strategy='currenttobest1exp', maxiter=10, tol=0.00001, disp=True)
        x = ret.x
        fun = ret.fun
    elif algorithm == 'shgo':
        ret = optimize.shgo(scanCompute, bounds, args=tuple(params), n=64, iters=3)
        x = ret.x
        fun = ret.fun
    elif algorithm == 'dual_annealing':
        ret = optimize.dual_annealing(scanCompute, bounds, args=params, maxiter=300)
        x = ret.x
        fun = ret.fun
    else:
        raise NameError('Parameter algorithm not set properly')

    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)
    sample = changeSampleParams(x, parameters, sample)


    #WriteSampleHDF5(fname, sample)
    for scan in scans:
        scanType = scan[1]
        name = scan[2]
        if scanType == 'Reflectivity':
            myDataScan = data[name]
            myData = list(myDataScan)
            E = myDataScan.attrs['Energy']
            pol = myDataScan.attrs['Polarization']
            qz = np.array(myData[0])

            Rdat = np.array(myData[2])
            qz, Rsim = sample.reflectivity(E, qz)
            Rsim = Rsim[pol]

            plt.figure()
            plt.suptitle('Fit: ' + name)
            plt.plot(qz,np.log10(Rdat))
            plt.plot(qz, np.log10(Rsim))
            plt.legend(['Data','Simulation'])
            plt.show()

        elif scanType == 'Energy':
            myDataScan = data[name]
            myData = list(myDataScan)
            Theta = myDataScan.attrs['Angle']
            E = np.array(myData[3])
            pol = myDataScan.attrs['Polarization']

            Rdat = np.array(myData[2])
            Rsim = sample.energy_scan(Theta, E)
            Rsim = Rsim[pol]



    f.close()
    return x, fun

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

    # ----------------------------------------- Parameter Selection -------------------------------------------------- #
    print()
    print('PARAMETER SELECTION')
    print()

    layer_formula = list()
    idx = 0
    for temp_layer in sample.structure:
        P = list()
        M = list()
        formula = ''
        for key in list(temp_layer.keys()):
            stoich = temp_layer[key].stoichiometry
            if stoich == 1:
                formula = formula + key
            else:
                formula = formula + key + str(stoich)
            if len(temp_layer[key].polymorph) > 0:
                P.append(key)
            if len(temp_layer[key].mag_density) > 0:
                M.append(key)

        layer_formula.append([idx, formula, P, M])
        idx = idx + 1

        # Sets up and prints scan table
    header = ['Layer', 'Formula', 'Polymorphs', 'Magnetic']
    tab = PrettyTable(header)
    tab.add_rows(layer_formula)
    print(tab)
    print()
    number_layers = len(sample.structure)

    cont = 'Y'
    while cont.upper() == 'Y':
        layer = input('Select layer you would like to optimize (0-'+str(number_layers-1)+"): ")
        if layer.upper() == 'EXIT':
            quit()
        while layer.upper() == 'SHOW':
            print(tab)
            print()
            layer = input('Select layer you would like to optimize (0-' + str(number_layers - 1) + "): ")
            if layer.upper() == 'EXIT':
                quit()

        while int(layer) < 0 or int(layer) > number_layers -1:
            layer = input('Select layer you would like to optimize (0-' + str(number_layers - 1) + "): ")
            if layer.upper() == 'EXIT':
                quit()
            while layer.upper() == 'SHOW':
                print(tab)
                print()
                layer = input('Select layer you would like to optimize (0-' + str(number_layers - 1) + "): ")
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

                    print()
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

                        print()
                        sample_params.append([layer,prop.upper(),mode.upper(),element, characteristic.upper()])
                        if len(char_list) > 0:
                            again = input('Would you liked to select another characteristic for '+ element+" (y/n): ")
                        else:
                            again = 'N'
                    val = input('Would you like to select another element (y/n)?')
        elif prop.upper() == "POLYMORPHOUS":
            b= []
            for ele in elements:
                if len(sample.structure[layer][ele].polymorph) > 0:
                    b.append(ele)
            poly_cont = 'Y'
            while poly_cont.upper() == 'Y' and len(b) > 0:
                poly_ele = input('Select the polymorph '+ str(b) +' you would liked to vary: ')
                if poly_ele.upper() == 'EXIT':
                    quit()

                num_poly = len(b)

                while poly_ele not in b:
                    poly_ele = input('Please select a polymorph ' + str(b) + ' : ')
                    if poly_ele.upper() == 'EXIT':
                        quit()


                polymorph = sample.structure[layer][poly_ele].polymorph
                whichPoly = input('Select polymorph (' + str(polymorph) + ') for which density you would like to vary?')
                if whichPoly.upper() == 'EXIT':
                    quit()
                while whichPoly not in polymorph:
                    whichPoly = input('Select polymorph (' + str(polymorph) + ') for which density you would like to vary?')
                    if whichPoly.upper() == 'EXIT':
                        quit()

                lw = input("Select the lower limit of the polymorph ratio (0-1): ")
                if lw.upper() == 'EXIT':
                    quit()
                up = input("Select the upper limit of the polymorph ratio (0-1): ")
                if up.upper() == 'EXIT':
                    quit()

                while(float(lw)>float(up) or float(lw)<0 or float(up) > 1):
                    if float(lw)>float(up):
                        lw = input("Make sure that your lower bound is smaller than the upper bound. Please select a new lower limit: ")
                        if lw.upper() == 'EXIT':
                            quit()
                        up = input("Please select a new upper bound: ")
                        if up.upper() == 'EXIT':
                            quit()
                    if float(lw) < 0:
                        lw = input("Please select a lower limit between 0 and 1: ")
                        if lw.upper() == 'EXIT':
                            quit()
                    if float(up) > 1:
                        up = input("Please select an upper limit between 0 and 1: ")
                        if up.upper() == 'EXIT':
                            quit()

                lowerbound.append(lw)
                upperbound.append(up)

                sample_params.append([layer, prop.upper(), poly_ele, whichPoly])
                # As of right now we will assume that we can have a maximum of 2 polymorphs
                print()
                if len(b) != 0:
                    b.remove(poly_ele)
                    poly_cont = input('Would you like to vary another polymorph in the same layer (y/n): ')
                    if poly_cont.upper() == 'EXIT':
                        quit()

        elif prop.upper() == 'MAGNETIC':
            my_mag = list()
            for ele in elements:
                if len(sample.structure[layer][ele].mag_density) > 0:
                    my_mag.append(ele)

            mag_ele = input('Select magnetic element you would like to vary (' + str(my_mag) + ': ')
            if mag_ele.upper() == 'EXIT':
                quit()
            while mag_ele not in my_mag:
                mag_ele = input('Please select an element element in (' + str(my_mag) + ': ')
                if mag_ele.upper() == 'EXIT':
                    quit()

            mag_poly = list(sample.structure[layer][mag_ele].polymorph)
            mag_cont = 'Y'
            while(mag_cont.upper() == 'Y' and len(my_mag)>0):
                if len(mag_poly) > 0:
                    mag_poly_cont = 'Y'
                    while(mag_poly_cont.upper() == 'Y' and len(mag_poly)>0):
                        whichMag = input('Select which polymorph ('+str(mag_poly)+') for which you would like to vary the magnetic density: ')
                        if whichMag.upper() == 'EXIT':
                            quit()
                        while whichMag not in mag_poly:
                            whichMag = input('Select which polymorph ('+str(mag_poly)+') for which you would like to vary the magnetic density: ')
                            if whichMag.upper() == 'EXIT':
                                quit()

                        lw = input('Enter the lower bound of the magnetic density of order mol/cm^3: ')
                        if lw.upper() == 'EXIT':
                            quit()
                        up = input('Enter the lower bound of the magnetic density of order mol.cm^3: ')
                        if up.upper() == 'EXIT':
                            quit()
                        while float(lw) > float(up):
                            lw = input('Make sure lower bound is smaller than upper bound. Enter lower bound again: ')
                            if lw.upper() == 'EXIT':
                                quit()
                            up = input('Enter upper bound: ')
                            if up.upper() == 'EXIT':
                                quit()
                        lowerbound.append(lw)
                        upperbound.append(up)

                        sample_params.append([layer, prop.upper(), mag_ele, whichMag])
                        print()
                        if len(mag_poly) != 0:
                            mag_poly_cont = input('Would you like to vary another polymorph magnetic density (y/n)?')
                            if mag_poly_cont.upper() == 'EXIT':
                                quit()
                            mag_poly.remove(whichMag)

                else:
                    lw = input('Enter the lower bound of the magnetic density of order mol/cm^3: ')
                    if lw.upper() == 'EXIT':
                        quit()
                    up = input('Enter the lower bound of the magnetic density of order mol.cm^3: ')
                    if up.upper() == 'EXIT':
                        quit()
                    while float(lw) > float(up):
                        lw = input('Make sure lower bound is smaller than upper bound. Enter lower bound again: ')
                        if lw.upper() == 'EXIT':
                            quit()
                        up = input('Enter upper bound: ')
                        if up.upper() == 'EXIT':
                            quit()
                    lowerbound.append(lw)
                    upperbound.append(up)
                    print()
                    sample_params.append([layer, prop.upper(), mag_ele])
                if len(my_mag) != 0:
                    my_mag.remove(mag_ele)
                    mag_cont = input('Would you like to select another magnetic element to vary (y/n)?')
                    if mag_cont.upper() == 'EXIT':
                        quit()

        cont = input('Would you liked to select another layer (y/n): ')


    # printing the chosen elements
    my_params = list()

    for sp in range(len(sample_params)):
        temp_list = list()
        params = sample_params[sp]
        lw = lowerbound[sp]
        up = upperbound[sp]

        temp_list.append(params[0])  # Layer info
        temp_list.append(params[1])  # Property

        if params[1] == 'STRUCTURAL':
            mode = params[2]  # mode
            if mode == 'COMPOUND':
                formula = layer_formula[sp][1]  # formula
                characteristic = params[3]
                temp_list.append(formula)
                temp_list.append('N/A')
                temp_list.append(characteristic)
            else:
                element = params[3]
                characteristic = params[4]
                temp_list.append(element)
                temp_list.append('N/A')
                temp_list.append(characteristic)
        elif params[1] == 'POLYMORPHOUS':
            temp_list.append(params[2])  # poly element
            temp_list.append(params[3])  # which poly
            temp_list.append('DENSITY RATIO')
        elif params[1] == 'MAGNETIC':
            if len(params) == 3:
                temp_list.append(params[2])  # element
                temp_list.append('N/A')
                temp_list.append('MAGNETIC DENSITY')
            else:
                temp_list.append(params[2])  # poly element
                temp_list.append(params[3])  # polymorph
                temp_list.append('MAGNETIC DENSITY')

        temp_list.append(lw)
        temp_list.append(up)
        my_params.append(temp_list)

    print()
    print('The list of scans you requested are:' + str(scans))
    print()
    print(my_params)
    header = ['Layer', 'Property', 'Element(s)', 'Polymorph', 'Characteristic', 'Upper Bound', 'Lower Bound']
    Ntab = PrettyTable(header)
    Ntab.add_rows(my_params)
    print(Ntab)

    return data_info[scans], data, sims, sample_params, lowerbound, upperbound


if __name__ == "__main__":
    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density=[0.027904, 0.027904, 0.083712], roughness=[7.58207, False, 5.77093])
    sample.addlayer(1, 'SrTiO3', 6, density=[0, 0.027904, 0], roughness=[7.58207, 4.03102, 5.77093])

    sample.addlayer(2, 'LaMnO3', 4, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2],
                    linked_roughness=[False, 0.5, False])
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 17.8, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 9, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 2.5, density=[0.025, 0.024, 0.05], roughness=[0.25, 0.25, 2])
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 4.5, density=[0.025, 0.042, 0.04], roughness=[0.25, 0.25, 2])
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.0053, 0], ['Co', 'Ni'])

    sample.addlayer(7, 'CCO', 11.1, density=[0.05, 0.05, 0.01], roughness=2, linked_roughness=[3, 1.5, False])

    fname = 'Pim10uc.h5'

    #sample.plot_density_profile(1)
    #plt.show()
    #WriteSampleHDF5(fname, sample)
    #print(ReadDataHDF5(fname))


    parameters = [[1, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  [2, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  [3, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  [4, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  ['SCATTERING FACTOR', 'STRUCTURAL', 'La']]


    lw = [3.5,3.5,17.3,8.5,-0.5]
    up = [6.5,6.5,19.8,11.5,0.5]
    bounds = list(zip(lw, up))
    scans = [0,1,2,3,4,5,6]


    start = time.time()
    x, fun = global_optimization(fname, scans, parameters, bounds, algorithm='differential_evolution')
    end = time.time()
    print(end-start)

    #x_expected = np.array([5,5, 18.8437, 10, 3, 4, 10.1373])

    #Var = sum((x-x_expected)**2)/len(x_expected)
    #print(Var)

