"""
Library: data_structure
Version: 0.3
Author: Lucas Korol
Institution: University of Saskatchewan
Last Updated: March 28nd, 2023
Python: version 3.7

Purpose: This python file contains all the reading and writing functions

---------------------------------------------------------------------------------------------------------------------------------------
Imported Libraries

numpy (version 1.21.4) - used for array manipulation

material_model (version 0.1) - Part of the 'name' software package and is used to retrieve the form factors and
                               calculate the optical constants

material_structure (version 0.1) - Part of the 'name' software package and is used to retrieve to calculate the reflectivity spectra

h5py (version 2.90) - Used to save files in hdf5 form

--------------------------------------------------------------------------------------------------------------------------------------
Note: Any additional parameters that need to be saved needs to be included into all functions.
In particular, if additional global optimization algorithms are included this will need to be added all reading and
writing functions (saveAsFileHDF5, saveFileHDF5, newFileHDF5, WriteDataHDF5, ReadAlgorithmHDF5). View saveAsFile for
explanation of how to include a new global optimization algorithm.

Warning: Once you've made a change to the HDF5 file format and you have made changes to the reading functions then you will no
longer be able to read in older files. Be careful when working with the older files while testing as this
may corrupt the files. I've done this a few times and lost precious work that I need to redo.

Version: I've included a version attribute for the file for the HDF5 file types and automatically save the
GO-RXR version. There is not current use for this, but there are numerous updates that will include different
global optimization algorithms, sample definitions, and fitting parameters. This means that the loading and saving
functions will need to change accordingly. However, there is a problem with the loading functions as a a workspace file
that was saved using an older version of GO-RXR will not be able to load properly. So the idea behind the
version attribute was to be able to identify which version this file was last saved so that the appropriate loading
sequence can be done.
"""

import matplotlib.pyplot as plt
import os

import numpy as np

from UTILS.material_structure import *
from UTILS.material_model import *
from time import *
import ast
import h5py


def getTitleInfo(title):
    """
    Purpose: Retrieves important information in the scan title from ReMagX
    :param title: title/label of the scan
    :return:
    """

    title = title.split('_')  # split title into an array

    angle = None
    energy = ''
    scan_number = title[0]
    polarization = 'S'

    # Determines the scan type
    if title[1] == 'A':
        energy = title[2]
        scanType = 'Reflectivity'
        if title[3] == 'S':
            polarization = 'S'
        elif title[3] == 'P':
            polarization = 'P'
        elif title[3] == 'L':
            polarization = 'LC'
        elif title[3] == 'R':
            polarization = 'RC'
        elif 'R' in title[3] and 'L' in title[3]:
            polarization = 'AC'
        elif 'S' in title[3] and 'P' in title[3]:
            polarization = 'AL'
    else:
        energy = title[1]
        scanType = 'Energy'
        angle = title[2].replace('Th','')
        if title[3] == 'S':
            polarization = 'S'
        elif title[3] == 'P':
            polarization = 'P'
        elif title[3] == 'L':
            polarization = 'LC'
        elif title[3] == 'R':
            polarization = 'RC'
        elif 'R' in title[3] and 'L' in title[3]:
            polarization = 'AC'
        elif 'S' in title[3] and 'P' in title[3]:
            polarization = 'AL'
    return scan_number, scanType, energy, polarization, angle

def Read_ReMagX(fname):
    """
    Purpose: Read ReMagX files and convert to a data and simulation dictionary
    :param fname: ReMagX filename
    :return: data_dict and sim_dict
    """
    data_dict = dict()
    data_info = []
    if fname.endswith('.all'):  # checks to make sure selected file is a ReMagX file type
        # good to go
        is_data = False
        new_data = True

        current_scan = 0  # keeps track of which scan we are on!

        # keeps track of the data
        qz_list = list()
        theta_list = list()
        R_list = list()
        E_list = list()
        scanType = 'Reflectivity'
        title = ''
        name = ''

        with open(fname, 'r') as f:
            for line in f.readlines():
                line = line.split()

                # determines if we can start reading in the data
                if '#' in line and 'measurement' in line and 'data' in line:
                    is_data = True

                if '#' in line and 'dataset' in line and str(current_scan+1) in line:
                    is_data = False

                if len(line) == 3:
                    if line[1] == '=':  # makes sure we are only evaluating what we want to be

                        if line[0] == 'datasetnumber':
                            current_scan = int(line[2])
                            new_data = True

                        # new data
                        if new_data:
                            if current_scan > 1:  # where we can set the data
                                if len(E_list) != 0:
                                    data = [qz_list, theta_list, R_list, E_list]
                                    data_info.append([current_scan, 'Energy', name])
                                else:
                                    data = [qz_list, theta_list, R_list]
                                    data_info.append([current_scan, 'Reflectivity', name])
                                data_dict[name]['Data'] = np.array(data)

                            qz_list = list()
                            theta_list = list()
                            R_list = list()
                            E_list = list()
                            new_data = False


                        if line[0] == 'datasettitle':
                            scan_number, scanType, energy, pol, angle = getTitleInfo(line[2])
                            name = ''
                            if scanType == 'Reflectivity':
                                name = scan_number +'_' + energy + '_' + pol
                                if pol == 'AL' or pol == 'AC':
                                    name = name + '_Asymm'
                                data_dict[name] = dict()
                                data_dict[name]['Polarization'] = pol
                                data_dict[name]['DatasetNumber'] = current_scan
                                data_dict[name]['Background Shift'] = 0
                                data_dict[name]['Scaling Factor'] = 1


                            elif scanType == 'Energy':
                                name = scan_number + '_' + energy + '_Th'+ angle  + '_' + pol
                                if pol == 'AL' or pol == 'AC':
                                    name = name + '_Asymm'

                                data_dict[name] = dict()
                                data_dict[name]['Polarization'] = pol
                                data_dict[name]['Angle'] = float(angle)
                                data_dict[name]['DatasetNumber'] = current_scan
                                data_dict[name]['Background Shift'] = 0
                                data_dict[name]['Scaling Factor'] = 1


                        elif line[0] == 'datasetenergy':
                            data_dict[name]['Energy'] = float(line[2])
                            E = float(line[2])
                        elif line[0] == 'datasetpoints':
                            data_dict[name]['DataPoints'] = int(line[2])
                            # do something else
                        elif line[0] == 'dataset_qz':
                            qz = float(line[2])
                            try:
                                theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi
                            except RuntimeWarning:
                                theta=0


                            qz_list.append(qz)
                            theta_list.append(theta)
                        elif line[0] == 'dataset_A':
                            R_list.append(float(line[2]))
                        elif line[0] == 'dataset_R0':
                            R_list.append(float(line[2]))
                        elif line[0] == 'dataset_eng':
                            E_list.append(float(line[2]))

    if len(E_list) != 0:
        data = [qz_list, theta_list, R_list, E_list]
        data_info.append([current_scan, 'Energy', name])
    else:
        data = [qz_list, theta_list, R_list]
        data_info.append([current_scan, 'Reflectivity', name])
    data_dict[name]['Data'] = np.array(data)
    f.close()


    return data_info, data_dict

def evaluate_list(string):
    """
    Purpose: This is used to evaluate scan boundaries
    :param string: Scan boundaries in string form
    :return: Boundary in list format
    """

    # check to make sure it is a list
    if type(string) is not np.ndarray and type(string) is not list:
        if string[0] != '[' and string[-1] != ']':
            raise TypeError('String is not a list.')

        final_array = []
        if string != '[]':
            string = string.strip('[')  # remove brackets
            string = string.strip(']')  # remove brackets

            array1 = string.split(',')  # get all the components of the list
            array2 = string.split()
            #array2.remove(',')

            my_array = []
            if len(array1) > len(array2):  # string has commas
                my_array = array1
            elif len(array1)< len(array2):  # string has spaces
                my_array = array2
            elif len(array1) == len(array2):  # string has spaces and commas
                my_array = array1

            if ',' in my_array:
                my_array = list(filter(lambda item: item != ',', my_array))

            final_array = [float(s.replace(',','')) for s in my_array]  # transforms in back into a list of strings
    else:
        final_array = string

    return final_array


def find_parameter_bound(string):
    """
    Purpose: Convert parameter boundaries from string type to list type
    :param string: Parameter boundary as string type
    :return: Parameter boundary as list type
    """
    my_list = []

    my_string = ''
    count_open = 0
    count_closed = 0

    for s in string:
        if s != '\'':
            my_string = my_string + s  # building string
        if s == '[':
            count_open = count_open + 1
        if s == ']':
            count_closed = count_closed + 1
        if count_open == 2 and count_closed == 2:
            my_idx = my_string.find('[')
            my_string = my_string[my_idx:]
            my_list.append(my_string)
            count_open = 0
            count_closed = 0
            my_string = ''

    return my_list

def find_closing_bracket(string, open_index):
    stack = []  # Use a stack to keep track of nested brackets
    for i in range(open_index, len(string)):
        if string[i] == '[':
            stack.append('[')
        elif string[i] == ']':
            if not stack:
                return -1  # No matching opening bracket found
            stack.pop()
            if not stack:
                return i  # Found the matching closing bracket

    return -1  # No matching closing bracket found

def find_each_bound(string):
    #string = string[1:-1]
    my_bounds = []

    while string != '':
        if string[0] == ',' or string[0] == ' ':
            string = string[1:]
        else:
            idx = find_closing_bracket(string,0)
            my_bounds.append(string[0:idx+1])
            string = string[idx+2:]

    return my_bounds

def find_the_bound(string):

    bound_list = []

    translation_table = str.maketrans('', '', "'")
    string = string.translate(translation_table)


    # gets string in expected form
    removeFront = True
    if string[0] != ' ':
        removeFront = False
    s_idx = 0
    while removeFront:
        if string[s_idx] == ' ':
            removeFront = False
        s_idx = s_idx + 1
    string = string[s_idx+1:-1]

    open_count = 0
    closed_count = 0
    my_string = ''
    skip = False
    for s in string:
        if s == '[':
            open_count = open_count + 1
        if s ==']':
            closed_count = closed_count + 1

        if open_count == 0 and closed_count == 0 and ( s == ' ' or s == ','):
           skip = False
        else:
            skip = True


        if s != '[' and s != ']' and skip:
            my_string = my_string + s

        if open_count == 1 and closed_count == 1:
            open_count = 0
            closed_count = 0
            my_string = my_string.split(',')
            bound_list.append(my_string)

            my_string = ''

    for i, bound in enumerate(bound_list):
        bound_list[i] = [float(bound[0]), float(bound[1])]

    return bound_list

def evaluate_bounds(bounds):
    # format = [[[],[]],[[],[],[]],[[],[],[],[],[],[],[]]]
    # list of list of lists
    if type(bounds) is not np.ndarray and type(bounds) is not list:
        if bounds[0] != '[' and bounds[-1] != ']':
            raise TypeError('String is not a list.')

        final_array = []
        if bounds != '[]':
            bounds = bounds[1:-1]
            bounds = find_each_bound(bounds)

            for bound in bounds:
                my_bound = find_the_bound(bound)
                final_array.append(my_bound)
    else:
        final_array = bounds

    return final_array

def evaluate_weights(weights):
    # format = [[],[],[]]
    weight_list = []
    final_array = []

    if type(weights) is not np.ndarray and type(weights) is not list:
        if weights != '[]':
            # removes uneccesary characters and spaces
            translation_table = str.maketrans('', '', "'")
            weights = weights.translate(translation_table)

            # gets string in expected form
            removeFront = True
            if weights[0] != ' ':
                removeFront = False
            s_idx = 0
            while removeFront:
                if weights[s_idx] == ' ':
                    removeFront = False
                s_idx = s_idx + 1
            weights = weights[s_idx + 1:-1]


            # finding the list of weights
            open_count = 0
            closed_count = 0
            my_string = ''
            skip = False
            for s in weights:


                if s == '[':
                    open_count = open_count + 1
                if s == ']':
                    closed_count = closed_count + 1

                if open_count == 0 and closed_count == 0 and (s == ' ' or s == ','):
                    skip = False
                else:
                    skip = True

                if s != '[' and s != ']' and skip:
                    my_string = my_string + s

                if open_count == 1 and closed_count == 1:
                    open_count = 0
                    closed_count = 0
                    my_string = my_string.split(',')
                    weight_list.append(my_string)

                    my_string = ''

        for weight in weight_list:
            temp_array = [float(w) for w in weight]

            final_array.append(temp_array)
    else:
        final_array = weights

    return final_array


def find_parameter_values(string):
    my_string = ''
    for s in string:
        if s != '[' and s != ']':
            my_string = my_string + s

    my_array = my_string.split(',')
    final_array = [float(my_array[0]), [float(my_array[1]), float(my_array[2])]]

    return final_array

def evaluate_parameters(parameters):
    # format = [[val, [lower,upper]],[val, [lower, upper]]]
    if type(parameters) is not np.ndarray and type(parameters) is not list:
        if parameters[0] != '[' and parameters[-1] != ']':
            raise TypeError('String is not a list.')

        final_array = []
        if parameters != '[]':
            parameters = parameters[1:-1]
            parameters = find_parameter_bound(parameters)

            for param in parameters:
                temp_param = find_parameter_values(param)
                final_array.append(temp_param)
    else:
        final_array = parameters

    return final_array

def saveSimulationHDF5(fname, sim_dict, version):
    """
    Purpose: Save the simulations
    :param fname: Filename to save the simulations to
    :param sim_dict: Simulation dictionary
    :return:
    """
    f = h5py.File(fname, 'a')  # create fname hdf5 file
    f.attrs['Version'] = version
    simulated = f['Simulated_data']

    simR = simulated['Reflectivity_Scan']
    simE = simulated['Energy_Scan']

    for name in list(sim_dict.keys()):

        if 'Angle' in list(sim_dict[name].keys()):
            dset = simE[name]
            m = np.shape(np.array(sim_dict[name]['Data']))

            dset.resize(m)
            dset[...] = np.array(sim_dict[name]['Data'])


            dset.attrs['DatasetNumber'] = sim_dict[name]['DatasetNumber']
            dset.attrs['DataPoints'] = sim_dict[name]['DataPoints']
            dset.attrs['Energy'] = sim_dict[name]['Energy']
            dset.attrs['Angle'] = sim_dict[name]['Angle']
            dset.attrs['Polarization'] = sim_dict[name]['Polarization']
            dset.attrs['Background Shift'] = sim_dict[name]['Background Shift']
            dset.attrs['Scaling Factor'] = sim_dict[name]['Scaling Factor']
        else:
            dset = simR[name]
            n = len(sim_dict[name]['Data'][0])
            m = np.shape(np.array(sim_dict[name]['Data']))


            dset.resize(m)
            dset[...] = np.array(sim_dict[name]['Data'])


            dset.attrs['DatasetNumber'] = sim_dict[name]['DatasetNumber']
            dset.attrs['DataPoints'] = sim_dict[name]['DataPoints']
            dset.attrs['Energy'] = sim_dict[name]['Energy']
            dset.attrs['Polarization'] = sim_dict[name]['Polarization']
            dset.attrs['Background Shift'] = sim_dict[name]['Background Shift']
            dset.attrs['Scaling Factor'] = sim_dict[name]['Scaling Factor']

def saveAsFileHDF5(fname, sample, data_dict, sim_dict, fit, optimization, version):
    """
    Purpose: Save workspace information with a specfied filename
    :param fname: filename
    :param sample: slab class
    :param data_dict: data dictionary
    :param sim_dict: simulation dictionary
    :param fit: Fitting parameters
    :param optimization: Optimization parameters
    :return:
    """

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # clears the data file if it already exists
    if os.path.exists(fname):
        f.clear()

    """
        Purpose: saves the GUI information to the current file
        :param fname: file path
        :param sample: slab class
        :param data_dict: data dictionary
        :param fit: fitting parameters
        :param optimization: optimization parameters
        :return:
        """

    WriteSampleHDF5(fname, sample, version)  # save the sample information

    f = h5py.File(fname, "a")

    simulated = f.create_group('Simulated_data')
    simR = simulated.create_group('Reflectivity_Scan')
    simE = simulated.create_group('Energy_Scan')

    experiment = f.create_group('Experimental_data')
    reflScan = experiment.create_group('Reflectivity_Scan')
    energyScan = experiment.create_group('Energy_Scan')

    # save the data from data dict
    for name in list(data_dict.keys()):
        if 'Angle' in data_dict[name].keys():
            dat = data_dict[name]['Data']
            dat = np.array(dat)
            m = np.shape(dat)
            dset = energyScan.create_dataset(name, m, data=dat, maxshape=(4,None), chunks=True)
            dset.attrs['DatasetNumber'] = data_dict[name]['DatasetNumber']
            dset.attrs['DataPoints'] = data_dict[name]['DataPoints']
            dset.attrs['Energy'] = data_dict[name]['Energy']
            dset.attrs['Angle'] = data_dict[name]['Angle']
            dset.attrs['Polarization'] = data_dict[name]['Polarization']
            dset.attrs['Background Shift'] = data_dict[name]['Background Shift']
            dset.attrs['Scaling Factor'] = data_dict[name]['Scaling Factor']

            dat1 = sim_dict[name]['Data']
            dset1 = simE.create_dataset(name, m, data=dat1,maxshape=(4,None), chunks=True)
            dset1.attrs['DatasetNumber'] = sim_dict[name]['DatasetNumber']
            dset1.attrs['DataPoints'] = sim_dict[name]['DataPoints']
            dset1.attrs['Energy'] = sim_dict[name]['Energy']
            dset1.attrs['Angle'] = sim_dict[name]['Angle']
            dset1.attrs['Polarization'] = sim_dict[name]['Polarization']
            dset1.attrs['Background Shift'] = sim_dict[name]['Background Shift']
            dset1.attrs['Scaling Factor'] = sim_dict[name]['Scaling Factor']
        else:
            dat = data_dict[name]['Data']
            dat = np.array(dat)
            m = np.shape(dat)

            dset = reflScan.create_dataset(name, m, data=dat, maxshape=(3,None), chunks=True)
            dset.attrs['DatasetNumber'] = data_dict[name]['DatasetNumber']
            dset.attrs['DataPoints'] = data_dict[name]['DataPoints']
            dset.attrs['Energy'] = data_dict[name]['Energy']
            dset.attrs['Polarization'] = data_dict[name]['Polarization']
            dset.attrs['Background Shift'] = data_dict[name]['Background Shift']
            dset.attrs['Scaling Factor'] = data_dict[name]['Scaling Factor']

            dat1 = sim_dict[name]['Data']
            dset1 = simR.create_dataset(name, m, data=dat1, maxshape=(3,None), chunks=True)
            dset1.attrs['DatasetNumber'] = sim_dict[name]['DatasetNumber']
            dset1.attrs['DataPoints'] = sim_dict[name]['DataPoints']
            dset1.attrs['Energy'] = sim_dict[name]['Energy']
            dset1.attrs['Polarization'] = sim_dict[name]['Polarization']
            dset1.attrs['Background Shift'] = sim_dict[name]['Background Shift']
            dset1.attrs['Scaling Factor'] = sim_dict[name]['Scaling Factor']


    sfBsFitParams = fit[0]
    sfBsVal = fit[1]
    sampleParams = fit[2]
    sampleVal = fit[3]
    scans = fit[4]
    bounds = fit[5]
    weights = fit[6]
    x = fit[7]
    fun = fit[8]


    opt = f.create_group("Optimization")

    diff_ev = opt.create_group("Differential Evolution")
    shgo = opt.create_group("Simplicial Homology")
    dual = opt.create_group("Dual Annealing")
    least = opt.create_group('Least Squares')
    """
    # How to include a new optimization algorithm? Let's use the direct algorithm as an example:
    direct = opt.create_group('Direct')  # create a group for the direct algorithm
    
    # Now we can include the attributes to the algorithm
    direct.attrs['eps'] = optimization['direct'][0]
    direct.attrs['maxFun'] = optimization['direct'][1]
    ...
    ...
    ...
    direct.attrs['last parameters'] = optimization['direct'][n]
    
    # Note that optimization is a dictionary that contains the algorithm name as the key and a list of the function 
    # variables as the values. This initialization of this dictionary can be found in GUI_GO in the initialization of 
    # the GlobalOptimizationWidget widget class. I will provide a little detail in GUI_GO of to include a new algorithm.
    """
    # load in optimization parameters for differential evolution
    diff_ev.attrs['strategy'] = optimization['differential evolution'][0]
    diff_ev.attrs['maxIter'] = optimization['differential evolution'][1]
    diff_ev.attrs['popsize'] = optimization['differential evolution'][2]
    diff_ev.attrs['tol'] = optimization['differential evolution'][3]
    diff_ev.attrs['atol'] = optimization['differential evolution'][4]
    diff_ev.attrs['min_mutation'] = optimization['differential evolution'][5]
    diff_ev.attrs['max_mutation'] = optimization['differential evolution'][6]
    diff_ev.attrs['recombination'] = optimization['differential evolution'][7]
    diff_ev.attrs['polish'] = str(optimization['differential evolution'][8])
    diff_ev.attrs['init'] = optimization['differential evolution'][9]
    diff_ev.attrs['updating'] = optimization['differential evolution'][10]

    # load in optimization parameters for simplicial homology
    if optimization['simplicial homology'][0] == None:
        shgo.attrs['n'] = str(optimization['simplicial homology'][0])
    else:
        shgo.attrs['n'] = optimization['simplicial homology'][0]

    shgo.attrs['iter'] = optimization['simplicial homology'][1]
    shgo.attrs['sampling'] = optimization['simplicial homology'][2]

    # load in optimization parameters for simplicial homology
    dual.attrs['maxiter'] = optimization['dual annealing'][0]
    dual.attrs['initial_temp'] = optimization['dual annealing'][1]
    dual.attrs['restart_temp'] = optimization['dual annealing'][2]
    dual.attrs['visit'] = optimization['dual annealing'][3]
    dual.attrs['accept'] = optimization['dual annealing'][4]
    dual.attrs['maxfun'] = optimization['dual annealing'][5]
    dual.attrs['local_search'] = str(optimization['dual annealing'][6])

    least.attrs['jac'] = optimization['least squares'][0]
    least.attrs['method'] = optimization['least squares'][1]
    least.attrs['ftol'] = optimization['least squares'][2]
    least.attrs['xtol'] = optimization['least squares'][3]
    least.attrs['gtol'] = optimization['least squares'][4]
    least.attrs['x_scale'] = optimization['least squares'][5]
    least.attrs['loss'] = optimization['least squares'][6]
    least.attrs['f_scale'] = optimization['least squares'][7]
    least.attrs['diff_step'] = optimization['least squares'][8]
    least.attrs['max_nfev'] = optimization['least squares'][9]

    # saving the fitting parameters
    fitting_parameters = f.create_group('Fitting Parameters')

    sample_fit = fitting_parameters.create_group('Sample Fit')
    scan_fit = fitting_parameters.create_group('Scan Fit')
    results = fitting_parameters.create_group('Results')

    sample_fit.attrs['sfbsFit'] = str(sfBsFitParams)
    sample_fit.attrs['sfbsVal'] = np.array(sfBsVal)
    sample_fit.attrs['Sample Fit'] = str(sampleParams)
    sample_fit.attrs['Sample Val'] = str(sampleVal)
    sample_fit.attrs['Selected Scans'] = str(scans)

    scan_fit.attrs['Selected Scans'] = str(scans)
    scan_fit.attrs['Bounds'] = str(bounds)
    scan_fit.attrs['Weights'] = str(weights)

    results.attrs['Value'] = str(x)
    results.attrs['Chi'] = fun

    experiment = f['Experimental_data']
    reflScan = experiment['Reflectivity_Scan']
    energyScan = experiment['Energy_Scan']

    # save the scaling factors and backgrounds shifts
    for idx, param in enumerate(sfBsFitParams):
        if param[0] == 'BACKGROUND SHIFT':  # background shift case
            if param[1] == 'ALL SCANS':  # all scans change
                for key in list(data_dict.keys()):
                    if 'Angle' in data_dict[param[1]].keys():
                        energyScan[key].attrs['Background Shift'] = float(sfBsVal[idx][0])
                    else:
                        reflScan[key].attrs['Background Shift'] = float(sfBsVal[idx][0])
            else:  # one scan changes
                if 'Angle' in data_dict[param[1]].keys():
                    energyScan[param[1]].attrs['Background Shift'] = float(sfBsVal[idx][0])
                else:
                    reflScan[param[1]].attrs['Background Shift'] = float(sfBsVal[idx][0])

        elif param[0] == 'SCALING FACTOR':  # scaling factor case
            if param[1] == 'ALL SCANS':  # all scans change
                for key in list(data_dict.keys()):
                    if 'Angle' in data_dict[param[1]].keys():
                        energyScan[key].attrs['Scaling Factor'] = float(sfBsVal[idx][0])
                    else:
                        reflScan[key].attrs['Scaling Factor'] = float(sfBsVal[idx][0])
            else:  # one scan changes
                if 'Angle' in data_dict[param[1]].keys():
                    energyScan[param[1]].attrs['Scaling Factor'] = float(sfBsVal[idx][0])
                else:
                    reflScan[param[1]].attrs['Scaling Factor'] = float(sfBsVal[idx][0])

    f.close()



def saveFileHDF5(fname, sample, data_dict, fit, optimization, version):
    """
    Purpose: saves the GUI information to the current file
    :param fname: file path
    :param sample: slab class
    :param data_dict: data dictionary
    :param fit: fitting parameters
    :param optimization: optimization parameters
    :return:
    """

    WriteSampleHDF5(fname, sample, version)  # save the sample information

    sfBsFitParams = fit[0]
    sfBsVal = fit[1]
    sampleParams = fit[2]
    sampleVal = fit[3]
    scans = fit[4]
    bounds = fit[5]
    weights = fit[6]
    x = fit[7]
    fun = fit[8]

    f = h5py.File(fname, "a")
    # initializes the optimization information
    opt = f["Optimization"]

    diff_ev = opt["Differential Evolution"]
    shgo = opt["Simplicial Homology"]
    dual = opt["Dual Annealing"]
    least = opt['Least Squares']

    # load in optimization parameters for differential evolution
    diff_ev.attrs['strategy'] = optimization['differential evolution'][0]
    diff_ev.attrs['maxIter'] = optimization['differential evolution'][1]
    diff_ev.attrs['popsize'] = optimization['differential evolution'][2]
    diff_ev.attrs['tol'] = optimization['differential evolution'][3]
    diff_ev.attrs['atol'] = optimization['differential evolution'][4]
    diff_ev.attrs['min_mutation'] = optimization['differential evolution'][5]
    diff_ev.attrs['max_mutation'] = optimization['differential evolution'][6]
    diff_ev.attrs['recombination'] = optimization['differential evolution'][7]
    diff_ev.attrs['polish'] = str(optimization['differential evolution'][8])
    diff_ev.attrs['init'] = optimization['differential evolution'][9]
    diff_ev.attrs['updating'] = optimization['differential evolution'][10]

    # load in optimization parameters for simplicial homology
    if optimization['simplicial homology'][0] == None:
        shgo.attrs['n'] = str(optimization['simplicial homology'][0])
    else:
        shgo.attrs['n'] = optimization['simplicial homology'][0]

    shgo.attrs['iter'] = optimization['simplicial homology'][1]
    shgo.attrs['sampling'] = optimization['simplicial homology'][2]

    # load in optimization parameters for simplicial homology
    dual.attrs['maxiter'] = optimization['dual annealing'][0]
    dual.attrs['initial_temp'] = optimization['dual annealing'][1]
    dual.attrs['restart_temp'] = optimization['dual annealing'][2]
    dual.attrs['visit'] = optimization['dual annealing'][3]
    dual.attrs['accept'] = optimization['dual annealing'][4]
    dual.attrs['maxfun'] = optimization['dual annealing'][5]
    dual.attrs['local_search'] = str(optimization['dual annealing'][6])

    least.attrs['jac'] = optimization['least squares'][0]
    least.attrs['method'] = optimization['least squares'][1]
    least.attrs['ftol'] = optimization['least squares'][2]
    least.attrs['xtol'] = optimization['least squares'][3]
    least.attrs['gtol'] = optimization['least squares'][4]
    least.attrs['x_scale'] = optimization['least squares'][5]
    least.attrs['loss'] = optimization['least squares'][6]
    least.attrs['f_scale'] = optimization['least squares'][7]
    least.attrs['diff_step'] = optimization['least squares'][8]
    least.attrs['max_nfev'] = optimization['least squares'][9]

    # saving the fitting parameters
    fitting_parameters = f['Fitting Parameters']

    sample_fit = fitting_parameters['Sample Fit']
    scan_fit = fitting_parameters['Scan Fit']
    results = fitting_parameters['Results']

    sample_fit.attrs['sfbsFit'] = str(sfBsFitParams)
    sample_fit.attrs['sfbsVal'] = str(sfBsVal)
    sample_fit.attrs['Sample Fit'] = str(sampleParams)
    sample_fit.attrs['Sample Val'] = str(sampleVal)
    sample_fit.attrs['Selected Scans'] = str(scans)

    scan_fit.attrs['Selected Scans'] = str(scans)
    scan_fit.attrs['Bounds'] = str(bounds)
    scan_fit.attrs['Weights'] = str(weights)

    results.attrs['Value'] = str(x)
    results.attrs['Chi'] = fun

    experiment = f['Experimental_data']
    reflScan = experiment['Reflectivity_Scan']
    energyScan = experiment['Energy_Scan']

    # save the scaling factors and backgrounds shifts
    for idx, param in enumerate(sfBsFitParams):
        if param[0] == 'BACKGROUND SHIFT':  # background shift case
            if param[1] == 'ALL SCANS':  # all scans change
                for key in list(data_dict.keys()):
                    if 'Angle' in data_dict[param[1]].keys():
                        energyScan[key].attrs['Background Shift'] = float(sfBsVal[idx][0])
                    else:
                        reflScan[key].attrs['Background Shift'] = float(sfBsVal[idx][0])
            else:  # one scan changes
                if 'Angle' in data_dict[param[1]].keys():
                    energyScan[param[1]].attrs['Background Shift'] = float(sfBsVal[idx][0])
                else:
                    reflScan[param[1]].attrs['Background Shift'] = float(sfBsVal[idx][0])

        elif param[0] == 'SCALING FACTOR':  # scaling factor case
            if param[1] == 'ALL SCANS':  # all scans change
                for key in list(data_dict.keys()):
                    if 'Angle' in data_dict[param[1]].keys():
                        energyScan[key].attrs['Scaling Factor'] = float(sfBsVal[idx][0])
                    else:
                        reflScan[key].attrs['Scaling Factor'] = float(sfBsVal[idx][0])
            else:  # one scan changes
                if 'Angle' in data_dict[param[1]].keys():
                    energyScan[param[1]].attrs['Scaling Factor'] = float(sfBsVal[idx][0])
                else:
                    reflScan[param[1]].attrs['Scaling Factor'] = float(sfBsVal[idx][0])

    f.close()

def newFileHDF5(fname, sample, version):
    """
        Purpose: Take in data and sample model and save it as an hdf5 file
        :param fname: File name with .hdf5 file extension
        :param AScans: Reflectivity scan data from ProcessRXR
        :param AInfo: Reflectivity scan information from ProcessRXR
        :param EScans: Energy scan data from Process RXR
        :param EInfo: Energy scan information from Process RXR
        :param sample: sample object
        :return:
        """

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # clears the data file if it already exists
    if os.path.exists(fname):
        f.clear()

    f.attrs['Version'] = version
    # creating group that will contain the sample information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)

    scattering_factor = sample.eShift
    mag_scattering_factor = sample.mag_eShift

    dsLayer = 0
    for my_layer in sample.structure:

        # Layer information
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            # Element information
            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # May be changed in the future as layermagnetized also contains this information
            # Original implemented to avoid problem of trying to load in the magnetic data that does not exist
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor



            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    grp1.attrs['FormFactors'] = str(scattering_factor)
    grp1.attrs['MagFormFactors'] = str(mag_scattering_factor)

    grp1.attrs['ffScale'] = str(sample.ff_scale)
    grp1.attrs['ffmScale'] = str(sample.ffm_scale)


    # initializes the experimental and simulation data (needs to be loaded later)
    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    grp4 = f.create_group("Optimization")

    diff_ev = grp4.create_group("Differential Evolution")
    shgo = grp4.create_group("Simplicial Homology")
    dual = grp4.create_group("Dual Annealing")
    least = grp4.create_group("Least Squares")

    # load in optimization parameters for differential evolution
    diff_ev.attrs['strategy'] = 'currenttobest1bin'
    diff_ev.attrs['maxIter'] = 150
    diff_ev.attrs['popsize'] = 15
    diff_ev.attrs['tol'] = 1e-6
    diff_ev.attrs['atol'] = 0
    diff_ev.attrs['min_mutation'] = 0.5
    diff_ev.attrs['max_mutation'] = 1
    diff_ev.attrs['recombination'] = 0.7
    diff_ev.attrs['polish'] = 'True'
    diff_ev.attrs['init'] = 'latinhypercube'
    diff_ev.attrs['updating'] = 'immediate'

    # load in optimization parameters for simplicial homology
    shgo.attrs['n'] = 'None'
    shgo.attrs['iter'] = 1
    shgo.attrs['sampling'] = 'simplicial'

    # load in optimization parameters for simplicial homology
    dual.attrs['maxiter'] = 150
    dual.attrs['initial_temp'] = 5230.0
    dual.attrs['restart_temp'] = 2e-5
    dual.attrs['visit'] = 2.62
    dual.attrs['accept'] = 5.0
    dual.attrs['maxfun'] = 10000000.0
    dual.attrs['local_search'] = 'False'

    least.attrs['jac'] = '2-point'
    least.attrs['method'] = 'trf'
    least.attrs['ftol'] = 1e-8
    least.attrs['xtol'] = 1e-8
    least.attrs['gtol'] = 1e-8
    least.attrs['x_scale'] = 1
    least.attrs['loss'] = 'linear'
    least.attrs['f_scale'] = 1
    least.attrs['diff_step'] = 'None'
    least.attrs['max_nfev'] = 'None'

    grp5 = f.create_group('Fitting Parameters')

    sample_fit = grp5.create_group('Sample Fit')
    scan_fit = grp5.create_group('Scan Fit')
    results = grp5.create_group('Results')

    sample_fit.attrs['sfbsFit'] = str([])
    sample_fit.attrs['sfbsVal'] = str([])
    sample_fit.attrs['Sample Fit'] = str([])
    sample_fit.attrs['Sample Val'] = str([])
    sample_fit.attrs['Selected Scans'] = str([])

    scan_fit.attrs['Selected Scans'] = str([])
    scan_fit.attrs['Bounds'] = str([])
    scan_fit.attrs['Weights'] = str([])

    results.attrs['Value'] = str([])
    results.attrs['Chi'] = 0


    f.close()

def createDataFileHDF5(fname, AScans, AInfo, EScans, EInfo):
    """
        Purpose: Take in data and sample model and save it as an hdf5 file
        :param fname: File name with .hdf5 file extension
        :param AScans: Reflectivity scan data from ProcessRXR
        :param AInfo: Reflectivity scan information from ProcessRXR
        :param EScans: Energy scan data from Process RXR
        :param EInfo: Energy scan information from Process RXR
        :param sample: sample object
        :return:
        """

    # Checking if fname already exists
    cwd = os.getcwd()
    path = cwd + '/' + fname
    if os.path.exists(path):
        raise OSError("HDF5 file already exists. To write new file remove the old file from the current working directory.")

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # Loading in the experimental data and simulated data
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp1 = f.create_group("Experimental_data")

    grpR = grp1.create_group("Reflectivity_Scan")
    grpE = grp1.create_group("Energy_Scan")


    # Loading reflectivity scans
    dsNum = 1
    for i in range(len(AScans)):
        qz = AScans[i][:, 2]  # momentum transfer
        R0 = AScans[i][:, 3]  # reflectivity
        energy = float(AInfo[i][3])  # energy of scan
        Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / np.pi  # initial angle
        dat = np.array([qz, Theta, R0])

        energy = float(AInfo[i][3])
        polarization = 'S'
        if (AInfo[i][1] == "S"):
            polarization = 'S'
        elif (AInfo[i][1] == "P"):
            polarization = 'P'
        elif (AInfo[i][1] == "L"):
            polarization = 'LC'
        elif (AInfo[i][1] == "R"):
            polarization = 'RC'
        datasetpoints = len(qz)

        name = str(AInfo[i][0]) + "_" + str(np.round(energy, 2)) + "_" + polarization

        dset = grpR.create_dataset(name, data=dat)

        dset.attrs['Energy'] = float(energy)

        # newly added
        dset.attrs['Background Shift'] = float(0)

        # newly added
        dset.attrs['Scaling Factor'] = float(1)

        dset.attrs['Polarization'] = str(polarization)

        dset.attrs['DataPoints'] = int(datasetpoints)

        dset.attrs['DatasetNumber'] = int(dsNum)

        dsNum = dsNum + 1
        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):

                qz = AScans[i][:, 2]
                A = (AScans[i - 1][:, 3] - AScans[i][:, 3]) / (AScans[i - 1][:, 3]) + AScans[i][:, 3]
                energy = float(AInfo[i][3])
                Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / np.pi  # initial angle
                dat = np.array([qz, Theta, A])

                datasetpoints = len(qz)

                if (AInfo[i - 1][1] == "S" or AInfo[i - 1][1] == "P"):
                    polarization = "AL"
                elif (AInfo[i - 1][1] == "L" or AInfo[i - 1][1] == "R"):
                    polarization = "AC"
                name = str(AInfo[i - 1][0]) + "-" + str(AInfo[i][0]) + "_" + str(
                    np.round(energy, 2)) + "_" + polarization + "_Asymm"
                dset = grpR.create_dataset(name, data=dat)

                dset.attrs['Energy'] = float(energy)

                # newly added
                dset.attrs['Background Shift'] = float(0)

                # newly added
                dset.attrs['Scaling Factor'] = float(1)

                dset.attrs['Polarization'] = polarization

                dset.attrs['DataPoints'] = int(datasetpoints)

                dset.attrs['DatasetNumber'] = int(dsNum)
                dsNum = dsNum + 1

    for i in range(len(EScans)):

        qz = EScans[i][:, 2]
        R0 = EScans[i][:, 3]
        E = EScans[i][:, 0]

        Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
        datasetpoints = len(qz)
        dat = np.array([qz, Theta, R0, E])

        energy = float(EInfo[i][3])
        polarization = 'S'
        if (EInfo[i][1] == "S"):
            polarization = 'S'
        elif (EInfo[i][1] == "P"):
            polarization = 'P'
        elif (EInfo[i][1] == "L"):
            polarization = 'LC'
        elif (EInfo[i][1] == "R"):
            polarization = 'RC'
        angle = float(EInfo[i][4])

        name = str(EInfo[i][0]) + "_E" + str(np.round(energy, 2)) + "_Th" + str(np.round(angle, 2)) + "_" + polarization

        dset = grpE.create_dataset(name, data=dat)

        dset.attrs['Energy'] = float(energy)

        # newly added
        dset.attrs['Background Shift'] = float(0)

        # newly added
        dset.attrs['Scaling Factor'] = float(1)

        dset.attrs['Polarization'] = polarization

        dset.attrs['DataPoints'] = int(datasetpoints)

        dset.attrs['Angle'] = float(angle)

        dset.attrs['DatasetNumber'] = int(dsNum)
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
                qz = EScans[i][:, 2]
                E = EScans[i][:, 0]
                A = (EScans[i - 1][:, 3] - EScans[i][:, 3]) / (EScans[i - 1][:, 3] + EScans[i][:, 3])
                Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
                dat = np.array([qz, Theta, A, E])
                energy = float(EInfo[i][3])

                if (EInfo[i - 1][1] == "S" or EInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                elif (EInfo[i - 1][1] == "L" or EInfo[i - 1][1] == "R"):
                    polarization = 'AC'

                angle = float(EInfo[i][4])
                name = str(EInfo[i - 1][0]) + "-" + str(EInfo[i][0]) + "_E" + str(np.round(energy, 2)) + "_Th" + str(
                    np.round(angle, 2)) + "_" + polarization + "_Asymm"

                dset = grpE.create_dataset(name, data=dat)

                dset.attrs['Energy'] = float(energy)

                # newly added
                dset.attrs['Background Shift'] = float(0)

                # newly added
                dset.attrs['Scaling Factor'] = float(1)

                dset.attrs['Polarization'] = polarization

                dset.attrs['Angle'] = float(angle)

                dset.attrs['DataPoints'] = int(datasetpoints)

                dset.attrs['DatasetNumber'] = int(dsNum)

                dsNum = dsNum + 1

    f.close()

def WriteDataHDF5(fname, AScans,AInfo, EScans, EInfo, sample):
    """
    Purpose: Take in data and sample model and save it as an hdf5 file
    :param fname: File name with .hdf5 file extension
    :param AScans: Reflectivity scan data from ProcessRXR
    :param AInfo: Reflectivity scan information from ProcessRXR
    :param EScans: Energy scan data from Process RXR
    :param EInfo: Energy scan information from Process RXR
    :param sample: sample object
    :return:
    """

    # Checking if fname already exists
    cwd = os.getcwd()
    path = cwd + '/' + fname
    if os.path.exists(path):
        raise OSError("HDF5 file already exists. To write new file remove the old file from the current working directory.")

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # creating group that will contain the sample information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)


    dsLayer = 0
    for my_layer in sample.structure:

        # Layer information
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            # Element information
            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # May be changed in the future as layermagnetized also contains this information
            # Original implemented to avoid problem of trying to load in the magnetic data that does not exist
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor



            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    grp1.attrs['FormFactors'] = str(sample.eShift)
    grp1.attrs['MagFormFactors'] = str(sample.mag_eShift)

    grp1.attrs['ffScale'] = str(sample.ff_scale)
    grp1.attrs['ffmScale'] = str(sample.ffm_scale)

    # Loading in the experimental data and simulated data
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    #Loading reflectivity scans
    dsNum = 1
    for i in range(len(AScans)):
        qz = AScans[i][:,2]  # momentum transfer
        R0 = AScans[i][:,3]  # reflectivity
        energy = float(AInfo[i][3])  # energy of scan
        Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / np.pi  # initial angle
        dat = np.array([qz, Theta, R0])

        energy = float(AInfo[i][3])
        polarization = 'S'
        if (AInfo[i][1] == "S"):
            polarization = 'S'
        elif (AInfo[i][1] == "P"):
            polarization = 'P'
        elif (AInfo[i][1] == "L"):
            polarization = 'LC'
        elif (AInfo[i][1] == "R"):
            polarization = 'RC'
        datasetpoints = len(qz)

        name = str(AInfo[i][0]) + "_" + str(np.round(energy,2)) + "_" + polarization
        qz, R = sample.reflectivity(energy, qz)
        sim = np.array([qz, Theta, R[polarization]])
        dset = grpR.create_dataset(name, data=dat)
        dset1 = subR.create_dataset(name, data=sim)

        dset.attrs['Energy'] = float(energy)
        dset1.attrs['Energy'] = float(energy)

        # newly added
        dset.attrs['Background Shift'] = float(0)
        dset1.attrs['Background Shift'] = float(0)

        # newly added
        dset.attrs['Scaling Factor'] = float(1)
        dset1.attrs['Scaling Factor'] = float(1)

        dset.attrs['Polarization'] = str(polarization)
        dset1.attrs['Polarization'] = str(polarization)

        dset.attrs['DataPoints'] = int(datasetpoints)
        dset1.attrs['DataPoints'] = int(datasetpoints)

        dset.attrs['DatasetNumber'] = int(dsNum)
        dset1.attrs['DatasetNumber'] = int(dsNum)

        dsNum = dsNum + 1
        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):

                qz = AScans[i][:,2]
                A = (AScans[i - 1][:,3] - AScans[i][:,3]) / (AScans[i - 1][:,3]) + AScans[i][:,3]
                energy = float(AInfo[i][3])
                Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / np.pi  # initial angle
                dat = np.array([qz, Theta, A])

                datasetpoints = len(qz)

                if (AInfo[i - 1][1] == "S" or AInfo[i - 1][1] == "P"):
                    polarization = "AL"
                elif (AInfo[i - 1][1] == "L" or AInfo[i - 1][1] == "R"):
                    polarization = "AC"
                name = str(AInfo[i-1][0])+ "-" + str(AInfo[i][0]) + "_" + str(np.round(energy,2)) +"_"+ polarization + "_Asymm"
                dset = grpR.create_dataset(name, data=dat)
                qz, R = sample.reflectivity(energy, qz)
                sim = np.array([qz,Theta,R[polarization]])
                dset1 = subR.create_dataset(name,data = sim)

                dset.attrs['Energy'] = float(energy)
                dset1.attrs['Energy'] = float(energy)

                # newly added
                dset.attrs['Background Shift'] = float(0)
                dset1.attrs['Background Shift'] = float(0)

                # newly added
                dset.attrs['Scaling Factor'] = float(1)
                dset1.attrs['Scaling Factor'] = float(1)

                dset.attrs['Polarization'] = polarization
                dset1.attrs['Polarization'] = polarization

                dset.attrs['DataPoints'] = int(datasetpoints)
                dset1.attrs['DataPoints'] = int(datasetpoints)

                dset.attrs['DatasetNumber'] = int(dsNum)
                dset1.attrs['DatasetNumber'] = int(dsNum)
                dsNum = dsNum + 1


    for i in range(len(EScans)):

        qz = EScans[i][:,2]
        R0 = EScans[i][:,3]
        E = EScans[i][:,0]

        Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
        datasetpoints = len(qz)
        dat = np.array([qz,Theta,R0,E])

        energy = float(EInfo[i][3])
        polarization = 'S'
        if (EInfo[i][1] == "S"):
            polarization = 'S'
        elif (EInfo[i][1] == "P"):
            polarization = 'P'
        elif (EInfo[i][1] == "L"):
            polarization = 'LC'
        elif (EInfo[i][1] == "R"):
            polarization = 'RC'
        angle = float(EInfo[i][4])

        name = str(EInfo[i][0]) + "_E" + str(np.round(energy,2)) + "_Th" + str(np.round(angle,2)) + "_" + polarization
        E, R = sample.energy_scan(angle,E)
        sim = np.array([qz, Theta, R[polarization], E])

        dset = grpE.create_dataset(name, data=dat)
        dset1 = subE.create_dataset(name, data=sim)

        dset.attrs['Energy'] = float(energy)
        dset1.attrs['Energy'] = float(energy)

        # newly added
        dset.attrs['Background Shift'] = float(0)
        dset1.attrs['Background Shift'] = float(0)

        # newly added
        dset.attrs['Scaling Factor'] = float(1)
        dset1.attrs['Scaling Factor'] = float(1)

        dset.attrs['Polarization'] = polarization
        dset1.attrs['Polarization'] = polarization

        dset.attrs['DataPoints'] = int(datasetpoints)
        dset1.attrs['DataPoints'] = int(datasetpoints)

        dset.attrs['Angle'] = float(angle)
        dset1.attrs['Angle'] = float(angle)

        dset.attrs['DatasetNumber'] = int(dsNum)
        dset1.attrs['DatasetNumber'] = int(dsNum)
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
                qz = EScans[i][:,2]
                E = EScans[i][:,0]
                A = (EScans[i-1][:,3]-EScans[i][:,3])/(EScans[i-1][:,3]+EScans[i][:,3])
                Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
                dat = np.array([qz, Theta, A,E])
                energy = float(EInfo[i][3])

                if (EInfo[i - 1][1] == "S" or EInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                elif (EInfo[i - 1][1] == "L" or EInfo[i - 1][1] == "R"):
                    polarization = 'AC'

                angle = float(EInfo[i][4])
                name = str(EInfo[i-1][0])+"-"+ str(EInfo[i][0]) + "_E" + str(np.round(energy,2)) + "_Th" + str(np.round(angle, 2)) + "_" + polarization + "_Asymm"
                E, R = sample.energy_scan(angle, E)
                sim = np.array([qz, Theta, R[polarization], E])

                dset = grpE.create_dataset(name, data=dat)
                dset1 = subE.create_dataset(name, data=sim)

                dset.attrs['Energy'] = float(energy)
                dset1.attrs['Energy'] = float(energy)

                # newly added
                dset.attrs['Background Shift'] = float(0)
                dset1.attrs['Background Shift'] = float(0)

                # newly added
                dset.attrs['Scaling Factor'] = float(1)
                dset1.attrs['Scaling Factor'] = float(1)

                dset.attrs['Polarization'] = polarization
                dset1.attrs['Polarization'] = polarization

                dset.attrs['Angle'] = float(angle)
                dset1.attrs['Angle'] = float(angle)

                dset.attrs['DataPoints'] = int(datasetpoints)
                dset1.attrs['DataPoints'] = int(datasetpoints)

                dset.attrs['DatasetNumber'] = int(dsNum)
                dset1.attrs['DatasetNumber'] = int(dsNum)

                dsNum = dsNum + 1

    grp4 = f.create_group("Optimization")

    diff_ev = grp4.create_group("Differential Evolution")
    shgo = grp4.create_group("Simplicial Homology")
    dual = grp4.create_group("Dual Annealing")
    least = grp4.create_group("Least Squares")

    # load in optimization parameters for differential evolution
    diff_ev.attrs['strategy'] = 'currenttobest1bin'
    diff_ev.attrs['maxIter'] = 150
    diff_ev.attrs['popsize'] = 15
    diff_ev.attrs['tol'] = 1e-6
    diff_ev.attrs['atol'] = 0
    diff_ev.attrs['min_mutation'] = 0.5
    diff_ev.attrs['max_mutation'] = 1
    diff_ev.attrs['recombination'] = 0.7
    diff_ev.attrs['polish'] = 'True'
    diff_ev.attrs['init'] = 'latinhypercube'
    diff_ev.attrs['updating'] = 'immediate'

    # load in optimization parameters for simplicial homology
    shgo.attrs['n'] = 'None'
    shgo.attrs['iter'] = 1
    shgo.attrs['sampling'] = 'simplicial'

    # load in optimization parameters for simplicial homology
    dual.attrs['maxiter'] = 150
    dual.attrs['initial_temp'] = 5230.0
    dual.attrs['restart_temp'] = 2e-5
    dual.attrs['visit'] = 2.62
    dual.attrs['accept'] = 5.0
    dual.attrs['maxfun'] = 10000000.0
    dual.attrs['local_search'] = 'False'

    least.attrs['jac'] = '2-point'
    least.attrs['method'] = 'trf'
    least.attrs['ftol'] = 1e-8
    least.attrs['xtol'] = 1e-8
    least.attrs['gtol'] = 1e-8
    least.attrs['x_scale'] = 1
    least.attrs['loss'] = 'linear'
    least.attrs['f_scale'] = 1
    least.attrs['diff_step'] = 'None'
    least.attrs['max_nfev'] = 'None'

    grp5 = f.create_group('Fitting Parameters')

    sample_fit = grp5.create_group('Sample Fit')
    scan_fit = grp5.create_group('Scan Fit')
    results = grp5.create_group('Results')

    sample_fit.attrs['sfbsFit'] = str([])
    sample_fit.attrs['sfbsVal'] = str([])
    sample_fit.attrs['Sample Fit'] = str([])
    sample_fit.attrs['Sample Val'] = str([])
    sample_fit.attrs['Selected Scans'] = str([])

    scan_fit.attrs['Selected Scans'] = str([])
    scan_fit.attrs['Bounds'] = str([])
    scan_fit.attrs['Weights'] = str([])

    results.attrs['Value'] = str([])
    results.attrs['Chi'] = 0

    f.close()

def ReadFitHDF5(fname):
    """
        Purpose: Read in the algorithm parameters
        :param fname: File name
        :return: the algorithms parameters
        """

    f = h5py.File(fname, 'r')

    fitting_parameters = f['Fitting Parameters']

    sample_fit = fitting_parameters['Sample Fit']
    scan_fit = fitting_parameters['Scan Fit']
    results = fitting_parameters['Results']

    sfbsFit = ast.literal_eval(sample_fit.attrs['sfbsFit'])
    sfbsVal = evaluate_parameters(sample_fit.attrs['sfbsVal'])
    sampleFit = ast.literal_eval(sample_fit.attrs['Sample Fit'])
    sampleVal = evaluate_parameters(sample_fit.attrs['Sample Val'])
    selected_scans = ast.literal_eval(sample_fit.attrs['Selected Scans'])
    bounds = evaluate_bounds(scan_fit.attrs['Bounds'])
    weights = evaluate_weights(scan_fit.attrs['Weights'])

    x = evaluate_list(results.attrs['Value'])

    chi = results.attrs['Chi']

    f.close()
    return [sfbsFit, sfbsVal, sampleFit, sampleVal, selected_scans, bounds, weights, x, chi]


def ReadAlgorithmHDF5(fname):
    """
    Purpose: Read in the algorithm parameters
    :param fname: File name
    :return: the algorithms parameters
    """

    f = h5py.File(fname, 'r')

    parameters = dict()
    optimization = f['Optimization']

    diff_ev = optimization["Differential Evolution"]
    shgo = optimization["Simplicial Homology"]
    dual = optimization["Dual Annealing"]
    least = optimization['Least Squares']

    # load in optimization parameters for differential evolution
    eStrategy = diff_ev.attrs['strategy']
    eMaxIter = diff_ev.attrs['maxIter']
    ePopsize = diff_ev.attrs['popsize']
    eTol = diff_ev.attrs['tol']
    eAtol = diff_ev.attrs['atol']
    eMinMutation = diff_ev.attrs['min_mutation']
    eMaxMutation = diff_ev.attrs['max_mutation']
    eRecombination = diff_ev.attrs['recombination']
    ePolish = diff_ev.attrs['polish']
    eInit = diff_ev.attrs['init']
    eUpdating = diff_ev.attrs['updating']

    if ePolish == 'True':
        ePolish = True
    elif ePolish == 'False':
        ePolish = False

    parameters['differential evolution'] = [eStrategy, eMaxIter, ePopsize, eTol, eAtol, eMinMutation, eMaxMutation,
                                            eRecombination, ePolish, eInit, eUpdating]

    # load in optimization parameters for simplicial homology
    sN = shgo.attrs['n']
    sIter = shgo.attrs['iter']
    sSampling = shgo.attrs['sampling']

    if sN == 'None':
        sN = None

    parameters['simplicial homology'] = [sN, sIter, sSampling]

    # load in optimization parameters for simplicial homology
    dMaxiter = dual.attrs['maxiter']
    dInitTemp = dual.attrs['initial_temp']
    dRestTemp = dual.attrs['restart_temp']
    dVisit = dual.attrs['visit']
    dAccept = dual.attrs['accept']
    dMaxfun = dual.attrs['maxfun']
    dLs = dual.attrs['local_search']

    if dLs == 'True':
        dLs = True
    elif dLs == 'False':
        dLs = False

    parameters['dual annealing'] = [dMaxiter, dInitTemp, dRestTemp, dVisit, dAccept, dMaxfun, dLs]

    jac = least.attrs['jac']
    method = least.attrs['method']
    ftol = least.attrs['ftol']
    xtol = least.attrs['xtol']
    gtol = least.attrs['gtol']
    x_scale = least.attrs['x_scale']
    loss = least.attrs['loss']
    f_scale = least.attrs['f_scale']
    diff_step = least.attrs['diff_step']
    max_nfev = least.attrs['max_nfev']

    parameters['least squares'] = [jac, method, ftol, xtol, gtol, x_scale, loss, f_scale,
                                   diff_step, max_nfev]

    f.close()
    return parameters

def ReadSampleHDF5(fname):

    """
    Purpose: Read the sample info from hdf5 file and recreate the sample object
    :param fname: File name
    :return: sample - the recreated sample object
    """
    f = h5py.File(fname, 'r')

    S = f['Sample']

    # Retieves the general information of the sample
    m = int(S.attrs['NumberLayers'])
    sample = slab(m)

    if S.attrs['PolyElements'] == '{}':
        sample.poly_elements = dict()
    else:
        sample.poly_elements = ast.literal_eval(S.attrs['PolyElements'])

    if S.attrs['MagElements'] == '{}':
        sample.mag_elements = dict()
    else:
        sample.mag_elements = ast.literal_eval(S.attrs['MagElements'])

    sample.mag_elements = dict()

    sample.layer_magnetized = S.attrs['LayerMagnetized']
    sample.eShift = ast.literal_eval(S.attrs['FormFactors'])
    sample.mag_eShift = ast.literal_eval(S.attrs['MagFormFactors'])

    sample.ff_scale = ast.literal_eval(S.attrs['ffScale'])
    sample.ffm_scale = ast.literal_eval(S.attrs['ffmScale'])
    sample.find_sf = ast.literal_eval(S.attrs['findFF'])


    # Retrieves the general layer information
    for lay_key in S.keys():
        layer = S[lay_key]
        formula = layer.attrs['Formula']
        lay_num = int(layer.attrs['LayerNumber'])
        sample.addlayer(lay_num, formula,20, density=1)  #pre-initialize parameters to random numbers for each layer

        # retrieves the information for each element
        for ele_key in layer.keys():
            element = layer[ele_key]
            sample.structure[lay_num][ele_key].molar_mass = element.attrs['MolarMass']
            sample.structure[lay_num][ele_key].density = element.attrs['Density']
            sample.structure[lay_num][ele_key].thickness = element.attrs['Thickness']
            sample.structure[lay_num][ele_key].roughness = element.attrs['Roughness']
            sample.structure[lay_num][ele_key].linked_roughness = element.attrs['LinkedRoughness']
            sample.structure[lay_num][ele_key].poly_ratio = element.attrs['PolyRatio']
            sample.structure[lay_num][ele_key].polymorph = element.attrs['Polymorph']
            sample.structure[lay_num][ele_key].gamma = element.attrs['Gamma']
            sample.structure[lay_num][ele_key].phi = element.attrs['Phi']
            if element.attrs['Magnetic']:

                if ele_key not in list(sample.mag_elements.keys()):
                    if ele_key not in list(sample.poly_elements):
                        sample.mag_elements[ele_key] = ele_key
                    else:
                        sample.mag_elements[ele_key] = element.attrs['Polymorph']
                sample.structure[lay_num][ele_key].mag_density = element.attrs['MagDensity']
                sample.structure[lay_num][ele_key].mag_scattering_factor = element.attrs['MagScatteringFactor']


            if element.attrs['ScatteringFactor'][0] == '[':
                sample.structure[lay_num][ele_key].scattering_factor = ast.literal_eval(element.attrs['ScatteringFactor'])
            else:
                sample.structure[lay_num][ele_key].scattering_factor = element.attrs['ScatteringFactor']
            sample.structure[lay_num][ele_key].position = element.attrs['Position']


    return sample

def ReadDataHDF5(fname):
    """
    Purpose: Reads in the experimental and simulated data from hdf5 file and then plots spectrum chosen by user
    :param fname: File name
    :return:
    """

    f = h5py.File(fname, 'r')
    experiment = f['Experimental_data']
    simulation = f['Simulated_data']

    RS = experiment['Reflectivity_Scan']
    Rsim = simulation['Reflectivity_Scan']

    ES = experiment['Energy_Scan']
    Esim = simulation['Energy_Scan']

    # Collects data information to print to terminal
    data = list()
    data_dict = dict()
    sim_dict = dict()

    for Rkey in RS.keys():
        mydata = RS[Rkey]
        data_dict[Rkey] = RS[Rkey]
        sim_dict[Rkey] = Rsim[Rkey]
        Rdat = [int(mydata.attrs['DatasetNumber']),'Reflectivity', Rkey]
        data.append(Rdat)

    for Ekey in ES.keys():
        mydata = ES[Ekey]
        data_dict[Ekey] = ES[Ekey]
        sim_dict[Ekey] = Esim[Ekey]
        Edat = [int(mydata.attrs['DatasetNumber']),'Energy', Ekey]
        data.append(Edat)

    # Sorts data in appropriate order
    data = np.array(data)

    sort_idx = np.argsort(data[:,0].astype(int))
    data = data[sort_idx]



    data_dict = hdf5ToDict(data_dict)
    sim_dict = hdf5ToDict(sim_dict)
    f.close()

    return data, data_dict, sim_dict



def LoadDataHDF5(fname):
    """
    Purpose: Reads in the experimental and simulated data from hdf5 file and then plots spectrum chosen by user
    :param fname: File name
    :return:
    """

    f = h5py.File(fname, 'r')
    experiment = f['Experimental_data']
    simulation = f['Simulated_data']

    RS = experiment['Reflectivity_Scan']
    ES = experiment['Energy_Scan']

    RSim = simulation['Reflectivity_Scan']
    ESim = simulation['Energy_Scan']
    # Collects data information to print to terminal
    data = list()
    data_dict = dict()
    sim_dict = dict()

    for Rkey in RS.keys():
        mydata = RS[Rkey]
        data_dict[Rkey] = RS[Rkey]
        sim_dict[Rkey] = RSim[Rkey]
        Rdat = [int(mydata.attrs['DatasetNumber']),'Reflectivity', Rkey]
        data.append(Rdat)

    for Ekey in ES.keys():
        mydata = ES[Ekey]
        data_dict[Ekey] = ES[Ekey]
        sim_dict[Ekey] = ESim[Ekey]
        Edat = [int(mydata.attrs['DatasetNumber']),'Energy', Ekey]
        data.append(Edat)

    # Sorts data in appropriate order
    data = np.array(data)
    sort_idx = np.argsort(data[:,0].astype(int))
    data = data[sort_idx]



    data_dict = hdf5ToDict(data_dict)
    sim_dict = hdf5ToDict(sim_dict)
    f.close()

    return data, data_dict, sim_dict

def WriteSampleHDF5(fname, sample, version):
    """
    Purpose: Write a new sample to the hdf5 file fname
    :param fname: File name
    :param sample: new sample information as sample type
    :return:
    """

    # deletes the previous sample information
    with h5py.File(fname, "a") as f:
        if 'Sample' in f:
            del f['Sample']
    f.close()

    f = h5py.File(fname, "a")
    f.attrs['Version'] = version
    grp1 = f.create_group('Sample')

    for key in list(sample.poly_elements):
        sample.poly_elements[key] = list(sample.poly_elements)
    for key in list(sample.mag_elements):
        sample.mag_elements[key] = list(sample.mag_elements)

    m = len(sample.structure)
    grp1.attrs['findFF'] = str(sample.find_sf)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)

    scattering_factor = sample.eShift
    mag_scattering_factor = sample.mag_eShift

    grp1.attrs['FormFactors'] = str(scattering_factor)
    grp1.attrs['MagFormFactors'] = str(mag_scattering_factor)

    grp1.attrs['ffScale'] = str(sample.ff_scale)
    grp1.attrs['ffmScale'] = str(sample.ffm_scale)

    # Sets the information for each layer
    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)


        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            # Remove the numbers in the element definition
            new_ele = ''
            for c in ele:
                if not(c.isdigit()):
                    new_ele = new_ele + c



            if stoich == 1:
                formula = formula + new_ele
            else:
                formula = formula + new_ele + str(stoich)



        layer.attrs['Formula'] = formula

        # Sets the information for each element
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True

                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = list(my_layer[ele].scattering_factor)

            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    f.close()


def WriteHDF5Simulation(sample):
    """
        Purpose: Write a new sample to the hdf5 file fname
        :param fname: File name
        :param sample: new sample information as sample type
        :return:
        """

    # deletes the previous sample information
    with h5py.File(fname, "a") as f:
        if 'Sample' in f:
            del f['Sample']
    f.close()

    f = h5py.File(fname, "a")
    grp1 = f.create_group('Sample')

    # Sets the general sample information
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)

    scattering_factor = sample.eShift
    mag_scattering_factor = sample.mag_eShift

    grp1.attrs['FormFactors'] = str(scattering_factor)
    grp1.attrs['MagFormFactors'] = str(mag_scattering_factor)

    grp1.attrs['ffScale'] = str(sample.ff_scale)
    grp1.attrs['ffmScale'] = str(sample.ffm_scale)

    # Sets the information for each layer
    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)
        ### Need to change back to previous version

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        # Sets the information for each element
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f["Simulated_data"]

    grpR = grp2["Reflectivity_Scan"]
    grpE = grp2["Energy_Scan"]

    # Recalculate the simulated reflectivity scan data
    for key in list(grpR.keys()):
        dset = grpR[key]  # retrieves the old data set
        Rdata = list(dset)  # gets the dataset information
        qz = np.array(Rdata[0])  # retrieves the momentum transfer
        E = float(dset.attrs['Energy'])  # retrieves the energy
        polarization = dset.attrs['Polarization']  # retrieves the polarization

        qz, R = sample.reflectivity(E, qz)  # computes reflecticity
        Rdata[2] = R[polarization]  # retrieves appropriate reflectivity for the correct polarization
        dset[...] = Rdata  # overwrites previous dataset with new data

    # Recalculates the simulated energy scan data
    for key in list(grpE.keys()):
        dset = grpE[key]  # retrieves the old data set
        Edata = list(dset)  # gets the dataset information
        E = Edata[3]  # retrieves the numpy energy array
        theta = dset.attrs['Angle']  # retrieves the angle
        polarization = dset.attrs['Polarization']  # retrieves the polarization
        E, R = sample.energy_scan(theta, E)  # recomputes the energy scan
        Edata[2] = R[polarization]  # retrieves appropriate reflectivity with the correct polarization
        dset[...] = Edata  # overwrites previous data with new data

    f.close()

def updateHDF5Data(fname, scalingFactor, backgroundShift):

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    f = h5py.File(fname, 'r')
    experiment = f['Experimental_data']

    RS = experiment['Reflectivity_Scan']

    ES = experiment['Energy_Scan']


    # Recalculate the simulated reflectivity scan data
    idx = 0
    for key in list(RS.keys()):
        RS[key].attrs['ScalingFactor'] = scalingFactor[idx]
        RS[key].attrs['BackgroundShift'] = backgroundShift[idx]
        idx = idx + 1

    # Recalculates the simulated energy scan data
    for key in list(ES.keys()):
        ES[key].attrs['ScalingFactor'] = scalingFactor[idx]
        ES[key].attrs['BackgroundShift'] = backgroundShift[idx]
        idx = idx + 1

    f.close()


def WriteSampleASCII(file,sample):
    """
    Purpose: Write the sample as an ASCII file
    :param file: File name
    :param sample: sample information as a slab object
    :return:
    """

    file.write("# Structure \n")  # header defining that the sample information is starting
    n = len(sample.structure)

    # General information for the sample model
    file.write("numberlayers = %s \n" % str(n))
    file.write("polyelements = %s \n" % str(sample.poly_elements))
    file.write("magelements = %s \n" % str(sample.mag_elements))
    file.write("layermagnetized = %s \n" % sample.layer_magnetized)
    file.write("scalingfactor = %s \n" % str(sample.scaling_factor))
    file.write("backgroundshift = %s \n\n" % str(sample.background_shift))

    # writing the layer and element information
    num_lay = 0
    for layer in sample.structure:

        # General layer information
        file.write("layer = %s \n" % str(num_lay))

        # Reconstructing the chemical formula
        formula = ''
        for ele in list(layer.keys()):
            stoich = layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        file.write("formula = %s \n\n" % formula)

        # writing the element information
        for ele in layer.keys():
            file.write("element = %s \n" % ele)
            file.write("molarmass = %f \n" % layer[ele].molar_mass)
            file.write("density = %f \n" % layer[ele].density)
            file.write("thickness = %f \n" % layer[ele].thickness)
            file.write("roughness = %f \n" % layer[ele].roughness)
            file.write("linkedroughness = %f \n" % layer[ele].linked_roughness)
            file.write("scatteringfactor = %s \n" % layer[ele].scattering_factor)
            file.write("polymorph = %s \n" % layer[ele].polymorph)

            poly_ratio = layer[ele].poly_ratio
            if type(poly_ratio) != int:
                poly_ratio = [str(poly) for poly in poly_ratio]
            file.write("polyratio = %s \n" % poly_ratio)


            file.write("gamma = %f \n" % layer[ele].gamma)
            file.write("phi = %f \n" % layer[ele].phi)

            mag_density = layer[ele].mag_density
            mag_density = [str(x) for x in mag_density]
            file.write("magdensity = %s \n" % mag_density)
            sfm = layer[ele].mag_scattering_factor
            file.write("magscatteringfactor = %s \n" % sfm)
            file.write("position = %s \n" % layer[ele].position)
            file.write("\n")

        num_lay = num_lay + 1


def WriteExperimentalDataASCII(file, AScans,AInfo,EScans,EInfo):
    """
    Purpose: Write experimental data to .all file for program use
    :param fname: Name of data file
    :param AScans: Reflection scan experimental data
    :param AInfo: Info about reflection data scans
    :param EScans: Energy scan experimental data
    :param EInfo: Energy scan info
    :param header: Additional header
    :return:
    """


    file.write("# Experimental_Data \n")  # header that states writing experimental data

    # Writing the reflectivity scan data
    dsNum = 1
    for i in range(len(AScans)):

        # Writing the general information of each scan
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % AInfo[i][3])
        if (AInfo[i][1] == "S"):
            file.write("polarization = S \n")
        elif (AInfo[i][1] == "P"):
            file.write("polarization = P \n")
        elif (AInfo[i][1] == "L"):
            file.write("polarization = LC \n")
        elif (AInfo[i][1] == "R"):
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

        # writing the data of each scan
        for j in range(len(AScans[i][:, 0])):
            file.write("dataset_qz = %f \n" % AScans[i][j][2])  # momentum transfer
            file.write("dataset_R0 = %e \n" % AScans[i][j][3])  # reflectivity
            # file.write("dataset_eng = %f \n" % AScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):
                # General info of each scan
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(AInfo[i - 1][0]) + "-" + str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i - 1][
                    1] + "-" + AInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % AInfo[i][3])
                if (AInfo[i-1][1] == "S" or AInfo[i-1][1] =="P"):
                    file.write("polarization = AL \n")
                elif (AInfo[i-1][1] == "L" or AInfo[i-1][1] == "R"):
                    file.write("polarization = AC \n")
                file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

                # writing data
                for j in range(len(AScans[i][:, 0])):
                    file.write("dataset_qz = %f \n" % AScans[i][j][2])  # momentum transfer
                    file.write("dataset_A = %e \n" % (
                                (AScans[i - 1][j][3] - AScans[i][j][3]) / (AScans[i - 1][j][3] + AScans[i][j][3])))

                file.write("\n\n")
                dsNum = dsNum + 1

    # Writing energy scan data
    for i in range(len(EScans)):
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]), 2)) + "_Th" + str(round(float(EInfo[i][4]), 2)) + "_" + EInfo[i - 1][
                   1] + "-" + EInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % EInfo[i][3])
        if (EInfo[i][1] == "S"):
            file.write("polarization = S \n")
        elif (EInfo[i][1] == "P"):
            file.write("polarization = P \n")
        elif (EInfo[i][1] == "L"):
            file.write("polarization = LC \n")
        elif (EInfo[i][1] == "R"):
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(EScans[i][:, 0]))

        # Writing the data for each energy scan
        for j in range(len(EScans[i][:, 0])):
            file.write("dataset_qz = %f \n" % EScans[i][j][2])
            file.write("dataset_R0 = %e \n" % EScans[i][j][3])
            file.write("dataset_eng = %f \n" % EScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):

                # General information
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(EInfo[i - 1][0]) + "-" + str(EInfo[i][0]) + "_E" + str(
                    round(float(EInfo[i][3]), 2)) + "_Th" + str(round(float(EInfo[i][4]), 2)) + "_" + EInfo[i - 1][
                           1] + "-" + EInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % EInfo[i][3])
                if (EInfo[i-1][1] == "S" or EInfo[i-1][1] =="P"):
                    file.write("polarization = AL \n")
                elif (EInfo[i-1][1] == "L" or EInfo[i-1][1] == "R"):
                    file.write("polarization = AC \n")

                # Writing the data for each energy scan
                for j in range(len(EScans[i][:, 0])):
                    file.write("dataset_qz = %f \n" % EScans[i][j][2])  # momentum transfer
                    file.write("dataset_A = %e \n" % (
                                (EScans[i - 1][j][3] - EScans[i][j][3]) / (EScans[i - 1][j][3] + EScans[i][j][3])))
                    file.write("dataset_eng = %f \n" % EScans[i][j][0])

                file.write("\n\n")
                dsNum = dsNum + 1



def WriteSimulationASCII(file, AScans,AInfo,EScans,EInfo, sample):
    """
        Purpose: Write experimental data to .all file for program use
        :param fname: Name of data file
        :param AScans: Reflection scan experimental data
        :param AInfo: Info about reflection data scans
        :param EScans: Energy scan experimental data
        :param EInfo: Energy scan info
        :param header: Additional header
        :return:
        """

    file.write("# Simulation \n")  # header that states start of simulation data
    qz = list()

    # Writing reflectivity data
    dsNum = 1
    polarization = 'S'
    for i in range(len(AScans)):
        # General information
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % AInfo[i][3])
        if (AInfo[i][1] == "S"):
            polarization = 'S'
            file.write("polarization = S \n")
        elif (AInfo[i][1] == "P"):
            polarization = 'P'
            file.write("polarization = P \n")
        elif (AInfo[i][1] == "L"):
            polarization = 'LC'
            file.write("polarization = LC \n")
        elif (AInfo[i][1] == "R"):
            polarization = 'RC'
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

        qz = AScans[i][:,2]  # retrieves momentum transfer
        qz, R = sample.reflectivity(float(AInfo[i][3]), qz)  # computes the reflectivity
        R = R[polarization]  # retrieves the reflectivity of the correct polarization

        # Writes the simulated data
        for j in range(len(qz)):
            file.write("dataset_qz = %f \n" % qz[j])  # momentum transfer
            file.write("dataset_R0 = %e \n" % R[j])  # reflectivity

        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):
                # General information
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(AInfo[i - 1][0]) + "-" + str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i - 1][
                    1] + "-" + AInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % AInfo[i][3])
                if (AInfo[i - 1][1] == "S" or AInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                    file.write("polarization = AL \n")
                elif (AInfo[i - 1][1] == "L" or AInfo[i - 1][1] == "R"):
                    polarization = 'AC'
                    file.write("polarization = AC \n")
                file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

                qz = AScans[i][:,2]  # retrieves the momentum transfer
                qz, R = sample.reflectivity(float(AInfo[i][3]), qz)  # computes reflectivity
                R = R[polarization]  # retrieves the reflectivity of the correct polarization

                # Writes the simulated data
                for j in range(len(qz)):
                    file.write("dataset_qz = %f \n" % qz[j])
                    # print(AScans[i-1][j][3]+AScans[i][j][3])
                    file.write("dataset_A = %e \n" % R[j])
                    # file.write("dataset_eng = %f \n" % AScans[i][j][0])
                file.write("\n\n")
                dsNum = dsNum + 1

    # Writing the energy scan simulation data
    for i in range(len(EScans)):
        # General information
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]), 2)) + "_Th" + str(
            round(float(EInfo[i][4]), 2)) + "_" + EInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % EInfo[i][3])
        if (EInfo[i][1] == "S"):
            polarization = 'S'
            file.write("polarization = S \n")
        elif (EInfo[i][1] == "P"):
            polarization = 'P'
            file.write("polarization = P \n")
        elif (EInfo[i][1] == "L"):
            polarization = 'S'
            file.write("polarization = LC \n")
        elif (EInfo[i][1] == "R"):
            polarization = 'S'
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(EScans[i][:, 0]))

        E = EScans[i][:,0]  # numpy energy array
        E, R = sample.energy_scan(float(EInfo[i][4]), E)  # computes the energy scan
        R = R[polarization]  # retrieves the energy scan for the correct polarization

        # Writes the energy scan data
        for j in range(len(E)):
            file.write("dataset_R0 = %e \n" % R[j])  # reflectivity
            file.write("dataset_eng = %f \n" % E[j])  # energy
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
                #General information
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(EInfo[i - 1][0]) + "-" + str(EInfo[i][0]) + "_E" + str(
                    round(float(EInfo[i][3]), 2)) + "_Th" + str(round(float(EInfo[i][4]), 2)) + "_" + EInfo[i - 1][
                           1] + "-" + EInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % EInfo[i][3])
                if (EInfo[i - 1][1] == "S" or EInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                    file.write("polarization = AL \n")
                elif (EInfo[i - 1][1] == "L" or EInfo[i - 1][1] == "R"):
                    polarization = 'AC'
                    file.write("polarization = AC \n")

                E = EScans[i][:,0]  # numpy energy array
                E, R = sample.energy_scan(float(EInfo[i][4]), E)  # compute energy scan
                R = R[polarization]  # retrieve energy scan with correct polarization

                # write energy scan data
                for j in range(len(E)):
                    file.write("dataset_A = %e \n" % R[j])  # reflectivity
                    file.write("dataset_eng = %f \n" % E[j])  # energy

                file.write("\n\n")
                dsNum = dsNum + 1


def WriteDataASCII(fname,AScans,AInfo,EScans,EInfo, sample):
    """
    Purpose: Write experimental, simulated data, and sample info to an ASCII file with .all extension
    :param fname: File name
    :param AScans: Array that contains the reflectivity scan experimental data
    :param AInfo: Array that contains the information related to each reflecitivty scan
    :param EScans: Array that contains the energy scan experimental data
    :param EInfo: Array that contains the information related to each energy scan
    :param sample: Slab object containing sample information
    :return:
    """
    # If you decide to use your own functions to create an ASCII file you must include the following headers before
    # you write the data
    #  - # Structure
    #  - # Experimental_Data
    #  - # Simulation

    # Information for the data fitting and neural network parameters will be added later on

    file = open(fname, "w")  # creating a new file with name fname
    WriteSampleASCII(file, sample)  # writes the sample information
    WriteExperimentalDataASCII(file, AScans,AInfo,EScans,EInfo)  # writes the experimental data
    WriteSimulationASCII(file, AScans, AInfo, EScans, EInfo, sample)  # writes the simulation data
    f.close()

def getScanInfo(title):
    """
    Purpose: Retrieves important information in the scan title
    :param title: title/label of the scan
    :return:
    """
    title = title.split('_')  # split title into an array

    angle = None
    scan_number = title[0]

    # Determines the scan type
    if title[1] == 'A':
        scanType = 'Reflectivity'
    else:
        scanType = 'Energy'
        angle = title[2].replace('Th','')

    return scan_number, scanType, angle

def createNewDict():
    """
    Purpose: Initializes a dicitonary with the required keys and values
    :return:
    """
    my_dict = dict()
    my_dict['scanNumber'] = None
    my_dict['dataNumber'] = None
    my_dict['scanType'] = None
    my_dict['angle'] = None
    my_dict['energy'] = None
    my_dict['polarization'] = None
    my_dict['numberPoints'] = None
    return my_dict

def ConvertASCIItoHDF5(fascii, fhdf5):
    """
    Purpose: Converts and ASCII file to an hdf5 file with the proper format
    :param fname: File name of the ASCII file
    :return:
    Current implementation does not allow user to change the hdf5 file name from the ASCII file name
    """

    # checks to make sure file type of ASCII and HDF5 files are correct
    if not(fascii.endswith('.all')):
        raise NameError('File name must be a .all file type')
    if not(fhdf5.endswith('.h5')):
        raise NameError('File name must be a .all file type')

    Sinfo, Sscan, SimInfo, SimScan, sample = ReadDataASCII(fascii)  # retrieves information from ASCII file

    # checking to make sure that the HDF5 file asked for does not already exist
    cwd = os.getcwd()
    path = cwd + '/' + fhdf5
    if os.path.exists(path):
        raise OSError(
            "HDF5 file already exists. To write new HDF5 file remove the old file from the current working directory.")

    f = h5py.File(fhdf5, 'a')  # creates a new HDF5 file

    # Sample conversion of general information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)



    # general layer information
    dsLayer = 0
    for my_layer in sample.structure:

        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        # reconstructing the chemical formula input by user
        formula = ''
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)
            stoich = int(my_layer[ele].stoichiometry)
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # Magnetic element check
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        layer.attrs['Formula'] = formula
        dsLayer = dsLayer + 1

    # Experimental and simulated data conversion

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    name = ''

    dsNum = 1
    for i in range(len(Sinfo)):
        info = Sinfo[i]
        data = Sscan[i]
        simData = SimScan[i]

        scanType = info['scanType']

        # Reflectivity scans
        if scanType == "Reflectivity":
            qz = np.array(data[0])
            R0 = np.array(data[1])
            R0sim = np.array(simData[1])
            energy = info['energy']
            Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / np.pi  # initial angle
            dat = [qz, Theta, R0]
            sim = [qz, Theta, R0sim]

            polarization = info['polarization']

            if polarization == 'S' or polarization == 'P' or polarization == 'LC' or polarization == 'RC':
                name = str(info['dataNumber']) + "_" + str(np.round(energy,2)) + "_" + info['polarization']
            elif polarization == 'AL' or polarization == 'AC':
                name = str(info['dataNumber']) + "_" +str(np.round(energy,2)) + "_" + info['polarization'] +"_Asymm"

            datasetpoints = len(qz)

            dset = grpR.create_dataset(name, data=dat)
            dset1 = subR.create_dataset(name, data=sim)


            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = polarization
            dset1.attrs['Polarization'] = polarization

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)

            dsNum = dsNum + 1

        elif scanType == 'Energy':
            E = np.array(data[0])
            R = np.array(data[1])
            Rsim = np.array(simData[1])
            qz = np.array(data[2])
            energy = info['energy']
            angle = info['angle']
            Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / np.pi  # initial angle

            dat = [qz, Theta, R, E]
            sim = [qz, Theta, Rsim, E]

            polarization = info['polarization']

            if polarization == 'S' or polarization == 'P' or polarization == 'LC' or polarization == 'RC':
                name = str(info['dataNumber']) + "_E" +str(np.round(energy,2)) + "_Th" + str(np.round(angle,2)) + "_" + info['polarization']
            elif polarization == 'AL' or polarization == 'AC':
                name = str(info['dataNumber']) + "_E" +str(np.round(energy,2)) + "_Th" + str(np.round(angle,2)) + "_" + info['polarization'] + "_Asymm"

            datasetpoints = len(E)
            dset = grpE.create_dataset(name, data=dat)
            dset1 = subE.create_dataset(name, data=sim)

            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = polarization
            dset1.attrs['Polarization'] = polarization

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['Angle'] = float(angle)
            dset1.attrs['Angle'] = float(angle)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)
            dsNum = dsNum + 1


    f.close()


def ReadDataASCII(fname):
    """
    Purpose: Read .all file that contains sample and experimental information
    :param fname: Name of file contatining related information
    :return: Sscan - list of energy/reflectivity data
             Sinfo - contains a dictionary of relevant scan information
             sample - pre-initialized slab class
    """

    file = open(fname)
    Structure = False
    experimental_data = False
    simulation = False
    for line in file:
        line = line.split()
        if len(line) != 0:
            if line[0] == '#':
                if line[1] == 'Structure':
                    Structure = True
                    experimental_data = False
                    simulation = False
                    layer = 0

                    new_element = False
                    element = ''

                elif line[1] == 'Experimental_Data':

                    if simulation:
                        SimInfo[idx]['scanNumber'] = scan_number
                        SimInfo[idx]['scanType'] = scanType
                        SimInfo[idx]['angle'] = angle
                        SimInfo[idx]['energy'] = energy
                        SimInfo[idx]['polarization'] = polarization
                        SimInfo[idx]['numberPoints'] = numberPoints
                        if scanType == 'Energy':
                            SimScan.append([x_axis, y_axis, energy_qz])
                        elif scanType == 'Reflectivity':
                            SimScan.append([x_axis, y_axis])

                    Structure = False
                    experimental_data = True
                    simulation = False

                    idx = 0  # index of scan
                    Sinfo = []  # will contain information of sample
                    Sinfo.append(createNewDict())
                    Sscan = []

                    # initialization of parameters
                    x_axis = list()
                    y_axis = list()
                    energy_qz = list()
                    scan_number = 0
                    scanType = 0
                    angle = 0
                    energy = 0
                    polarization = 0
                    numberPoints = 0

                    #NewScan = False  # determines if we have a new scan
                    first = True
                    # Read in each line one at a time
                elif line[1] == 'Simulation':
                    # Attaches last experiment element
                    if experimental_data:
                        Sinfo[idx]['scanNumber'] = scan_number
                        Sinfo[idx]['scanType'] = scanType
                        Sinfo[idx]['angle'] = angle
                        Sinfo[idx]['energy'] = energy
                        Sinfo[idx]['polarization'] = polarization
                        Sinfo[idx]['numberPoints'] = numberPoints
                        if scanType == 'Energy':
                            Sscan.append([x_axis, y_axis, energy_qz])
                        elif scanType == 'Reflectivity':
                            Sscan.append([x_axis, y_axis])

                    Structure = False
                    experimental_data = False
                    simulation = True

                    idx = 0  # index of scan
                    SimInfo = []  # will contain information of sample
                    SimInfo.append(createNewDict())
                    SimScan = []

                    # initialization of parameters
                    x_axis = list()
                    y_axis = list()
                    energy_qz = list()
                    scan_number = 0
                    scanType = 0
                    angle = 0
                    energy = 0
                    polarization = 0
                    numberPoints = 0

                    first = True
                    # Read in each line one at a time
            else:
                if Structure:
                    if '=' not in line:
                        raise SyntaxError('Data file is improperly initialized.')
                    line.remove('=')  # removes the equal sign

                    # initializing the slab with correct number of layers
                    if line[0] == 'numberlayers':
                        sample = slab(int(line[1]))
                    elif line[0] == 'polyelements':
                        line.pop(0)
                        line = ''.join(line)
                        polyelements = ast.literal_eval(line)
                        sample.poly_elements = polyelements
                    elif line[0] == 'magelements':
                        line.pop(0)
                        line = ''.join(line)
                        magelements = ast.literal_eval(line)
                        sample.mag_elements = magelements

                    elif line[0] == 'layermagnetized':
                        line.pop(0)
                        line = ''.join(line)
                        layermagnetized = ast.literal_eval(line)
                        sample.layer_magnetized = layermagnetized
                    elif line[0] == 'scalingfactor':
                        scaling_factor = float(line[1])
                        sample.scaling_factor = scaling_factor
                    elif line[0] == 'backgroundshift':
                        background_shift = float(line[1])
                        sample.background_shift = background_shift
                    elif line[0] == 'layer':
                        layer = int(line[1])
                    elif line[0] == 'formula':
                        formula = line[1]

                        thickness = 20  # temporary assignment of the layer thickness
                        sample.addlayer(layer, formula, thickness)  # add layer

                    elif line[0] == 'element':
                        element = line[1]
                        new_element = True
                        sample.structure[layer][element].name = element

                    if new_element:
                        if line[0] == 'molarmass':
                            molarmass= float(line[1])
                            sample.structure[layer][element].molar_mass = molarmass

                        elif line[0] == 'density':
                            density = float(line[1])
                            sample.structure[layer][element].density = density
                        elif line[0] == 'thickness':
                            thickness = float(line[1])
                            sample.structure[layer][element].thickness = thickness
                        elif line[0] == 'roughness':
                            roughness = float(line[1])
                            sample.structure[layer][element].roughness = roughness
                        elif line[0] == 'linkedroughness':
                            linked_roughness = float(line[1])
                            sample.structure[layer][element].linked_roughness = linked_roughness
                        elif line[0] == 'scatteringfactor':
                            line.pop(0)
                            if len(line) == 1:
                                scatteringfactor = line[0]
                            else:
                                scatteringfactor = ast.literal_eval(''.join(line))
                            sample.structure[layer][element].scattering_factor = scatteringfactor

                        elif line[0] == 'polymorph':
                            line.pop(0)
                            polymorph = ast.literal_eval(''.join(line))
                            sample.structure[layer][element].polymorph = polymorph
                        elif line[0] == 'polyratio':
                            line.pop(0)
                            if len(line) == 1:
                                polyratio = 1
                            else:
                                polyratio = ast.literal_eval(''.join(line))
                                polyratio = np.array([float(poly) for poly in polyratio])
                            sample.structure[layer][element].poly_ratio = np.array(polyratio)
                        elif line[0] == 'gamma':
                            gamma = float(line[1])
                            sample.structure[layer][element].gamma = gamma
                        elif line[0] == 'phi':
                            phi = float(line[1])
                            sample.structure[layer][element].phi = phi
                        elif line[0] == 'magdensity':
                            line.pop(0)
                            magdensity = ast.literal_eval(''.join(line))
                            if len(magdensity) != 0:
                                magdensity = np.array([float(mag) for mag in magdensity])
                            sample.structure[layer][element].mag_density = magdensity

                        elif line[0] == 'magscatteringfactor':
                            line.pop(0)
                            if len(line)==0:
                                magscatteringfactor = None
                            else:
                                magscatteringfactor = ast.literal_eval(''.join(line))
                            sample.structure[layer][element].mag_scattering_factor = magscatteringfactor
                        elif line[0] == 'position':
                            position = int(line[1])
                            sample.structure[layer][element].position = position
                            new_element = False  # make sure we do not enter this if statement

                elif experimental_data:

                    if '=' not in line:
                        raise SyntaxError('Data file is improperly initialized.')
                    line.remove('=')  # removes the equal sign

                    info = line[0]  # data identifier
                    data = line[1]  # data value

                    # retrieves the datasetnumber
                    if info == 'datasetnumber':
                        if first:
                            data = int(data)
                            Sinfo[idx]['dataNumber'] = data
                            first = False
                        else:
                            Sinfo[idx]['scanNumber'] = scan_number
                            Sinfo[idx]['scanType'] = scanType
                            Sinfo[idx]['angle'] = angle
                            Sinfo[idx]['energy'] = energy
                            Sinfo[idx]['polarization'] = polarization
                            Sinfo[idx]['numberPoints'] = numberPoints

                            # Attaches appropriate dataset for different scan types
                            if scanType == 'Energy':
                                Sscan.append([x_axis, y_axis, energy_qz])
                            elif scanType == 'Reflectivity':
                                Sscan.append([x_axis, y_axis])


                            # resets all parameters when new scan
                            Sinfo.append(createNewDict())
                            x_axis = list()
                            y_axis = list()
                            energy_qz = list()
                            scan_number = None
                            scanType = None
                            angle = None
                            energy = None
                            numberPoints = None
                            idx = idx + 1
                            data = int(data)
                            Sinfo[idx]['dataNumber'] = data

                    # retrieves the data set title
                    if info == 'datasettitle':
                        scan_number, scanType, angle = getScanInfo(data)
                        if angle != None:
                            angle = float(angle)

                    # sets parameters based on scan type
                    if scanType == 'Energy':
                        if info == 'datasetenergy':
                            energy = float(data)
                        if info == 'polarization':
                            polarization = data
                        if info == 'datasetpoints':
                            numberPoints = int(data)
                        if info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_eng':
                            data = float(data)
                            x_axis.append(data)
                        if info == 'dataset_qz':
                            data = float(data)
                            energy_qz.append(data)
                    elif scanType == 'Reflectivity':
                        if info == 'datasetenergy':
                            energy = float(data)
                        elif info == 'polarization':
                            polarization = data
                        elif info == 'datasetpoints':
                            numberPoints = int(data)
                        elif info == 'dataset_qz':
                            data = float(data)
                            x_axis.append(data)
                        elif info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        elif info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)

                elif simulation:

                    if '=' not in line:
                        raise SyntaxError('Data file is improperly initialized.')
                    line.remove('=')  # removes the equal sign

                    info = line[0]  # data identifier
                    data = line[1]  # data value

                    # retrieves the datasetnumber
                    if info == 'datasetnumber':

                        if first:
                            data = int(data)
                            SimInfo[idx]['dataNumber'] = data
                            first = False
                        else:
                            SimInfo[idx]['scanNumber'] = scan_number
                            SimInfo[idx]['scanType'] = scanType
                            SimInfo[idx]['angle'] = angle
                            SimInfo[idx]['energy'] = energy
                            SimInfo[idx]['polarization'] = polarization
                            SimInfo[idx]['numberPoints'] = numberPoints

                            # Attaches different elements for different scan types
                            if scanType == 'Energy':
                                SimScan.append([x_axis, y_axis, energy_qz])
                            elif scanType == 'Reflectivity':
                                SimScan.append([x_axis, y_axis])


                            # resets all parameters when new scan
                            SimInfo.append(createNewDict())
                            x_axis = list()
                            y_axis = list()
                            energy_qz = list()
                            scan_number = None
                            scanType = None
                            angle = None
                            energy = None
                            numberPoints = None
                            idx = idx + 1
                            data = int(data)
                            SimInfo[idx]['dataNumber'] = data

                    # retrieves the data set title
                    if info == 'datasettitle':
                        scan_number, scanType, angle = getScanInfo(data)
                        if angle != None:
                            angle = float(angle)

                    # sets parameters based on scan type
                    if scanType == 'Energy':
                        if info == 'datasetenergy':
                            energy = float(data)
                        if info == 'polarization':
                            polarization = data
                        if info == 'datasetpoints':
                            numberPoints = int(data)
                        if info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_eng':
                            data = float(data)
                            x_axis.append(data)
                        if info == 'dataset_qz':
                            data = float(data)
                            energy_qz.append(data)
                    elif scanType == 'Reflectivity':
                        if info == 'datasetenergy':
                            energy = float(data)
                        elif info == 'polarization':
                            polarization = data
                        elif info == 'datasetpoints':
                            numberPoints = int(data)
                        elif info == 'dataset_qz':
                            data = float(data)
                            x_axis.append(data)
                        elif info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        elif info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)


    if experimental_data:
        Sinfo[idx]['scanNumber'] = scan_number
        Sinfo[idx]['scanType'] = scanType
        Sinfo[idx]['angle'] = angle
        Sinfo[idx]['energy'] = energy
        Sinfo[idx]['polarization'] = polarization
        Sinfo[idx]['numberPoints'] = numberPoints
        if scanType == 'Energy':
            Sscan.append([x_axis, y_axis, energy_qz])
        elif scanType == 'Reflectivity':
            Sscan.append([x_axis, y_axis])

    elif simulation:
        SimInfo[idx]['scanNumber'] = scan_number
        SimInfo[idx]['scanType'] = scanType
        SimInfo[idx]['angle'] = angle
        SimInfo[idx]['energy'] = energy
        SimInfo[idx]['polarization'] = polarization
        SimInfo[idx]['numberPoints'] = numberPoints
        if scanType == 'Energy':
            SimScan.append([x_axis, y_axis, energy_qz])
        elif scanType == 'Reflectivity':
            SimScan.append([x_axis, y_axis])

    f.close()

    return Sinfo, Sscan, SimInfo, SimScan, sample



def hdf5ToDict(hform):

    mydict = dict()
    for key in hform.keys():
        mydict[key] = dict()
        mydict[key]['Data'] = list(hform[key])

        for attrskey,val in hform[key].attrs.items():
            mydict[key][attrskey] = val

    return mydict

def saveNewFile(fname, info, data_dict, sample):
    """
        Purpose: Take in data and sample model and save it as an hdf5 file
        :param fname: File name with .hdf5 file extension
        :param AScans: Reflectivity scan data from ProcessRXR
        :param AInfo: Reflectivity scan information from ProcessRXR
        :param EScans: Energy scan data from Process RXR
        :param EInfo: Energy scan information from Process RXR
        :param sample: sample object
        :return:
        """

    # Checking if fname already exists
    cwd = os.getcwd()
    path = cwd + '/' + fname
    if os.path.exists(path):
        raise OSError(
            "HDF5 file already exists. To write new file remove the old file from the current working directory.")

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # creating group that will contain the sample information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)
    grp1.attrs['ScalingFactor'] = float(sample.scaling_factor)
    grp1.attrs['BackgroundShift'] = float(sample.background_shift)

    dsLayer = 0
    for my_layer in sample.structure:

        # Layer information
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            # Element information
            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # May be changed in the future as layermagnetized also contains this information
            # Original implemented to avoid problem of trying to load in the magnetic data that does not exist
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    # Loading in the experimental data and simulated data
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    for inf in info:
        name = inf[2]

        if inf[1] == 'Energy':
            energy = data_dict[name]['Energy']
            angle = data_dict[name]['Angle']
            polarization = data_dict[name]['Polarization']
            datasetpoints = data_dict[name]['DataPoints']
            dsNum = data_dict[name]['DatasetNumber']
            data = data_dict[name]['Data']
            theta = data[1]
            E = data[3]

            E, R = sample.energy_scan(angle, E)
            R = R[polarization]

            sim = np.array([data[0], data[1], R, data[3]])
            dat = np.array([data[0], data[1], data[2], data[3]])

            dset = grpE.create_dataset(name, data=dat)
            dset1 = subE.create_dataset(name, data=sim)

            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = polarization
            dset1.attrs['Polarization'] = polarization

            dset.attrs['Angle'] = float(angle)
            dset1.attrs['Angle'] = float(angle)

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)

        elif inf[1] == 'Reflectivity':

            energy = data_dict[name]['Energy']
            polarization = data_dict[name]['Polarization']
            datasetpoints = data_dict[name]['DataPoints']
            dsNum = data_dict[name]['DatasetNumber']
            data = data_dict[name]['Data']
            dat = np.array([data[0], data[1], data[2]])
            qz, R = sample.reflectivity(energy, data[0])
            R = R[polarization]

            sim = np.array([data[0],data[1],R])

            dset = grpR.create_dataset(name, data=dat)
            dset1 = subR.create_dataset(name, data=sim)

            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = str(polarization)
            dset1.attrs['Polarization'] = str(polarization)

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)

    f.close()
    return

def QUAD_to_data_dict(fname):
    data_dict = dict()

    # read the data file and convert each line into a list format
    with open(fname) as f:
        lines = [line for line in f]

    f.close()
    start_new = True
    data = [[],[],[],[]]
    name = ''
    scan_type = ''
    for idx,line in enumerate(lines):
        if idx == 0 or start_new:  # initialization
            line = line.split()
            scan_number = line[0]
            pol = line[1]
            scan_type = line[2]
            energy = float(line[3])
            angle = float(line[4])

            temp_energy = str(round(float(line[3]),2))
            temp_angle = str(round(float(line[4]),2))

            name = ''
            if scan_type == 'A':
                name = line[0] + '_' + temp_energy + '_' + line[1]
            else:
                name = line[0] + '_E' + temp_energy + '_Th' + temp_angle + '_' + line[1]

            data_dict[name] = dict()
            data_dict[name]['Background Shift'] = 0
            data_dict[name]['Scaling Factor'] = 1
            data_dict[name]['DatasetNumber'] = scan_number
            data_dict[name]['Polarization'] = pol
            data_dict[name]['Energy'] = energy
            if scan_type == 'E':
                data_dict[name]['Angle'] = angle

            start_new = False
            data = [[],[],[],[]]

        elif line[0] == '=':
            # terminate and start new
            n = len(data[0])
            if scan_type == 'A':
                data[0] = np.array(data[0])
                data[1] = np.array(data[1])
                data[2] = np.array(data[2])
                data.pop()
            else:
                data[0] = np.array(data[0])
                data[1] = np.array(data[1])
                data[2] = np.array(data[2])
                data[3] = np.array(data[3])

            data_dict[name]['DataPoints'] = n
            data_dict[name]['Data'] = data
            start_new = True
        else:
            line = line.split()
            # save the data
            if scan_type == 'A':
                data[0].append(float(line[2]))  # qz
                data[1].append(float(line[1]))  # Theta
                data[2].append(float(line[3]))  # R
            else:
                data[0].append(float(line[2]))  # qz
                data[1].append(float(line[1]))  # Theta
                data[2].append(float(line[3]))  # R
                data[3].append(float(line[0]))  # E
    return data_dict

def createDataHDF5fromDict(filename, data_dict_list):
    """
    Purpose: Create HDF5 data file from data_dict
    :param filename: Name of file (must have .h5 extension)
    :param data_dict_list: Dictionary or list of dictionaries
    :return: N/A
    """

    if filename.endswith('.h5'):

        f = h5py.File(filename, "a")
        if type(data_dict_list) is dict:
            data_dict_list = [data_dict_list]

        simulated = f.create_group('Simulated_data')
        simR = simulated.create_group('Reflectivity_Scan')
        simE = simulated.create_group('Energy_Scan')

        experiment = f.create_group('Experimental_data')
        reflScan = experiment.create_group('Reflectivity_Scan')
        energyScan = experiment.create_group('Energy_Scan')

        for data_dict in data_dict_list:  # loops through all data_dicts
            for name in data_dict.keys(): # loop through all data
                if 'Angle' in data_dict[name].keys():
                    dat = data_dict[name]['Data']
                    dat = np.array(dat)
                    m = np.shape(dat)
                    dset = energyScan.create_dataset(name, m, data=dat, maxshape=(4, None), chunks=True)
                    dset.attrs['DatasetNumber'] = data_dict[name]['DatasetNumber']
                    dset.attrs['DataPoints'] = data_dict[name]['DataPoints']
                    dset.attrs['Energy'] = data_dict[name]['Energy']
                    dset.attrs['Angle'] = data_dict[name]['Angle']
                    dset.attrs['Polarization'] = data_dict[name]['Polarization']
                    dset.attrs['Background Shift'] = data_dict[name]['Background Shift']
                    dset.attrs['Scaling Factor'] = data_dict[name]['Scaling Factor']

                    dat1 = data_dict[name]['Data']
                    dset1 = simE.create_dataset(name, m, data=dat1, maxshape=(4, None), chunks=True)
                    dset1.attrs['DatasetNumber'] = data_dict[name]['DatasetNumber']
                    dset1.attrs['DataPoints'] = data_dict[name]['DataPoints']
                    dset1.attrs['Energy'] = data_dict[name]['Energy']
                    dset1.attrs['Angle'] = data_dict[name]['Angle']
                    dset1.attrs['Polarization'] = data_dict[name]['Polarization']
                    dset1.attrs['Background Shift'] = data_dict[name]['Background Shift']
                    dset1.attrs['Scaling Factor'] = data_dict[name]['Scaling Factor']
                else:
                    dat = data_dict[name]['Data']
                    dat = np.array(dat)
                    m = np.shape(dat)

                    dset = reflScan.create_dataset(name, m, data=dat, maxshape=(3, None), chunks=True)
                    dset.attrs['DatasetNumber'] = data_dict[name]['DatasetNumber']
                    dset.attrs['DataPoints'] = data_dict[name]['DataPoints']
                    dset.attrs['Energy'] = data_dict[name]['Energy']
                    dset.attrs['Polarization'] = data_dict[name]['Polarization']
                    dset.attrs['Background Shift'] = data_dict[name]['Background Shift']
                    dset.attrs['Scaling Factor'] = data_dict[name]['Scaling Factor']

                    dat1 = data_dict[name]['Data']
                    dset1 = simR.create_dataset(name, m, data=dat1, maxshape=(3, None), chunks=True)
                    dset1.attrs['DatasetNumber'] = data_dict[name]['DatasetNumber']
                    dset1.attrs['DataPoints'] = data_dict[name]['DataPoints']
                    dset1.attrs['Energy'] = data_dict[name]['Energy']
                    dset1.attrs['Polarization'] = data_dict[name]['Polarization']
                    dset1.attrs['Background Shift'] = data_dict[name]['Background Shift']
                    dset1.attrs['Scaling Factor'] = data_dict[name]['Scaling Factor']

        f.close()


class DataFile:
    # This class is used to save the data files in the appropriate form for GO-RXR
    def __init__(self):
        self.experimental_data = dict()

    def addEnergyScan(self, dnum, min_energy, theta, polarization, EList, RList, thetaList=[], background_shift=0, scaling_factor=1):
        """
        Purpose: Add energy scan to experimental data
        :param dnum: data scan number (entered as an integer type)
        :param min_energy: Minimum energy (entered as a float)
        :param theta: Constant grazing angle in degrees (entered as a float/integer)
        :param EList: List of energies (energies must be float values)
        :param RList: List of reflectivity values (reflectivity must be float values)
        :param thetaList: List of measured theta values (default is to use constant theta value)
        :param background_shift: Shift value to apply to experimental data (default set to 0)
        :param scaling factor: Scaling factor applied to experimental data (default set to 1)
        """
        qzList = []
        # checks all variables that they are the correct type
        if type(dnum) != int:
            raise TypeError('Dataset number must be an integer type.')
        if type(min_energy) != float and type(min_energy) != int:
            raise TypeError('Minimum energy must be an integer or float type.')
        if type(theta) != float and type(theta) != int:
            raise TypeError('Grazing angle must be a float or integer type.')
        if polarization not in ['S', 'P', 'LC', 'RC']:
            raise TypeError('Polarization must be S, P, LC, or RC.')
        if type(EList) != list and type(EList) != np.ndarray:
            raise TypeError('List of energies must be a list or numpy array.')
        if type(RList) != list and type(RList) != np.ndarray:
            raise TypeError('List of reflectivity must be a list or numpy array.')
        if type(thetaList) != list and type(thetaList) != np.ndarray:
            raise TypeError('List of grazing angles must be a list or numpy array.')
        if type(background_shift) != float and type(background_shift) != int:
            raise TypeError('Background shift must be an integer or float type.')
        if type(scaling_factor) != float and type(scaling_factor) != int:
            raise TypeError('Scaling factor must be an integer or float type.')


        if len(thetaList) == 0:
            if len(EList) != len(RList):
                raise ValueError('Elist and Rlist must be the same length')
            thetaList = [theta for i in range(len(EList))]  # creates thetaList if not set
            qzList = np.sin(np.array(thetaList) * np.pi / 180) * (np.array(EList) * 0.001013546143)
        else:
            if len(EList) != len(RList) and len(EList) != len(thetaList):
                raise ValueError('Elist, Rlist, and thetaList must be the same length')
            print(thetaList)
            qzList = np.sin(np.array(thetaList) * np.pi / 180) * (np.array(EList) * 0.001013546143)

        num_points = len(EList)
        # create energy scan name
        temp_energy = str(round(float(min_energy), 2))
        temp_theta = str(round(float(theta), 2))
        name = str(dnum) + '_E' + str(temp_energy) + '_Th' + str(temp_theta) + '_' + polarization

        # create dictionary for input data
        self.experimental_data[name] = dict()
        self.experimental_data[name]['Background Shift'] = background_shift
        self.experimental_data[name]['Scaling Factor'] = scaling_factor
        self.experimental_data[name]['DatasetNumber'] = str(dnum)
        self.experimental_data[name]['Polarization'] = polarization
        self.experimental_data[name]['Energy'] = min_energy
        self.experimental_data[name]['Angle'] = theta
        self.experimental_data[name]['DataPoints'] = num_points

        # set the data in the appropriate form
        data = [np.array(qzList), np.array(thetaList), np.array(RList), np.array(EList)]
        self.experimental_data[name]['Data'] = data

    def addReflectivityScan(self,dnum, energy, polarization, qzList, RList, thetaList=[], background_shift=0,
                      scaling_factor=1):
        """
        Purpose: Add energy scan to experimental data
        :param dnum: data scan number (entered as an integer type)
        :param energy: Minimum energy (entered as a float)
        :param qzList: List of energies (energies must be float values)
        :param RList: List of reflectivity values (reflectivity must be float values)
        :param thetaList: List of measured theta values (default is to use constant theta value)
        :param background_shift: Shift value to apply to experimental data (default set to 0)
        :param scaling factor: Scaling factor applied to experimental data (default set to 1)
        """

        # checks all variables that they are the correct type
        if type(dnum) != int:
            raise TypeError('Dataset number must be an integer type.')
        if type(energy) != float and type(energy) != int:
            raise TypeError('Constant energy must be an integer or float type.')
        if polarization not in ['S', 'P', 'LC', 'RC']:
            raise TypeError('Polarization must be S, P, LC, or RC.')
        if type(qzList) != list and type(qzList) != np.ndarray:
            raise TypeError('List of energies must be a list or numpy array.')
        if type(RList) != list and type(RList) != np.ndarray:
            raise TypeError('List of reflectivity must be a list or numpy array.')
        if type(thetaList) != list and type(thetaList) != np.ndarray:
            raise TypeError('List of grazing angles must be a list or numpy array.')
        if type(background_shift) != float and type(background_shift) != int:
            raise TypeError('Background shift must be an integer or float type.')
        if type(scaling_factor) != float and type(scaling_factor) != int:
            raise TypeError('Scaling factor must be an integer or float type.')

        if len(thetaList) == 0:
            if len(qzList) != len(RList):
                raise ValueError('Elist and Rlist must be the same length')
            thetaList = np.arcsin(np.array(qzList)/(energy * 0.001013546143))*180/np.pi
        else:
            if len(qzList) != len(RList) and len(thetaList) != len(thetaList):
                raise ValueError('Elist, Rlist, and thetaList must be the same length')


        num_points = len(qzList)
        # create energy scan name
        temp_energy = str(round(float(energy), 2))

        name = str(dnum) + '_' + str(temp_energy) + '_' + polarization

        # create dictionary for input data
        self.experimental_data[name] = dict()
        self.experimental_data[name]['Background Shift'] = background_shift
        self.experimental_data[name]['Scaling Factor'] = scaling_factor
        self.experimental_data[name]['DatasetNumber'] = str(dnum)
        self.experimental_data[name]['Polarization'] = polarization
        self.experimental_data[name]['Energy'] = energy
        self.experimental_data[name]['DataPoints'] = num_points

        # set the data in the appropriate form
        data = [np.array(qzList), np.array(thetaList), np.array(RList)]
        self.experimental_data[name]['Data'] = data

    def addAsymmetrySpectrum(self, dnum, min_energy, theta, polarization, EList, R1List, R2List, thetaList=[], background_shift=0,
                      scaling_factor=1):
        """
        Purpose: Add energy scan to experimental data
        :param dnum: data scan number (entered as an integer type)
        :param min_energy: Minimum energy (entered as a float)
        :param theta: Constant grazing angle in degrees (entered as a float/integer)
        :param EList: List of energies (energies must be float values)
        :param R1List: List of reflectivity values (reflectivity must be float values)
        :param R2List: List of reflectivity values (reflectivity must be float values)
        :param thetaList: List of measured theta values (default is to use constant theta value)
        :param background_shift: Shift value to apply to experimental data (default set to 0)
        :param scaling factor: Scaling factor applied to experimental data (default set to 1)
        """
        # Note: For the asymmetry calculation
        #       -> A = (R1-R2)/(R1+R2)
        qzList = []
        # checks all variables that they are the correct type
        if type(dnum) != int:
            raise TypeError('Dataset number must be an integer type.')
        if type(min_energy) != float and type(min_energy) != int:
            raise TypeError('Minimum energy must be an integer or float type.')
        if type(theta) != float and type(theta) != int:
            raise TypeError('Grazing angle must be a float or integer type.')
        if polarization not in ['AC', 'AL']:
            raise TypeError('Polarization must be AL or AC.')
        if type(EList) != list and type(EList) != np.ndarray:
            raise TypeError('List of energies must be a list or numpy array.')
        if type(R1List) != list and type(R1List) != np.ndarray:
            raise TypeError('List of reflectivity must be a list or numpy array.')
        if type(R2List) != list and type(R2List) != np.ndarray:
            raise TypeError('List of reflectivity must be a list or numpy array.')
        if type(thetaList) != list and type(thetaList) != np.ndarray:
            raise TypeError('List of grazing angles must be a list or numpy array.')
        if type(background_shift) != float and type(background_shift) != int:
            raise TypeError('Background shift must be an integer or float type.')
        if type(scaling_factor) != float and type(scaling_factor) != int:
            raise TypeError('Scaling factor must be an integer or float type.')

        if len(thetaList) == 0:
            if len(EList) != len(R1List):
                raise ValueError('Elist and R1list must be the same length')
            if len(R2List) != len(R1List):
                raise ValueError('R1list and R2list must be the same length')
            thetaList = [theta for i in range(len(EList))]  # creates thetaList if not set
            qzList = np.sin(np.array(thetaList) * np.pi / 180) * (np.array(EList) * 0.001013546143)
        else:
            if len(EList) != len(R1List) and len(EList) != len(thetaList) and len(R1List) != len(R2List):
                raise ValueError('Elist, R1list, R2List, and thetaList must be the same length')

            qzList = np.sin(np.array(thetaList) * np.pi / 180) * (np.array(EList) * 0.001013546143)

        num_points = len(EList)
        # create energy scan name
        temp_energy = str(round(float(min_energy), 2))
        temp_theta = str(round(float(theta), 2))
        name = str(dnum) + '_E' + str(temp_energy) + '_Th' + str(temp_theta) + '_' + polarization

        # create dictionary for input data
        self.experimental_data[name] = dict()
        self.experimental_data[name]['Background Shift'] = background_shift
        self.experimental_data[name]['Scaling Factor'] = scaling_factor
        self.experimental_data[name]['DatasetNumber'] = str(dnum)
        self.experimental_data[name]['Polarization'] = polarization
        self.experimental_data[name]['Energy'] = min_energy
        self.experimental_data[name]['Angle'] = theta
        self.experimental_data[name]['DataPoints'] = num_points

        # set the data in the appropriate form
        A = (np.array(R1List)-np.array(R2List))/(np.array(R1List)+np.array(R2List))
        data = [np.array(qzList), np.array(thetaList), A, np.array(EList)]
        self.experimental_data[name]['Data'] = data

    def addAsymmetryCurve(self, dnum, energy, polarization, qzList, R1List, R2List, thetaList=[], background_shift=0,
                            scaling_factor=1):
        """
        Purpose: Add energy scan to experimental data
        :param dnum: data scan number (entered as an integer type)
        :param energy: Minimum energy (entered as a float)
        :param qzList: List of energies (energies must be float values)
        :param R1List: List of reflectivity values (reflectivity must be float values)
        :param R2List: List of reflectivity values (reflectivity must be float values)
        :param thetaList: List of measured theta values (default is to use constant theta value)
        :param background_shift: Shift value to apply to experimental data (default set to 0)
        :param scaling factor: Scaling factor applied to experimental data (default set to 1)
        """
        # Note: The asymmetry calculation is calculated as follows
        #           -> A = (R1-R2)/(R1+R2)

        # checks all variables that they are the correct type
        if type(dnum) != int:
            raise TypeError('Dataset number must be an integer type.')
        if type(energy) != float and type(energy) != int:
            raise TypeError('Constant energy must be an integer or float type.')
        if polarization not in ['AL', 'AC']:
            raise TypeError('Polarization must be AL or AC.')
        if type(qzList) != list and type(qzList) != np.ndarray:
            raise TypeError('List of energies must be a list or numpy array.')
        if type(R1List) != list and type(R1List) != np.ndarray:
            raise TypeError('List of reflectivity must be a list or numpy array.')
        if type(R2List) != list and type(R2List) != np.ndarray:
            raise TypeError('List of reflectivity must be a list or numpy array.')
        if type(thetaList) != list and type(thetaList) != np.ndarray:
            raise TypeError('List of grazing angles must be a list or numpy array.')
        if type(background_shift) != float and type(background_shift) != int:
            raise TypeError('Background shift must be an integer or float type.')
        if type(scaling_factor) != float and type(scaling_factor) != int:
            raise TypeError('Scaling factor must be an integer or float type.')

        if len(thetaList) == 0:
            if len(qzList) != len(R1List):
                raise ValueError('qzList and R1List must be the same length')
            if len(R2List) != len(R1List):
                raise ValueError('R1List and R2List must be the same length')
            thetaList = np.arcsin(np.array(qzList) / (energy * 0.001013546143)) * 180 / np.pi
        else:
            if len(qzList) != len(R1List) and len(thetaList) != len(thetaList) and len(R2List) != len(R1List):
                raise ValueError('qzList, R1List, R2List, and thetaList must be the same length')

        num_points = len(qzList)
        # create energy scan name
        temp_energy = str(round(float(energy), 2))

        name = str(dnum) + '_' + str(temp_energy) + '_' + polarization

        # create dictionary for input data
        self.experimental_data[name] = dict()
        self.experimental_data[name]['Background Shift'] = background_shift
        self.experimental_data[name]['Scaling Factor'] = scaling_factor
        self.experimental_data[name]['DatasetNumber'] = str(dnum)
        self.experimental_data[name]['Polarization'] = polarization
        self.experimental_data[name]['Energy'] = energy
        self.experimental_data[name]['DataPoints'] = num_points

        # set the data in the appropriate form
        A = (np.array(R1List)- np.array(R2List))/(np.array(R1List) + np.array(R2List))
        data = [np.array(qzList), np.array(thetaList), A]
        self.experimental_data[name]['Data'] = data


    def return_data_dict(self):
        # returns the data dictonary
        return self.experimental_data

    def is_name_in_data(self, name):
        """
        Purpose: determine if some experimental data has already been included into the dataset
        :param name: The name of the experimental data
        :return: True if found in experimental data
        """
        if name in list(self.experimental_data.keys()):
            return True
        else:
            return False

    def saveData(self,fname):
        createDataHDF5fromDict(fname, self.experimental_data)

if __name__ == "__main__":

    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density=5.12, roughness=4)
    sample.addlayer(1, 'SrTiO3', 5.28602, density=5.12, roughness=2)

    sample.addlayer(2, 'LaMnO3', 18.84,density = 5, roughness=3.77)
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [0.9,0.1], sf =['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0.00023,0], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 7.73, roughness= 4.54)
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.9, 0.1], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.001, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 3.51, roughness=0.788804)
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [0.9, 0.1], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.024, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 2.43, roughness=1.3591)
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 3.86, roughness=0.8215)
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.5, 0.5], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.005, 0], ['Co', 'Ni'])

    sample.addlayer(7, 'CCO', 4, density = 2.5, roughness = 2)

    fname = "Pim10uc.h5"
    fnew = 'test.h5'
    #info, data_dict, sim_dict=ReadDataHDF5(fname)
    #print(len(data_dict['59_E429.58_Th5.0_S']['Data'][3]))

    #WriteSampleHDF5(fname, sample)

    #Read_ReMagX("Pim10uc.all")


