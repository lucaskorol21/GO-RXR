from scipy import optimize, signal
from material_structure import *
import numpy as np
from data_structure import *
import matplotlib.pyplot as plt
import time
from tkinter import *
from material_model import *
from tkinter import ttk
import multiprocessing as mp
import sys
import os
from PIL import Image, ImageTk
import functools
import copy

global x_vars
x_vars = []

def changeSampleParams(x, parameters, sample, backS, scaleF, script, use_script=False):
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

        if params[0] == "SCALING FACTOR":  # scale the reflectivity spectra
            name = params[1]
            if name == 'ALL SCANS':  # case where all scans have same scaling factor
                for key in list(scaleF.keys()):
                    scaleF[key] = str(x[p])
            else:  # each scan has separate scaling factor
                scaleF[name] = str(x[p])

            sample.scaling_factor = x[p]  # sets the sample scaling factor
        elif params[0] == "BACKGROUND SHIFT":
            name = params[1]
            if name == 'ALL SCANS':
                for key in list(backS.keys()):
                    backS[key] = "{:e}".format(x[p])
            else:
                backS[name] = "{:e}".format(x[p])

        elif params[0] == 'SCATTERING FACTOR':
            mode =params[1]
            element = params[2]
            dE = x[p]

            if mode == 'STRUCTURAL':
                sample.eShift[element] = dE
            elif mode == 'MAGNETIC':
                sample.mag_eShift[element] = dE

        else:
            layer = params[0]  # Determine which layer to change
            property = params[1]  # structural/polymorphous/magnetic
            if property == 'STRUCTURAL':
                mode = params[2]
                # Determine if user wants to use compound or element mode
                if mode == 'COMPOUND':
                    characteristic = params[3]  # thickness/density/roughness/linked roughness
                    ele_idx = params[4]

                    my_ele = list(sample.structure[layer].keys())[ele_idx]
                    diffT = x[p] - sample.structure[layer][my_ele].thickness
                    diffD = x[p] - sample.structure[layer][my_ele].density
                    diffR = x[p] - sample.structure[layer][my_ele].roughness
                    diffLR = x[p] - sample.structure[layer][my_ele].linked_roughness
                    #print(sample.structure)
                    # determine the difference parameter for compound mode (will always be the first element)
                    for ele in list(sample.structure[layer].keys()):
                        if characteristic == 'THICKNESS':

                            if ele == my_ele:
                                sample.structure[layer][ele].thickness = x[p]
                            else:
                                sample.structure[layer][ele].thickness = sample.structure[layer][ele].thickness + diffT
                        elif characteristic == 'DENSITY':
                            stoich = sample.structure[layer][ele].stoichiometry  # stoichiometry
                            molar_mass = sample.structure[layer][ele].molar_mass  # molar mass
                            if ele == my_ele:
                                sample.structure[layer][ele].density = x[p]  # set density
                            else:
                                sample.structure[layer][ele].density = sample.structure[layer][ele].density + diffD*stoich # set density
                        elif characteristic == 'ROUGHNESS':
                            if ele == my_ele:
                                sample.structure[layer][ele].roughness = x[p]
                            else:
                                sample.structure[layer][ele].roughness = sample.structure[layer][ele].roughness + diffR
                        elif characteristic == 'LINKED ROUGHNESS':
                            if ele == my_ele:
                                sample.structure[layer][ele].linked_roughness = x[p]
                            else:
                                sample.structure[layer][ele].linked_roughness = sample.structure[layer][ele].linked_roughness + diffLR

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
                element = params[2]  # determines the element that contains the polymorph
                polymorph = params[3]  # determines the polymorph to change

                ratio = 1 - x[p]  # Assumes only two possible polymorphs for now and computes other polymorph ratio

                #poly = np.where(sample.structure[layer][element].polymorph == polymorph)  # determines location of polymorph
                my_poly = list(sample.structure[layer][element].polymorph)
                num = len(my_poly)

                poly = my_poly.index(polymorph)

                if num == 2:
                    # sets poly_ratio value making sure sum equals to 1
                    if poly == 0:
                        sample.structure[layer][element].poly_ratio[0] = x[p]
                        sample.structure[layer][element].poly_ratio[1] = ratio
                    elif poly == 1:
                        sample.structure[layer][element].poly_ratio[1] = x[p]
                        sample.structure[layer][element].poly_ratio[0] = ratio
                else:  # more than one element variation

                    # maintain the element variation ratios (except the change variable)
                    x_temp = []
                    sum = 0

                    my_idx = [i for i in range(num) if i != poly]  #indices of polymorph that are not the varying element variation
                    poly_ratio = copy.copy(np.array([float(sample.structure[layer][element].poly_ratio[i]) for i in range(num)]))

                    sum = np.sum(poly_ratio[my_idx])
                    #for i in range(num):
                    #    x_temp.append(float(sample.structure[layer][element].poly_ratio[i]))
                    #    if i != poly:
                    #        sum = sum + float(sample.structure[layer][element].poly_ratio[i])

                    x_prime = ratio*poly_ratio/sum
                    for i in range(num):
                        if i == poly:
                            sample.structure[layer][element].poly_ratio[i] = x[p]
                        else:
                            sample.structure[layer][element].poly_ratio[i] = x_prime[i]


            elif property == 'MAGNETIC':

                element = params[2]  # determines the magnetic element to use

                # determines if magnetic element is polymorphous
                if len(params) == 3:
                    sample.structure[layer][element].mag_density[0] = x[p]  # non-polymorphous case
                else:
                    polymorph = params[3]  # polymorphous case
                    poly = list(sample.structure[layer][element].polymorph).index(polymorph)
                    #poly = np.where(sample.structure[layer][element].polymorph == polymorph)
                    sample.structure[layer][element].mag_density[poly] = x[p]

    # the purpose of this code is to run the script and change the values accordingly
    if use_script:
        sample = readScript(sample,script)
    # include th script functionality here
    return sample, backS, scaleF

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def readScript(sample, script):

    variables = dict()
    for line in script:
        if len(line) == 2:  # setting a variable as a value
            key = line[0].strip(' ')
            function = line[1].split('(')[0].strip(' ')
            variables[key] = 0

            params = line[1].strip(' ').strip(function)  # parameters

            params = params.strip('(')
            params = params.strip(')')
            params = params.strip(' ')
            params = params.split(',')
            if function.lower() == 'getroughness':
                layer = int(params[0].strip(' '))
                K = params[1].strip(' ')
                variables[key] = float(sample.getRoughness(layer, K))

            elif function.lower() == 'getdensity':
                layer = int(params[0].strip(' '))
                K = params[1].strip(' ')
                variables[key] = float(sample.getDensity(layer, K))

            elif function.lower() == 'getthickness':
                layer = int(params[0].strip(' '))
                K = params[1].strip(' ')
                variables[key] = float(sample.getThickness(layer, K))

            elif function.lower() == 'gettotalthickness':
                start_layer = int(params[0].strip(' '))
                end_layer = int(params[1].strip(' '))
                K = params[2].strip(' ')
                variables[key] = float(sample.getTotalThickness(start_layer,end_layer, K))

        elif len(line) == 1:  # setting functions
            function = line[0].split('(')[0].strip(' ')

            params = line[0].strip(' ').strip(function)  # parameters

            params = params.strip('(')
            params = params.strip(')')
            params = params.strip(' ')
            params = params.split(',')

            if function.lower() == 'setroughness':
                layer = int(params[0])
                identifier = params[1].strip(' ')
                key = params[2].strip(' ')
                if isfloat(key):
                    sample.setRoughness(layer, identifier, float(key))
                else:
                    sample.setRoughness(layer, identifier, variables[key])

            elif function.lower() == 'setdensity':
                layer = int(params[0])
                identifier = params[1].strip(' ')
                key = params[2].strip(' ')
                if isfloat(key):
                    sample.setDensity(layer, identifier, float(key))
                else:
                    sample.setDensity(layer, identifier, variables[key])

            elif function.lower() =='setthickness':
                layer = int(params[0])
                identifier = params[1].strip(' ')
                key = params[2].strip(' ')
                if isfloat(key):
                    sample.setThickness(layer, identifier, float(key))
                else:
                    sample.setThickness(layer, identifier, variables[key])

            elif function.lower() == 'setcombinedthickness':
                layer_start = int(params[0])
                layer_end = int(params[1])
                identifier = params[2].strip(' ')
                key = params[3].strip(' ')

                if isfloat(key):
                    sample.setCombinedThickness(layer_start, layer_end, identifier, float(key))
                else:
                    sample.setCombinedThickness(layer_start, layer_end, identifier, variables[key])

            elif function.lower() == 'setratio':
                layer = int(params[0])
                symbol = params[1].strip(' ')
                identifier1 = params[2].strip(' ')
                identifier2 = params[3].strip(' ')
                ratio = params[3].strip(' ')
                sample.setRatio(layer, symbol, identifier1, identifier2, ratio)

    return sample
    # script must be a list of lines


def smooth_function(R):
    """
    Purpose: Uses the Savgol filter to smooth out the data
    :param R: Reflectivity data
    :return: Smoothed reflectivity data
    """
    return signal.savgol_filter(R, window_length=11, polyorder=5, mode="nearest")

def rolling_average(R, window, iter=1):
    # Uses rolling average to smooth out the data
    for i in range(iter):
        if i == 0:
            kernel = np.ones(9) / 9
            R = np.convolve(R, kernel, mode='same')
        else:
            kernel = np.ones(5) / 5
            R = np.convolve(R, kernel, mode='same')
    return R

def total_variation(R, Rsim):
    """
    Purpose: Calculate the difference in the total variation between R and Rsim
    :param R: Smoothed reflectivity data
    :param Rsim: Reflectivity simulation
    :return: Difference in total variation
    """

    totVar = sum([abs(R[idx+1]-R[idx]) for idx in range(len(R)-1)])/len(R)  # total variation in fitted data (arc length)
    totVarSim = sum([abs(Rsim[idx+1]-Rsim[idx]) for idx in range(len(Rsim)-1)])/len(Rsim)  # total variation in simulation

    variation = abs(totVar-totVarSim)
    return variation


def scanCompute(x, *args):
    """
    Purpose: Calculate the cost function for the data fitting algorithms
    :param x: List of parameters values
    :param args: List of required parameters for cost function calculation
    :return:
    """

    fun = 0  # used to minimize the global optimization
    gamma = 0

    sample = args[0]  # slab class
    scans = args[1]  # list of scan names to use in cost function
    data = args[2]  # Data for all scans
    backS = args[3]  # background shift
    scaleF = args[4]  # scaling factor
    parameters = args[5]  # list of parameters to change
    sBounds = args[6]  # defines the bounds of the scans
    sWeights = args[7]  # defines the weights for each scan
    objective = args[8]  # defines methodology in cost function (chi-square, L1-norm, L2-norm)
    shape_weight = args[9]  # total variation weighting
    optimizeSave = args[10]  # save optimization values
    r_scale = args[11]  # defines how to transform R [log(x), ln(x), x, qz^4]
    smooth_dict = args[12]  # dictionary of already smoothed out data
    script = args[13]
    use_script = args[14]

    # determines if saving will be done in callback function
    if optimizeSave:
        x_vars.append(x)

    # change the sample parameters
    sample, backS, scaleF = changeSampleParams(x, parameters, sample, backS, scaleF,script, use_script=use_script)


    i = 0 # keeps track of which boundary to use
    for scan in scans:

        scanType = scan[1]  # retrieves scan type
        name = scan[2]  # name of scan

        Rsmooth = smooth_dict[name]['Data'][2]  # retrieve smoothed data scan

        xbound = sBounds[i]  # scan boundaries
        weights = sWeights[i]  # scan weights
        i = i + 1
        background_shift = float(backS[name])  # retrieves the background shift
        scaling_factor = float(scaleF[name])  # retrieves the scaling factor
        if scanType == 'Reflectivity':  # reflectivity scan case
            myDataScan = data[name]
            myData = myDataScan['Data']  # retrieves original data
            E = myDataScan['Energy']  # retrieves energy pf scan
            pol = myDataScan['Polarization']  # retrieves polarization of scan
            Rdat = np.array(myData[2])  # retrieves the reflectivity

            fun_val = 0  # for cost function calculation


            qz = np.array(myData[0])  # retrieves momentum transfer
            # calculate simulation
            qz, Rsim = sample.reflectivity(E, qz, bShift=background_shift, sFactor=scaling_factor)
            Rsim = Rsim[pol]

            j = [x for x in range(len(qz)) if qz[x] > xbound[0][0] and qz[x] < xbound[-1][-1]]
            qz = qz[j]
            if len(Rsim) != len(j):
                Rsim = Rsim[j]
            if len(Rdat) != len(j):
                Rdat = Rdat[j]
            if len(Rsmooth) != len(j):
                Rsmooth = Rsmooth[j]
            # transforms R depending on users selection
            if r_scale == 'log(x)':
                Rsim = np.log10(Rsim)
                Rdat = np.log10(Rdat)
                Rsmooth = np.log10(Rsmooth)
            elif r_scale == 'ln(x)':
                Rsim = np.log(Rsim)
                Rdat = np.log(Rdat)
                Rsmooth = np.log(Rsmooth)
            elif r_scale == 'qz^4':
                Rsim = np.multiply(Rsim, np.power(qz,4))
                Rdat = np.multiply(Rdat, np.power(qz,4))
                Rsmooth = np.multiply(Rsmooth, np.power(qz,4))
            elif r_scale == 'x':
                pass

            # calculate the cost function
            m = 0  # keeps track of total number of data values used for each scan
            for b in range(len(xbound)):
                lw = xbound[b][0]  # lower boundary
                up = xbound[b][1]  # upper boundary
                w = weights[b]  # weights

                # determines indices of qz values found in boundary
                idx = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] < up]  # determines index boundaries
                n = len(idx)
                m = m + n  # updates number of data points used

                # calculates cost function depending on objective function used
                if n != 0:
                    if objective == 'Chi-Square':
                        fun_val = fun_val + sum((Rdat[idx]-Rsim[idx])**2/abs(Rsim[idx]))*w
                    elif objective == 'L1-Norm':
                        fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                    elif objective == 'L2-Norm':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx])**2) * w

            fun = fun + fun_val/m  # updates cost function

            # calculates total variation over entire boundary
            var_idx = [x for x in range(len(qz)) if qz[x] >= xbound[0][0] and qz[x] < xbound[-1][1]]
            gamma = gamma + total_variation(Rsmooth[var_idx], Rsim[var_idx])/len(Rsmooth[var_idx])

        elif scanType == 'Energy':
            myDataScan = data[name]
            myData = myDataScan['Data']  # retrieves energy scan data
            Theta = myDataScan['Angle']  # retrieves angle in degrees
            Rdat = np.array(myData[2])  # obtains reflectivity data
            E = np.array(myData[3])  # retrieves list of energy
            pol = myDataScan['Polarization']  # retrieves polarization of scan
            fun_val = 0

            n = len(Rdat)

            # calculates simulation
            E, Rsim = sample.energy_scan(Theta, E)
            Rsim = Rsim[pol]

            j = [x for x in range(len(E)) if E[x] > xbound[0][0] and E[x] < xbound[-1][-1]]
            E = E[j]

            if len(Rsim) != len(j):
                Rsim = Rsim[j]
            if len(Rdat) != len(j):
                Rdat = Rdat[j]
            if len(Rsmooth) != len(j):
                Rsmooth = Rsmooth[j]
            # transforms data and simulation as specified by user
            if r_scale == 'log(x)':
                Rsim = np.log10(Rsim)
                Rdat = np.log10(Rdat)
                Rsmooth = np.log10(Rsmooth)
            elif r_scale == 'ln(x)':
                Rsim = np.log(Rsim)
                Rdat = np.log(Rdat)
                Rsmooth = np.log(Rsmooth)
            elif r_scale == 'qz^4':
                qz = np.sin(Theta*np.pi/180)*(E * 0.001013546143)
                Rsim = np.multiply(Rsim, np.power(qz, 4))
                Rdat = np.multiply(Rdat, np.power(qz, 4))
                Rsmooth = np.multiply(Rsmooth, np.power(qz,4))
            elif r_scale == 'x':
                pass

            m = 0  # keeps track of total number of points used in specified boundary
            for b in range(len(xbound)):
                lw = xbound[b][0]  # upper boundary
                up = xbound[b][1]  # lower boundary
                w = weights[b]  # weights

                idx = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]  # determines index boundaries

                n = len(idx)
                m = m + n  # update number of points used

                # determines which objective function to use
                if len(idx) != 0:
                    if objective == 'Chi-Square':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2 / abs(Rsim[idx])) * w
                    elif objective == 'L1-Norm':
                        fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                    elif objective == 'L2-Norm':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx])**2) * w

            fun = fun + fun_val/m  # calculates cost function

            # calculates the total variation for entire boundary
            var_idx = [x for x in range(len(E)) if E[x] >= xbound[0][0] and E[x] < xbound[-1][1]]
            gamma = gamma + total_variation(Rsmooth[var_idx], Rsim[var_idx])/len(Rsmooth[var_idx])

    fun = fun + gamma*shape_weight  # adds the total variation to the cost function

    return fun

def differential_evolution(sample, data_info, data,scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale, smooth_dict, script, use_script=False):
    # performs the differential evolution global optimization
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale, smooth_dict,script,use_script]  # required format for function scanCompute

    p=True
    if goParam[8] == 'True':
        p = True
    else:
        p = False
    # This line will be used to select and use different global optimization algorithms
    ret = optimize.differential_evolution(scanCompute, bounds, args=params, strategy=goParam[0], maxiter=int(goParam[1]),
                                          popsize=int(goParam[2]),tol=float(goParam[3]), atol=float(goParam[4]),
                                          mutation=(float(goParam[5]), float(goParam[6])), recombination=float(goParam[7]),
                                          polish=p, init=goParam[9], updating=goParam[10], disp=True,
                                          callback=cb.stop_evolution)
    x = ret.x
    fun = ret.fun


    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    return x, fun

def shgo(sample, data_info, data, scan, backS, scaleF, parameters, bounds, sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale, smooth_dict,script, use_script=False):
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale, smooth_dict, script, use_script]  # required format for function scanCompute

    p = None
    if goParam[0] == 'None' or goParam[0] == None:
        p = None
    else:
        p = int(goParam[0])
    ret = optimize.shgo(scanCompute, bounds, args=tuple(params), n=p, iters=int(goParam[1]),sampling_method=goParam[2],
                        options={'disp': True}, callback=cb.stop_simplicial)
    x = ret.x
    fun = ret.fun

    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun

def dual_annealing(sample, data_info, data, scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale, smooth_dict,script, use_script=False):
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale, smooth_dict,script, use_script]

    p = True
    if goParam[6] == 'True':
        p = False
    else:
        p = True

    ret = optimize.dual_annealing(scanCompute, bounds, args=params, maxiter=int(goParam[0]), initial_temp=float(goParam[1]),
                                  restart_temp_ratio=float(goParam[2]), visit=float(goParam[3]), accept=float(goParam[4]),
                                  maxfun=float(goParam[5]), no_local_search=p, callback=cb.stop_annealing)
    x = ret.x
    fun = ret.fun


    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun

def least_squares(x0, sample, data_info, data, scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale, smooth_dict, script, use_script=False):
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data, backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, True, r_scale, smooth_dict, script, use_script]

    diff = goParam[8]
    _max = goParam[9]

    if diff.upper() == 'NONE':
        diff = None
    else:
        diff = float(diff)
    if _max.upper() == 'NONE':
        _max = None
    else:
        _max = float(_max)

    if goParam[1] == 'lm':
        result = optimize.least_squares(scanCompute, x0, args=params, jac=goParam[0], method=goParam[1],
                                        ftol=float(goParam[2]), xtol=float(goParam[3]), gtol=float(goParam[4]),
                                        x_scale=float(goParam[5]), loss=goParam[6], f_scale=float(goParam[7]),
                                        diff_step=diff,
                                        max_nfev=_max)
    else:
        result = optimize.least_squares(scanCompute, x0,bounds=bounds, args=params, jac=goParam[0], method=goParam[1],
                                        ftol=float(goParam[2]), xtol=float(goParam[3]), gtol=float(goParam[4]),
                                        x_scale=float(goParam[5]), loss=goParam[6], f_scale=float(goParam[7]),
                                        diff_step=diff,
                                        max_nfev=_max)



    x = result.x
    fun = result.cost

    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun

def direct(sample, data_info, data,scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale, smooth_dict,script, use_script=False):
    # performs the differential evolution global optimization
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale, smooth_dict, script, use_script]  # required format for function scanCompute

    # checking if locally biased
    p=True
    if goParam[3] == 'True':
        p = True
    else:
        p = False
    # This line will be used to select and use different global optimization algorithms
    #ret = optimize.direct(scanCompute, bounds, )
    #x = ret.x
    #fun = ret.fun

    x =[]
    fun = []

    #print('Chi: ' + str(fun))
    #print('Fitting parameters: ', x)

    return x, fun

def return_x():
    global x_vars
    return x_vars

class MinimizeStopper(object):
    def __init__(self, max_sec=0.3):
        self.max_sec = max_sec
        self.start = time.time()

    def __call__(self, xk=None, convergence=None):
        elapsed = time.time() - self.start
        if elapsed > self.max_sec:
            print("Terminating optimization: time limit reached")
            return True
        else:
            # you might want to report other stuff here
            # print("Elapsed: %.3f sec" % elapsed)
            return False