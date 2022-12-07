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

global x_vars
x_vars = []

def changeSampleParams(x, parameters, sample, backS, scaleF):
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
            name = params[1]
            if name == 'ALL SCANS':
                for key in list(scaleF.keys()):
                    scaleF[key] = str(x[p])
            else:
                scaleF[name] = str(x[p])

            sample.scaling_factor = x[p]
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
                pass
            elif mode == 'MAGNETIC':
                pass

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

                poly = np.where(sample.structure[layer][element].polymorph == polymorph)  # determines location of polymorph

                # sets poly_ratio value making sure sum equals to 1
                if poly == 0:
                    sample.structure[layer][element].poly_ratio[0] = x[p]
                    sample.structure[layer][element].poly_ratio[1] = ratio
                elif poly == 1:
                    sample.structure[layer][element].poly_ratio[1] = x[p]
                    sample.structure[layer][element].poly_ratio[0] = ratio

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

    return sample, backS, scaleF



def smooth_function(R):

    return signal.savgol_filter(R, window_length=11, polyorder=5, mode="nearest")

def rolling_average(R, window, iter=1):

    for i in range(iter):
        if i == 0:
            kernel = np.ones(9) / 9
            R = np.convolve(R, kernel, mode='same')
        else:
            kernel = np.ones(5) / 5
            R = np.convolve(R, kernel, mode='same')

    return R

def total_variation(R, Rsim):
    # calculates the total variation (arc-length) of the function
    # this will be a good judge of whether or not the shape of the function is right
    #Rnew = R-Rsim
    totVar = sum([abs(R[idx+1]-R[idx]) for idx in range(len(R)-1)])/len(R)  # total variation in fitted data
    totVarSim = sum([abs(Rsim[idx+1]-Rsim[idx]) for idx in range(len(Rsim)-1)])/len(Rsim)  # total variation in simulation

    variation = abs(totVar-totVarSim)
    return variation


def scanCompute(x, *args):

    fun = 0  # used to minimize the global optimization
    gamma = 0

    sample = args[0]  # slab class
    scans = args[1]  # data info
    data = args[2]  # data dict
    backS = args[3]
    scaleF = args[4]
    parameters = args[5]  # defines which parameters to change
    sBounds = args[6]  # defines the bounds of the scans
    sWeights = args[7]
    objective = args[8]
    shape_weight = args[9]
    optimizeSave = args[10]
    r_scale = args[11]

    # determines if saving will be done in callback function

    if optimizeSave:
        x_vars.append(x)


    sample, backS, scaleF = changeSampleParams(x, parameters, sample, backS, scaleF)


    i = 0 # keeps track of which boundary to use
    for scan in scans:

        scanType = scan[1]
        name = scan[2]

        xbound = sBounds[i]
        weights = sWeights[i]
        i = i + 1
        background_shift = float(backS[name])
        scaling_factor = float(scaleF[name])
        if scanType == 'Reflectivity':
            myDataScan = data[name]
            myData = myDataScan['Data']
            E = myDataScan['Energy']
            pol = myDataScan['Polarization']
            Rdat = np.array(myData[2])
            fun_val = 0


            qz = np.array(myData[0])
            qz, Rsim = sample.reflectivity(E, qz, bShift=background_shift, sFactor=scaling_factor)
            Rsim = Rsim[pol]

            if r_scale == 'log(x)':
                Rsim = np.log10(Rsim)
                Rdat = np.log10(Rdat)
            elif r_scale == 'ln(x)':
                Rsim = np.log(Rsim)
                Rdat = np.log(Rdat)
            elif r_scale == 'qz^4':
                Rsim = np.multiply(Rsim, np.power(qz,4))
                Rdat = np.multiply(Rdat, np.power(qz,4))
            elif r_scale == 'x':
                pass



            window = 5
            R = rolling_average(Rdat, window)

            m = 0
            for b in range(len(xbound)):
                lw = xbound[b][0]
                up = xbound[b][1]
                w = weights[b]

                idx = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] < up]  # determines index boundaries
                n = len(idx)
                m = m + n
                if n != 0:
                    if objective == 'Chi-Square':
                        fun_val = fun_val + sum((Rdat[idx]-Rsim[idx])**2/abs(Rsim[idx]))*w
                    elif objective == 'L1-Norm':
                        fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                    elif objective == 'L2-Norm':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx])**2) * w

            fun = fun + fun_val/m
            gamma = gamma + total_variation(R[2 * window:], Rsim[2 * window:])/len(R[2*window:])

        elif scanType == 'Energy':
            myDataScan = data[name]
            myData = myDataScan['Data']
            Theta = myDataScan['Angle']
            Rdat = np.array(myData[2])
            E = np.array(myData[3])
            pol = myDataScan['Polarization']
            fun_val = 0

            n = len(Rdat)


            E, Rsim = sample.energy_scan(Theta, E)
            Rsim = Rsim[pol]

            if r_scale == 'log(x)':
                Rsim = np.log10(Rsim)
                Rdat = np.log10(Rdat)
            elif r_scale == 'ln(x)':
                Rsim = np.log(Rsim)
                Rdat = np.log(Rdat)
            elif r_scale == 'qz^4':
                qz = np.sin(Theta*np.pi/180)*(E * 0.001013546143)
                Rsim = np.multiply(Rsim, np.power(qz, 4))
                Rdat = np.multiply(Rdat, np.power(qz, 4))
            elif r_scale == 'x':
                pass

            window = 5
            R = rolling_average(Rdat, window)

            m = 0
            for b in range(len(xbound)):
                lw = xbound[b][0]
                up = xbound[b][1]
                w = weights[b]

                idx = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]  # determines index boundaries

                n = len(idx)
                m = m + n

                if len(idx) != 0:
                    if objective == 'Chi-Square':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2 / abs(Rsim[idx])) * w
                    elif objective == 'L1-Norm':
                        fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                    elif objective == 'L2-Norm':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx])**2) * w

            fun = fun + fun_val/m
            gamma = gamma + total_variation(R[2 * window:], Rsim[2 * window:])/len(R[2*window:])

    fun = fun + gamma*shape_weight

    return fun

def differential_evolution(sample, data_info, data,scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale):
    # performs the differential evolution global optimization
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale]  # required format for function scanCompute

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

def shgo(sample, data_info, data, scan, backS, scaleF, parameters, bounds, sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale):
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale]  # required format for function scanCompute

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

def dual_annealing(sample, data_info, data, scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale):
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale]

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

def least_squares(x0, sample, data_info, data, scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale):
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data, backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, True, r_scale]

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

def direct(sample, data_info, data,scan,backS, scaleF, parameters, bounds,sBounds, sWeights, goParam, cb, objective, shape_weight, r_scale):
    # performs the differential evolution global optimization
    global x_vars
    x_vars = []

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data,backS, scaleF, parameters, sBounds, sWeights, objective, shape_weight, False, r_scale]  # required format for function scanCompute

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