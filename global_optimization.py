from scipy import optimize
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
                # may need to change this implementation
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

    sample = args[0]  # slab class
    scans = args[1]  # data info
    data = args[2]  # data dict
    sims = args[3]  # simulation dict
    parameters = args[4]  # defines which parameters to change
    sBounds = args[5]  # defines the bounds of the scans
    sWeights = args[6]

    sample = changeSampleParams(x, parameters, sample)

    i = 0 # keeps track of which boundary to use
    for scan in scans:

        scanType = scan[1]
        name = scan[2]

        xbound = sBounds[i]
        weights = sWeights[i]
        i = i + 1

        if scanType == 'Reflectivity':
            myDataScan = data[name]
            myData = myDataScan['Data']
            E = myDataScan['Energy']
            pol = myDataScan['Polarization']
            qz = np.array(myData[0])
            Rdat = np.log10(np.array(myData[2]))
            qz, Rsim = sample.reflectivity(E, qz)

            # need to toggle between log10 and not depending on the polarization
            Rsim = np.log10(Rsim[pol])*sample.scaling_factor + np.ones(len(Rsim[pol]))*sample.background_shift


            for b in range(len(xbound)):
                lw = xbound[b][0]
                up = xbound[b][1]
                w = weights[b]

                idx = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] <= up]  # determines index boundaries

                if len(idx) != 0:
                    chi2 = chi2 + sum((Rdat[idx]-Rsim[idx])**2/abs(Rsim[idx]))*w


            #chi2 = chi2 + sum((Rdat-Rsim)**2/abs(Rsim))

        elif scanType == 'Energy':
            myDataScan = data[name]
            myData = myDataScan['Data']
            Theta = myDataScan['Angle']
            E = np.array(myData[3])
            pol = myDataScan['Polarization']

            Rdat = np.log(np.array(myData[2]))
            E, Rsim = sample.energy_scan(Theta, E)
            Rsim = np.log(Rsim[pol])

            chi2 = chi2 + sum((Rdat-Rsim)**2/abs(Rsim))

    return chi2

def differential_evolution(fname,scan, parameters, bounds,sBounds, sWeights, goParam):

    sample = ReadSampleHDF5(fname)  # import the sample information

    data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data, sims, parameters, sBounds, sWeights]  # required format for function scanCompute

    p=True
    if goParam[8] == 'True':
        p = True
    else:
        p = False
    # This line will be used to select and use different global optimization algorithms
    ret = optimize.differential_evolution(scanCompute, bounds, args=params, strategy=goParam[0], maxiter=int(goParam[1]),
                                          popsize=int(goParam[2]),tol=float(goParam[3]), atol=float(goParam[4]),
                                          mutation=(float(goParam[5]), float(goParam[6])), recombination=float(goParam[7]),
                                          polish=p, init=goParam[9], updating=goParam[10])
    x = ret.x
    fun = ret.fun


    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    return x, fun

def shgo(fname, scan,parameters, bounds, sBounds, sWeights, goParam):

    sample = ReadSampleHDF5(fname)  # import the sample information
    data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data, sims, parameters, sBounds, sWeights]  # required format for function scanCompute

    p = None
    if goParam[0] == 'None':
        p = None
    else:
        p = int(goParam[0])
    ret = optimize.shgo(scanCompute, bounds, args=tuple(params), n=p, iters=int(goParam[1]),sampling_method=goParam[2])
    x = ret.x
    fun = ret.fun

    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun

def dual_annealing(fname, scan, parameters, bounds,sBounds, sWeights, goParam):
    sample = ReadSampleHDF5(fname)  # import the sample information

    data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    scans = []
    for s, info in enumerate(data_info):
        if info[2] in scan:
            scans.append(info)

    params = [sample, scans, data, sims, parameters, sBounds, sWeights]

    p = True
    if goParam[6] == 'True':
        p = False
    else:
        p = True

    ret = optimize.dual_annealing(scanCompute, bounds, args=params, maxiter=int(goParam[0]), initial_temp=float(goParam[1]),
                                  restart_temp_ratio=float(goParam[2]), visit=float(goParam[3]), accept=float(goParam[4]),
                                  maxfun=float(goParam[5]), no_local_search=p)
    x = ret.x
    fun = ret.fun


    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun