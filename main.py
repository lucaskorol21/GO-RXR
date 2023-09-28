"""

Library: main.py
Version: 0.1
Author: Lucas Korol
Institution: University of Saskatchewan
Last Updated: April 6, 2023
Python: version 3.7

Purpose: This script contains functions that are used to update the form factor database

Imported Libraries

numpy (version 1.21.4) - used for array manipulation

Note - This python script is also used to test out function and analyze data
"""
import os

import numpy as np
import scipy.interpolate

import data_structure as ds
import matplotlib.pyplot as plt
from numba import *

def getFF(name):
    """
    Purpose: Retrieve the form factor information from the database
    :param name: Name of the form factor
    :return: data - the data
             source - the source of the form factor
             updated - the date of when the form factor was last updated
    """

    import pickle

    # opens the current form factor database
    with open('form_factor.pkl', 'rb') as f:
        my_dict = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    my_dict.close()
    data = ''
    source = ''
    updated = ''
    if name not in my_dict.keys():
        print(name, ' not in keys')
    else:
        data = my_dict[name]['Data']
        source = my_dict[name]['Source']
        updated = my_dict[name]['Updated']

    return data, source, updated

def getFFM(name):
    """
    Purpose: Retrieve the form factor information from the database
    :param name: Name of the form factor
    :return: data - the data
             source - the source of the form factor
             updated - the date of when the form factor was last updated
    """
    import pickle

    # opens the current form factor database
    with open('form_factor_magnetic.pkl', 'rb') as f:
        my_dict = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    my_dict.close()

    data = ''
    source = ''
    updated = ''
    if name not in my_dict.keys():
        print(name, ' not in keys')
    else:
        data = my_dict[name]['Data']
        source = my_dict[name]['Source']
        updated = my_dict[name]['Updated']

    return data, source, updated

def createFF(directory, name , source, date):
    """
    Purpose: Add a new or change form factor in form factor database
    :param directory: Path to the desired form factor textfile
    :param name: Desired name of form factor (e.g. A, Mn, Mn2+...)
    :param source: Name of place the form factor was retrieved (e.g. The Center for X-ray Optics)
    :param date: Date the change was made (entered as a sting)
    """
    import os
    import pickle

    # opens the current form factor database
    with open('form_factor.pkl', 'rb') as f:
        my_dict = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    my_dict.close()

    data = np.loadtxt(directory)

    # includes the new form factor into the database
    my_dict[name]['Data'] = data  # stores the form factor data
    my_dict[name]['Source'] = source  # stores the url of the source
    my_dict[name]['Updated'] = date  # stores the date last updated

    with open('form_factor.pkl', 'wb') as handle:
        pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def createFFM(directory, name ,source, date):
    """
        Purpose: Add a new or change form factor in form factor database
        :param directory: Path to the desired form factor textfile
        :param name: Desired name of form factor (e.g. A, Mn, Mn2+...)
        :param source: Name of place the form factor was retrieved (e.g. The Center for X-ray Optics)
        :param date: Date the change was made (entered as a sting)
        """
    import os
    import pickle

    # opens the current form factor database
    with open('form_factor_magnetic.pkl', 'rb') as f:
        my_dict = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    my_dict.close()

    data = np.loadtxt(directory)

    # includes the new form factor into the database
    my_dict[name]['Data'] = data  # stores the form factor data
    my_dict[name]['Source'] = source  # stores the url of the source
    my_dict[name]['Updated'] = date  # stores the date last updated

    with open('form_factor_magnetic.pkl', 'wb') as handle:
        pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    import time as time
    from scipy.interpolate import UnivariateSpline
    import material_structure as ms
    import material_model as mm
    from scipy.special import erf

    import data_structure as ds
    """
    #sample = ms.slab(3)
    #sample.addlayer(0,'SrTiO3', 10, roughness=[0.5,2.5,1.3])
    #sample.addlayer(1, 'LaMnO3', 52.6, roughness=[1.25,0.72,0.2])
    #sample.addlayer(2, 'CCO', [15,10,5], density=[0,0.082,0.042], roughness=[0,1.8,2.5])
    sample = ms.slab(12)
    sample.addlayer(0, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(0,'A', ['Sr','La'], [1,0],sf=['Sr','La'])
    sample.polymorphous(0, 'B', ['Ti', 'Mn'], [1, 0], sf=['Ti','Mn'])
    sample.addlayer(1, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(1, 'A', ['Sr', 'La'], [0.95, 0.05],sf=['Sr','La'])
    sample.polymorphous(1, 'B', ['Ti', 'Mn'], [1, 0], sf=['Ti','Mn'])
    sample.addlayer(2, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(2, 'A', ['Sr', 'La'], [0.75, 0.25],sf=['Sr','La'])
    sample.polymorphous(2, 'B', ['Ti', 'Mn'], [0.95, 0.05], sf=['Ti','Mn'])
    sample.addlayer(3, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(3, 'A', ['Sr', 'La'], [0.25, 0.75],sf=['Sr','La'])
    sample.polymorphous(3, 'B', ['Ti', 'Mn'], [0.78, 0.22], sf=['Ti','Mn'])
    sample.addlayer(4, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(4, 'A', ['Sr', 'La'], [0.05, 0.95],sf=['Sr','La'])
    sample.polymorphous(4, 'B', ['Ti', 'Mn'], [0.22, 0.78], sf=['Ti','Mn'])
    sample.addlayer(5, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(5, 'A', ['Sr', 'La'], [0, 1],sf=['Sr','La'])
    sample.polymorphous(5, 'B', ['Ti', 'Mn'], [0.05, 0.95], sf=['Ti','Mn'])
    sample.addlayer(6, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(6, 'A', ['Sr', 'La'], [0, 1],sf=['Sr','La'])
    sample.polymorphous(6, 'B', ['Ti', 'Mn'], [0, 1], sf=['Ti','Mn'])
    sample.addlayer(7, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(7, 'A', ['Sr', 'La'], [0, 1],sf=['Sr','La'])
    sample.polymorphous(7, 'B', ['Ti', 'Mn'], [0, 1], sf=['Ti','Mn'])
    sample.addlayer(8, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(8, 'A', ['Sr', 'La'], [0, 1],sf=['Sr','La'])
    sample.polymorphous(8, 'B', ['Ti', 'Mn'], [0, 1], sf=['Ti','Mn'])
    sample.addlayer(9, 'ABO3', 3.905, density=[0.028,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(9, 'A', ['Sr', 'La'], [0, 1],sf=['Sr','La'])
    sample.polymorphous(9, 'B', ['Ti', 'Mn'], [0, 1],sf=['Ti','Mn'])
    sample.addlayer(10, 'ABO3', 3.905, density=[0.005,0.028,0.084], roughness=[0,0,0])
    sample.polymorphous(10, 'A', ['Sr', 'La'], [0, 1], sf=['Sr','La'])
    sample.polymorphous(10, 'B', ['Ti', 'Mn'], [0, 1], sf=['Ti','Mn'])
    sample.addlayer(11, 'CCO', [15,10,5], density=[0,0.082,0.042], roughness=[0,1.8,2.5], linked_roughness=[0,0.5,0.5])


    sample.energy_shift()
    sample.plot_density_profile()
    plt.show()

    sample.eShift['Mn'] = -1.1
    sample.eShift['La'] = 0.5

    constant_energy = [400,500,455,642.0,700,833]
    number = [1,3,5,7,9,11]
    theta = np.linspace(0.1,89.9,134)

    example_file = ds.DataFile()

    for idx,E in enumerate(constant_energy):
        qz = np.sin(theta * np.pi / 180) * (E * 0.001013546143)
        qz, R = sample.reflectivity(E,qz)

        example_file.addReflectivityScan(number[idx], E, 'S', qz, R['S'])
        example_file.addReflectivityScan(number[idx]+1, E, 'P', qz, R['P'])


    constant_theta = [5,10,15]
    E_Mn = np.linspace(630,670,211)
    E_Ti = np.linspace(450,480, 211)
    E_La = np.linspace(830,860,211)

    number = [13,15,17]
    for idx,angle in enumerate(constant_theta):
        E, R = sample.energy_scan(angle, E_Ti)
        example_file.addEnergyScan(number[idx], 450, angle, 'S', E_Ti, R['S'])
        example_file.addEnergyScan(number[idx]+1, 450, angle, 'P', E_Ti, R['P'])

    number = [19, 21, 23]
    for idx, angle in enumerate(constant_theta):
        E, R = sample.energy_scan(angle, E_Mn)
        example_file.addEnergyScan(number[idx], 630, angle, 'S', E_Mn, R['S'])
        example_file.addEnergyScan(number[idx] + 1, 630, angle, 'P', E_Mn, R['P'])

    number = [25, 27, 29]
    for idx, angle in enumerate(constant_theta):
        E, R = sample.energy_scan(angle, E_La)
        example_file.addEnergyScan(number[idx], 830, angle, 'S', E_La, R['S'])
        example_file.addEnergyScan(number[idx] + 1, 830, angle, 'P', E_La, R['P'])

    example_file.saveData('example2.h5')
    """

    """
    ##fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_v4.h5"

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    #struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')
    import ast
    x = '[0.22,"0.65"]'
    try:
        result = ast.literal_eval(x)
        # Do something with the result
    except (ValueError, SyntaxError) as e:
        # Handle the exception (e.g., invalid input)
        print("Error:", e)

    for val in result:
        print(isinstance(val, str))
    # Global Minimum Example
    """
    """
    sample = ds.ReadSampleHDF5(fname)
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    sample.energy_shift()
    name = '70_E814.75_Th25.0_S'
    energy = data_dict[name]['Data'][3]

    bShift = 0
    sFactor = 1

    Elen = len(energy)
    R = {'S': np.zeros(Elen),
         'P': np.zeros(Elen),
         'AL': np.zeros(Elen),
         'LC': np.zeros(Elen),
         'RC': np.zeros(Elen),
         'AC': np.zeros(Elen)}

    thickness, density, density_magnetic = sample.density_profile(step=0.1)  # Computes the density profile
    # Magnetic Scattering Factor
    sfm = dict()
    sf = dict()

    # Non-Magnetic Scattering Factor
    for e in sample.find_sf[0].keys():
        dE = float(sample.eShift[sample.find_sf[0][e]])
        scale = float(sample.ff_scale[sample.find_sf[0][e]])
        sf[e] = mm.find_form_factor(sample.find_sf[0][e], energy + dE, False) * scale
    # Magnetic Scattering Factor
    for em in sample.find_sf[1].keys():
        dE = float(sample.mag_eShift[sample.find_sf[1][em]])
        scale = float(sample.ffm_scale[sample.find_sf[1][em]])
        sfm[em] = mm.find_form_factor(sample.find_sf[1][em], energy + dE, True) * scale

    d_len = len(thickness)
    delta, beta = mm.IoR(density, sf, energy)  # gets absorptive and dispersive components of refractive index

    delta_m, beta_m = mm.MOC(density_magnetic, sfm, energy,
                          d_len)  # absorptive and dispersive components for magnetic components

    epsilon = 1 - 2 * delta + 1j * beta * 2  # dielectric constant

    # definition as described in Lott Dieter Thesis
    Q = beta_m + 1j * delta_m  # magneto-optical constant
    epsilon_mag = Q * epsilon * (-2)  # magneto-optical permittivity
    # retrieves the slabs at each energy using list comprehension

    number_slabs = [1,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,30,60,90,120,150,300,900]
    my_slabs = []

    my_time = np.zeros(len(number_slabs))

    num = 10
    Theta = 25.0
    for i in range(num):
        for idx, n in enumerate(number_slabs):
            R = {'S': np.zeros(Elen),
                 'P': np.zeros(Elen),
                 'AL': np.zeros(Elen),
                 'LC': np.zeros(Elen),
                 'RC': np.zeros(Elen),
                 'AC': np.zeros(Elen)}

            all_slabs = [np.arange(0, len(thickness), n, dtype=int) for E in range(len(energy))]


            # initializes the object for reflectivity computation using list comprehension


            Alist = [
                ms.generate_structure(thickness, sample.structure, all_slabs[s], epsilon[s], epsilon_mag[s], sample.layer_magnetized,
                                   sample.transition) for s in range(len(all_slabs))]

            h = 4.135667696e-15  # Plank's Constant [eV s]
            c = 2.99792450e10  # Speed of light in vacuum [cm/s]

            wavelength = h * c / (energy * 1e-10)

            # reflectivity computation using list comprehension

            start = time.time_ns()
            R = [ms.energy_reflectivity(Alist[int(E)], Theta, wavelength[int(E)], R, int(E), backS=bShift, scaleF=sFactor) for E in
                 range(len(all_slabs))]
            end = time.time_ns()

            my_time[idx] = my_time[idx] + (end-start)

            if i == 0:
                my_slabs.append(len(all_slabs[0]))
        print(i)

    my_time = my_time[1:]/num/1e6

    plt.figure()
    plt.plot(my_slabs[1:], my_time, 'o')
    plt.show()
    """
    """
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_v4.h5"
    struct_names, mag_names = mm._use_given_ff("//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3")
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    sample = ds.ReadSampleHDF5(fname)
    name = '64_836.04_S'

    E = 836.04
    qz = data_dict[name]['Data'][0]
    
   
    
    thickness, density, density_magnetic = sample.density_profile(step=0.1)
    sf = dict()  # scattering factors of non-magnetic components
    sfm = dict()  # scattering factors of magnetic components

    # print(self.find_sf[1])

    # Non-Magnetic Scattering Factor
    for e in sample.find_sf[0].keys():
        dE = float(sample.eShift[sample.find_sf[0][e]])  # retrieve the energy shift of each scattering factor
        scale = float(sample.ff_scale[sample.find_sf[0][e]])  # retrieve scaling factor of each scattering factor
        sf[e] = mm.find_form_factor(sample.find_sf[0][e], E + dE,
                                    False) * scale  # find the scattering factor at energy E + dE
    # Magnetic Scattering Factor
    for em in sample.find_sf[1].keys():
        dE = float(sample.mag_eShift[sample.find_sf[1][em]])
        scale = float(sample.ffm_scale[sample.find_sf[1][em]])
        sfm[em] = mm.find_form_factor(sample.find_sf[1][em], E + dE, True) * scale

    delta, beta = mm.index_of_refraction(density, sf, E)  # calculates depth-dependent refractive index components
    delta_m, beta_m = mm.magnetic_optical_constant(density_magnetic, sfm,
                                                   E)  # calculates depth-dependent magnetic components

    # definition of magneto-optical constant as described in Lott Dieter Thesis
    n = 1 + np.vectorize(complex)(-delta, beta)  # complex index of refraction
    epsilon = n ** 2  # dielectric constant computation

    # magneto-optical constant as defined in Lott Dieter Thesis
    Q = np.vectorize(complex)(beta_m, delta_m)
    epsilon_mag = Q * epsilon * 2 * (-1)

    my_slabs = ms.ALS(epsilon.real, epsilon.imag, epsilon_mag.real, epsilon_mag.imag, 0.01)  # performs the adaptive layer segmentation using Numba

    my_thickness1 = []
    my_epsilon1 = []
    for idx in my_slabs:
        idx = int(idx)
        my_thickness1.append(thickness[idx])
        my_thickness1.append(thickness[idx])
        my_thickness1.append(thickness[idx])

        my_epsilon1.append(0)
        my_epsilon1.append(delta[idx])
        my_epsilon1.append(0)

    my_slabs = ms.ALS(epsilon.real, epsilon.imag, epsilon_mag.real, epsilon_mag.imag,
                      1e-20)  # performs the adaptive layer segmentation using Numba

    my_thickness2 = []
    my_epsilon2 = []
    for idx in my_slabs:
        idx = int(idx)
        my_thickness2.append(thickness[idx])
        my_thickness2.append(thickness[idx])
        my_thickness2.append(thickness[idx])

        my_epsilon2.append(0)
        my_epsilon2.append(delta[idx])
        my_epsilon2.append(0)

    my_slabs = ms.ALS(epsilon.real, epsilon.imag, epsilon_mag.real, epsilon_mag.imag,
                      0.001)  # performs the adaptive layer segmentation using Numba

    my_thickness3 = []
    my_epsilon3 = []
    for idx in my_slabs:
        idx = int(idx)
        my_thickness3.append(thickness[idx])
        my_thickness3.append(thickness[idx])
        my_thickness3.append(thickness[idx])

        my_epsilon3.append(0)
        my_epsilon3.append(delta[idx])
        my_epsilon3.append(0)

    fig, axes = plt.subplots(1,2)

    from matplotlib import ticker
    # Plot the data on the top-left subplot
    axes[0].plot(my_thickness1, my_epsilon1, linewidth=0.5)
    axes[0].plot(thickness, delta)
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[0].tick_params(axis='y', labelleft=False)
    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].set_ylabel(r'$\mathrm{\delta \left(E \right)}$', fontsize=20)


    # Plot the data on the top-right subplot
    axes[1].plot(my_thickness3, my_epsilon3, linewidth=0.5)
    axes[1].plot(thickness, delta)
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].tick_params(axis='y', labelleft=False)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    #axes[1].legend(frameon=False)

    axes[0].set_xlim(30,45)
    axes[1].set_xlim(30, 45)
    axes[0].text(0.02, 0.98, "(a)", transform=axes[0].transAxes, fontsize=14, fontweight='bold', va='top')
    axes[1].text(0.02, 0.98, "(b)", transform=axes[1].transAxes, fontsize=14, fontweight='bold', va='top')
    # Adjust the layout
    plt.tight_layout()

    # Display the figure



    plt.figure(3)
    plt.plot(my_thickness2, my_epsilon2,linewidth=0.5)
    plt.plot(thickness, delta)
    plt.xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=16)
    plt.ylabel(r'$\mathrm{\delta \left(E \right)}$', fontsize=20)

    plt.xlim([30,45])
    # Set minor ticks
    plt.minorticks_on()

    # Set tick parameters
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.text(0.02, 0.98, '(c)', transform=plt.gca().transAxes, fontsize=14, fontweight='bold', va='top', ha='left')
    #plt.text(0.02, 0.98, "(c)", transform=axes[1].transAxes, fontsize=14, fontweight='bold', va='top')
    plt.show()
    """


    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    struct_names, mag_names = mm._use_given_ff("//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3")  # look for form factors in directory

    # set form factors found in file directory


    name = '64_836.04_S'

    E = 836.04
    qz = data_dict[name]['Data'][0]
    qz,R = sample.reflectivity(E,qz, precision = 1e-2, s_min=0.1)
    qz, R1 = sample.reflectivity(E, qz, precision = 1e-3, s_min=0.1)
    qz, R2 = sample.reflectivity(E, qz, precision=1e-5, s_min=0.1)
    qz, R3 = sample.reflectivity(E, qz, precision=1e-20, s_min=0.1)
    R = R['S']
    R1 = R1['S']
    R2 = R2['S']

    plt.figure()
    plt.plot(qz, R,'b')
    plt.plot(qz, R1,'r')
    plt.plot(qz, R2, 'y')
    plt.plot(qz, R2, 'k--')
    plt.legend(['precision=0.01','precision=0.001','precision=1e-5','precision=1e-20'], frameon=False)
    plt.ylabel('Fractional Reflected Intensity ($\mathrm{R/R_{0}}$)', fontsize=16)
    plt.xlabel(r'Momentum Transfer, $q_{z}$ ($\mathrm{\AA^{-1}}$)', fontsize=16)
    # Set minor ticks
    plt.minorticks_on()
    # Set tick parameters
    plt.tick_params(which='both', direction='in', top=True, right=True)
    #plt.tick_params(axis='y', labelleft=False)
    plt.yscale('log')
    plt.show()


    """
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"

    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    struct_names, mag_names = mm._use_given_ff("//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3")
    sample = ds.ReadSampleHDF5(fname)
    name = '64_836.04_S'

    E = 836.04
    qz = data_dict[name]['Data'][0]

    precision_array = [1e-2,1e-2,5e-3,1e-3,2e-4,5e-4,1e-4,2e-5,5e-5,1e-5,2e-6, 5e-6,1e-6, 5e-7, 1e-7, 5e-8, 1e-8]
    time_array = np.zeros(len(precision_array))
    slab_array = []
    print(len(qz))
    for i in range(100):
        print(i)
        for idx, prec in enumerate(precision_array):

            start = time.time_ns()
            sample.reflectivity(E, qz, precision=prec)
            end = time.time_ns()
            time_array[idx] = time_array[idx] + (end-start)

            thickness, density, density_magnetic = sample.density_profile()
            sf = dict()  # scattering factors of non-magnetic components
            sfm = dict()  # scattering factors of magnetic components

            # print(self.find_sf[1])

            # Non-Magnetic Scattering Factor
            for e in sample.find_sf[0].keys():
                dE = float(sample.eShift[sample.find_sf[0][e]])  # retrieve the energy shift of each scattering factor
                scale = float(sample.ff_scale[sample.find_sf[0][e]])  # retrieve scaling factor of each scattering factor
                sf[e] = mm.find_form_factor(sample.find_sf[0][e], E + dE,
                                         False) * scale  # find the scattering factor at energy E + dE
            # Magnetic Scattering Factor
            for em in sample.find_sf[1].keys():
                dE = float(sample.mag_eShift[sample.find_sf[1][em]])
                scale = float(sample.ffm_scale[sample.find_sf[1][em]])
                sfm[em] = mm.find_form_factor(sample.find_sf[1][em], E + dE, True) * scale

            delta, beta = mm.index_of_refraction(density, sf, E)  # calculates depth-dependent refractive index components
            delta_m, beta_m = mm.magnetic_optical_constant(density_magnetic, sfm,
                                                        E)  # calculates depth-dependent magnetic components

            # definition of magneto-optical constant as described in Lott Dieter Thesis
            n = 1 + np.vectorize(complex)(-delta, beta)  # complex index of refraction
            epsilon = n ** 2  # dielectric constant computation

            # magneto-optical constant as defined in Lott Dieter Thesis
            Q = np.vectorize(complex)(beta_m, delta_m)
            epsilon_mag = Q * epsilon * 2 * (-1)

            my_slabs = ms.ALS(epsilon.real, epsilon.imag, epsilon_mag.real, epsilon_mag.imag, prec)  # performs the adaptive layer segmentation using Numba
            if i == 0:
                slab_array.append(len(my_slabs))
            #thickness, density, density_mag = sample.density_profile()

    time_array = time_array/100/1e6
    plt.figure()
    plt.plot(slab_array[1:], time_array[1:] ,'o')
    plt.xlabel('Number of Slices', fontsize=12)
    plt.ylabel('Average Execution Time (ms)', fontsize=12)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.show()
    """
    """
    thickness, density, mag_density = sample.density_profile(step=0.01)
    idx = [i for i in range(len(thickness)) if thickness[i] <= 60]
    plt.figure()
    plt.plot(thickness[idx], density['Sr'][idx])
    plt.plot(thickness[idx],density['La'][idx],'--')
    plt.legend(['Sr','La'])
    plt.xlabel(r'Depth ($\AA$)')
    plt.ylabel(r'Density ($mol/cm^{3}$)')
    plt.show()
    """
    """
    # Jesus Data
    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)

    E = 460.76  # energy

    name = '35_460.76_S'
    #name = '31_635.99_S'

    qz_min = data_dict[name]['Data'][0][0]
    qz_max = data_dict[name]['Data'][0][-1]
    qz = np.linspace(qz_min, qz_max, 2000)

    qz, R = sample.reflectivity(E, qz)
    R = R['S']

    my_data = np.array([qz, R])

    z = np.linspace(0,20,100)

    val =sample.error_function(z, 0,10, True)
    plt.figure()
    plt.plot(z,val)
    plt.show()
    np.savetxt('E1_460.76_best.txt', my_data.transpose())
    """
    """

    from matplotlib import ticker
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_complete.h5"
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    name1 = '47_E429.58_Th10.0_S'

    E = data_dict[name1]['Data'][3]
    R = data_dict[name1]['Data'][2]
    plt.figure()
    plt.plot(E, R)
    plt.xlabel('Energy, E (eV)', fontsize=12)
    plt.ylabel('Fractional Reflected Intensity ($\mathrm{R/R_{0}}$)', fontsize=12)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #plt.tick_params(axis='y', labelleft=False)
    plt.minorticks_on()
    plt.show()

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()
    name2 = '49_E600.18_Th10.0_S'

    E2 = data_dict[name2]['Data'][3]
    R2 = data_dict[name2]['Data'][2]


    Eb, Rsb = sample.energy_scan(10.0, E2)
    Rsb = Rsb['S']

    sample.eShift['Mn2'] = 1.1
    sample.eShift['Mn3'] = 1.1
    Ed, Rs2 = sample.energy_scan(10.0, E2)
    Rs2 = Rs2['S']


    fig, axes = plt.subplots(1,2)
    axes[1].plot(E2, R2, label=r'Exp ($\sigma$)')
    axes[1].plot(Eb, Rsb, label=r'Calc ($\sigma$)')
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].tick_params(axis='y', labelleft=False)
    axes[1].legend([r'Exp', r'Calc'], loc='upper right', frameon=False)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].legend(frameon=False)

    axes[0].plot(E2, R2, label=r'Exp ($\sigma$)')
    axes[0].plot(Eb, Rs2, label=r'Calc ($\sigma$)')
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[0].tick_params(axis='y', labelleft=False)
    axes[0].legend([r'Exp', r'Calc'], loc='upper right', frameon=False)
    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].legend(frameon=False)
    font_props = {'family': 'sans-serif', 'size': 12, 'weight': 'bold'}
    axes[0].text(0.02, 0.98, "(a)", transform=axes[0].transAxes, fontsize=14, fontweight='bold', va='top')
    axes[1].text(0.02, 0.98, "(b)", transform=axes[1].transAxes, fontsize=14, fontweight='bold', va='top')


    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel(r'Momentum Transfer, $q_{z}$ ($\mathrm{\AA^{-1}}$)', fontsize=12)

    axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_x_axes().join(axes[0], shared_x)
    axes[1].get_shared_x_axes().join(axes[1], shared_x)
    axes[0].set_ylabel('Fractional Reflected Intensity ($\mathrm{R/R_{0}}$)', fontsize=12)
    axes[0].ticklabel_format(axis='y', scilimits=[-3, 3])

    plt.tight_layout()
    plt.show()
    """
    """
    E = 833

    Theta = np.arange(0.1, 90, 0.1)
    qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)

    # thickness variation
    sample_t1 = ms.slab(2)
    sample_t1.addlayer(0, 'SrTiO3', 50, roughness=0)
    sample_t1.addlayer(1, 'LaMnO3', 5, roughness=0)
    sample_t1.energy_shift()

    sample_t2 = ms.slab(2)
    sample_t2.addlayer(0, 'SrTiO3', 50, roughness=0)
    sample_t2.addlayer(1, 'LaMnO3', 20, roughness=0)
    sample_t2.energy_shift()

    sample_t3 = ms.slab(2)
    sample_t3.addlayer(0, 'SrTiO3', 50, roughness=0)
    sample_t3.addlayer(1, 'LaMnO3', 40, roughness=0)
    sample_t3.energy_shift()

    qz1, R1 = sample_t1.reflectivity(E, qz)
    qz2, R2 = sample_t2.reflectivity(E, qz)
    qz3, R3 = sample_t3.reflectivity(E, qz)

    R1['S'] = np.log10(R1['S'])
    R2['S'] = np.log10(R2['S'])
    R3['S'] = np.log10(R3['S'])

    plt.figure(1)
    plt.plot(qz1, (R1['S']-min(R1['S']))/(max(R1['S'])-min(R1['S']))+1)
    plt.plot(qz1, (R2['S']-min(R2['S']))/(max(R2['S'])-min(R2['S']))+0.5)
    plt.plot(qz1, (R3['S']-min(R3['S']))/(max(R3['S'])-min(R3['S']))+0)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #plt.tick_params(axis='y', labelleft=False)
    plt.minorticks_on()
    plt.legend([r'5 $\mathrm{\AA}$', r'20 $\mathrm{\AA}$', r'40 $\mathrm{\AA}$'], frameon=False, fontsize=10)
    plt.xlabel(r'Momentum Transfer, $q_{z}$ ($\mathrm{\AA^{-1}}$)', fontsize=12)
    plt.ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    plt.show()



    sample1 = ms.slab(2)
    sample1.addlayer(0,'SrTiO3',10, roughness=0)
    sample1.addlayer(1,'LaMnO3',80, roughness=0)
    sample1.energy_shift()

    sample2 = ms.slab(2)
    sample2.addlayer(0, 'SrTiO3', 10, roughness=3)
    sample2.addlayer(1, 'LaMnO3', 80, roughness = 0)
    sample2.energy_shift()

    sample3 = ms.slab(2)
    sample3.addlayer(0, 'SrTiO3', 10, roughness=6)
    sample3.addlayer(1, 'LaMnO3', 80, roughness=0)
    sample3.energy_shift()

    qz1, R1i = sample1.reflectivity(E, qz)
    qz2, R2i = sample2.reflectivity(E, qz)
    qz3, R3i = sample3.reflectivity(E, qz)

    sample1 = ms.slab(2)
    sample1.addlayer(0, 'SrTiO3', 10, roughness=0)
    sample1.addlayer(1, 'LaMnO3', 80, roughness=0)
    sample1.energy_shift()

    sample2 = ms.slab(2)
    sample2.addlayer(0, 'SrTiO3', 10, roughness=0)
    sample2.addlayer(1, 'LaMnO3', 80, roughness=3)
    sample2.energy_shift()

    sample3 = ms.slab(2)
    sample3.addlayer(0, 'SrTiO3', 10, roughness=0)
    sample3.addlayer(1, 'LaMnO3', 80, roughness=6)
    sample3.energy_shift()

    from matplotlib import ticker
    qz1, R1s = sample1.reflectivity(E, qz)
    qz2, R2s = sample2.reflectivity(E, qz)
    qz3, R3s = sample3.reflectivity(E, qz)


    R1i['S'] = np.log10(R1i['S'])
    R2i['S'] = np.log10(R2i['S'])
    R3i['S'] = np.log10(R3i['S'])

    R1s['S'] = np.log10(R1s['S'])
    R2s['S'] = np.log10(R2s['S'])
    R3s['S'] = np.log10(R3s['S'])

    fig, axes = plt.subplots(1,2)
    axes[0].plot(qz1, (R1i['S']-min(R1i['S']))/(max(R1i['S'])-min(R1i['S']))+1, label=r'$\sigma_{i}$=0 $\mathrm{\AA}$')
    axes[0].plot(qz1, (R2i['S']-min(R2i['S']))/(max(R2i['S'])-min(R2i['S']))+0.5, label=r'$\sigma_{i}$=3 $\mathrm{\AA}$')
    axes[0].plot(qz1, (R3i['S']-min(R3i['S']))/(max(R3i['S'])-min(R1i['S']))+0, label=r'$\sigma_{i}$=6 $\mathrm{\AA}$')
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[0].tick_params(axis='y', labelleft=False)
    axes[0].legend([r'Exp', r'Calc'], loc='upper right', frameon=False)
    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].legend(frameon=False)
    print(min(R1i['S']))
    axes[1].plot(qz1, (R1s['S']-min(R1s['S']))/(max(R1s['S'])-min(R1s['S']))+1, label=r'$\sigma_{s}$=0 $\mathrm{\AA}$')
    axes[1].plot(qz1, (R2s['S']-min(R2s['S']))/(max(R2s['S'])-min(R2s['S']))+0.5, label=r'$\sigma_{s}$=3 $\mathrm{\AA}$')
    axes[1].plot(qz1, (R3s['S']-min(R3s['S']))/(max(R3s['S'])-min(R3s['S']))+0, label=r'$\sigma_{s}$=6 $\mathrm{\AA}$')
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].tick_params(axis='y', labelleft=False)
    axes[1].legend([r'Exp', r'Calc'], loc='upper right', frameon=False)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].legend(frameon=False)



    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel(r'Momentum Transfer, $q_{z}$ ($\mathrm{\AA^{-1}}$)', fontsize=12)

    axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_x_axes().join(axes[0], shared_x)
    axes[1].get_shared_x_axes().join(axes[1], shared_x)
    axes[0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    font_props = {'family': 'sans-serif', 'size': 12, 'weight': 'bold'}
    axes[0].text(0.01, 0.98, "(a)", transform=axes[0].transAxes, va='top', **font_props)
    axes[1].text(0.01, 0.98, "(b)", transform=axes[1].transAxes, va='top', **font_props)

    plt.tight_layout()
    plt.show()

    hello = np.loadtxt("//cabinet/work$/lsk601/Downloads/SrTiO3_attenuation.txt")
    E = hello[:,0]
    attenuation = hello[:,1]

    plt.figure()
    plt.plot(E, attenuation)
    plt.xlabel('X-Ray Energy (eV)', fontsize=12)
    plt.ylabel(r'Attenuation length ($\times 10^{-6}$ meters)', fontsize=12)
    plt.grid(True)
    plt.show()
    """
    """
    sample = ds.ReadSampleHDF5(fname)

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')  # look for form factors in directory
    sample.energy_shift()
    Theta = 15
    energy = 640.2

    data = []
    Theta = np.arange(1,50.5,0.5)
    qz = np.sin(Theta * np.pi / 180) * (energy * 0.001013546143)
    my_qz = []
    energy = np.arange(455,470+1,1)
    my_energy = []
    qz_use = []
    for angle in Theta:
        E, R = sample.energy_scan(angle, energy)
        qz = np.sin(angle *np.pi/180) *((energy * 0.001013546143))
        my_qz.append(qz)
        my_energy.append(energy)
        R = R['S']*np.power(qz, 4)
        data.append(R)
        qz_use.append(np.sin(angle * np.pi / 180) * (energy * 0.001013546143))
        print(angle)


    data = np.array(data)
    my_qz = np.array(my_qz)
    my_energy = np.array(my_energy)





    from scipy import interpolate

    interp_func = scipy.interpolate.interp2d(Theta, energy, data.transpose(), kind='linear')

    Theta_new = np.arange(0.1,50.1,0.01)
    E_new = np.arange(455,470+0.01,0.01)

    data_new = interp_func(Theta_new, E_new)

    X, Y = np.meshgrid(Theta_new, E_new)
    qz_use = np.array(qz_use)

    import matplotlib.pyplot as plt

    from matplotlib.ticker import FormatStrFormatter
    from mpl_toolkits.mplot3d import Axes3D


    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(plt.cm.jet(np.linspace(0, 1, 256)**0.65))
    # plot surface
    ax.plot_surface(X, Y, data_new, cmap=cmap, alpha=1)
    ax.grid(False)
    # hide z-axis and turn off mesh
    ax.set_zticks([])

    # Set the formatter for the x-axis
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Set the formatter for the y-axis
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.0f'))

    # Set labels
    ax.set_xlabel(r'Grazing Angle,  $\theta_g$ $(degrees)$')
    ax.set_ylabel('Energy (eV)')


    # Show plot
    plt.show()
    
   
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)

    keys = ['26_452.77_S' ,'35_460.76_S','19_500.71_S', '31_635.99_S','22_640.99_S','24_644.02_S','33_834.59_S',
            '9_642.12_LC' ,'10_642.12_RC', '9-10_642.12_AC_Asymm', '13_644.03_LC','14_644.03_RC','13-14_644.03_AC_Asymm',
            '16_653.06_LC', '17_653.06_RC', '16-17_653.06_AC_Asymm']

    for key in keys:
        E = data_dict[key]['Energy']
        pol = data_dict[key]['Polarization']
        qz = data_dict[key]['Data'][0]

        # use to create new number of points!
        qz_min = qz[0]
        qz_max = qz[-1]
        number_points = 1000

        qz_new = np.linspace(qz_min,qz_max,num=number_points)
        Theta = np.arcsin(qz_new / E / (0.001013546247)) * 180 / np.pi  # initial angle

        qz_new, R = sample.reflectivity(E,qz_new)
        R = R[pol]


        sim_dict[key]['Data'] = np.array([qz_new, Theta, R])


        print('Done - ', key)

    ds.saveSimulationHDF5(fname, sim_dict)
    
    


    thickness, density, mag_density = sample.density_profile()

    my_keys = ['Sr', 'Ti', 'La', 'Mn','O']
    my_keys = ['Ti','Mn']


    d = 37.550
    idx = [i for i in range(len(thickness)) if thickness[i] < d]


    electron = {'Sr':2,'La':3,'Ti':4,'Mn2':2,'Mn3':3, 'O':-2}
    density['Mn'] = density['Mn2'] + density['Mn3']

    from scipy import integrate

    total = integrate.trapezoid(mag_density['Mn3'], x=thickness)
    print(total/10)
    #print(total/10/4)

    x = [4,7,10]
    y = [0.011817737262843186/4, 0.06735279450392744/7, 0.10131622378913185/10]
    plt.figure(1)
    plt.bar(x,y,)
    plt.xlabel('Sample (uc)')
    plt.ylabel('Magnetic Density per units cell')


    plt.figure(2)
    for key in my_keys:
        plt.plot(thickness, density[key], c=np.random.rand(3,))
    
    plt.ylabel('Density (mol/cm^3)')
    plt.xlabel('Thickness (angstroms)')
    plt.legend(my_keys)
    plt.show()
    
    oxidation = np.zeros(len(thickness))
    total = np.zeros(len(thickness))
    for key in electron.keys():
        oxidation = oxidation + electron[key]*density[key]
        total = total+density[key]

    oxidation = oxidation/total

    plt.figure()
    plt.plot(thickness[idx], oxidation[idx])
    plt.show()
    
    x = [4,7,10]
    y = [0.0007044,0.00220607,0.00279876]

    plt.figure(3)
    plt.plot(x,y)
    plt.ylabel('Average Magnetic Density (mol/cm^2)')
    plt.xlabel('Unit Cells (uc)')


    
    
    #KK = variationKK(E_prime,E0,E2,E4)
    

    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"

    struct_names, mag_names = mm._use_given_ff("//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas")  # look for form factors in directory

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()


    
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)

    keys = ['26_452.77_S', '35_460.76_S', '19_500.71_S', '31_635.99_S', '22_640.99_S', '24_644.02_S', '33_834.59_S',
            '9_642.12_LC', '10_642.12_RC', '9-10_642.12_AC_Asymm', '13_644.03_LC', '14_644.03_RC',
            '13-14_644.03_AC_Asymm',
            '16_653.06_LC', '17_653.06_RC', '16-17_653.06_AC_Asymm']


    
    destination = '\\cabinet\work$\lsk601\My Documents\Data_for_Jesus\RXR-Twente-E1-150K-simulation-0.1eV'
    for i in range(0,300+1,1):
        E = 450 + i*0.1
        Theta = np.linspace(0.1, 60, num=2000)
        qz = np.sin(Theta * np.pi / 180) * E * (0.001013546247)

        qz, R = sample.reflectivity(E, qz)

        filename = 'E1_' + str(E)
        dat = np.transpose(np.array([qz, R['S']]))
        file = r"\\cabinet\work$\lsk601\My Documents\Data_for_Jesus\RXR-Twente-E1-150K-simulation-0.1eV"
        file = file + '\E1_' + str(E)
        np.savetxt(file, dat)



    #np.savetxt('E1_' + str(E), dat)
    

    #E = np.arange(10,3000,1)
    #imaginary = np.zeros(len(E))
    #real = np.zeros(len(E))

    #data = np.transpose(np.array([E, imaginary, real]))

    #np.savetxt('Z.txt', data,fmt='%.4s')

    
    dir= "C:/Users/lsk601/PycharmProjects/MaterialReflection/Magnetic_Scattering_Factor"
    import pickle
    import os

    #with open('form_factor.pkl', 'rb') as f:
    #    my_dict = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    #f.close()

    my_dict = dict()

    source = 'The Center for X-ray Optics'
    date = 'April 6, 2023'
    for filename in os.listdir(dir):
        f = os.path.join(dir, filename)
        if not (filename.startswith('README')):
            name = filename.split('.')[0]
            data = np.loadtxt(f)
            my_dict[name] = dict()
            my_dict[name]['Data'] = data  # stores the form factor data
            my_dict[name]['Source'] = source  # stores the url of the source
            my_dict[name]['Updated'] = date  # stores the date last updated


    import pickle
    with open('form_factor_magnetic.pkl', 'wb') as handle:
        pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

    
    new = np.loadtxt('833.csv')
    t = new[:,0]

    text = np.loadtxt('833_ReMagX.txt')



    #plt.figure(1)
    #plt.suptitle('Figure')
    #plt.plot(new[:, 0], new[:, 3])
    #plt.plot(text[:,0], np.log10(text[:,1]))
    #plt.legend(['GO-RXR', 'ReMagX'])
    #plt.xlabel('Momentum Transfer (A^{-1})')
    #plt.ylabel('Reflectivity, R')



    #print(sum(abs(np.log10(text[:,3])-new[:,1]))/len(text[:,0]))
    import material_structure as ms
    remagx = np.loadtxt('optical_profile_text.txt')
    d_d = remagx[:,0]
    delta_r = remagx[:,1]
    idx_d = [i for i in range(0, len(d_d),4)]
    d_b = remagx[:,2]
    idx = [i for i in range(len(d_b)) if d_b[i]!=0]
    d_b = d_b[idx]
    beta_r = remagx[:,3][idx]



    idx_b = [i for i in range(0,len(d_b),2)]
    thickness_r = d_b[idx_b]  # thickness I will be using
    n_r = 1 + np.vectorize(complex)(-delta_r[idx_d], beta_r[idx_b])  # index of refraction

    new_delta = []
    new_beta = []
    new_thickness = []
    for i in range(0,len(delta_r[idx_d])-1):
        new_thickness.append(thickness_r[i])
        new_thickness.append(thickness_r[i])

        new_delta.append(delta_r[idx_d][i])
        new_delta.append(delta_r[idx_d][i+1])

        new_beta.append(beta_r[idx_b][i])
        new_beta.append(beta_r[idx_b][i+1])

    epsilon_r = n_r ** 2  # epsilon for ReMagX

    import Pythonreflectivity as pr
    m = len(epsilon_r)  # number of slabs
    A = pr.Generate_structure(m)  # creates Pythonreflectivity object class

    E = 833.85
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s
    wavelength = h * c / (E * 1e-10)  # wavelength

    m_j = 0
    idx = 0

    for m_i in range(len(epsilon_r)):
        m_i = int(m_i)
        d = thickness_r[m_i] - thickness_r[m_j]  # computes thickness of slab

        # eps = (epsilon[m_i] + epsilon[m_j])/2  # computes the dielectric constant value to use
        eps = epsilon_r[m_i]


        # A[idx].seteps(eps)  # sets dielectric tensor for non-magnetic layer
        A[idx].seteps([eps, eps, eps, 0])

        if idx != 0:
            A[idx].setd(d)  # sets thickness of layer if and only if not substrate layer

        # move onto the next layer
        m_j = m_i
        idx = idx + 1

    qz = new[:, 0]
    plt.show()

    Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
    R = pr.Reflectivity(A, Theta, wavelength, MagneticCutoff=1e-20)  # Computes the reflectivity



    data = np.loadtxt('optical_profile.csv')
    delta = data[:,1]
    beta = data[:,3]
    n = 1 + np.vectorize(complex)(-delta, beta)  # index of refraction

    epsilon = n ** 2
    thickness = data[:,0]
    my_slabs = ALS(delta, beta, precision=1e-1000)

    my_thickness = []
    my_delta = []
    my_beta = []
    # Recent versions of Pythonreflectivity use the susceptibility instead of the dielectric constant
    import Pythonreflectivity as pr
    my_slabs = my_slabs[1:]  # removes first element
    import copy
    copy_slabs = copy.copy(my_slabs)

    m = len(my_slabs)  # number of slabs
    A = pr.Generate_structure(m)  # creates Pythonreflectivity object class

    E = 833.85
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s
    wavelength = h * c / (E * 1e-10)  # wavelength

    m_j=0
    idx = 0

    for m_i in my_slabs:
        m_i = int(m_i)
        d = thickness[m_i] - thickness[m_j]  # computes thickness of slab

        # eps = (epsilon[m_i] + epsilon[m_j])/2  # computes the dielectric constant value to use
        eps = epsilon[m_j]


        # A[idx].seteps(eps)  # sets dielectric tensor for non-magnetic layer
        A[idx].seteps([eps, eps, eps, 0])

        if idx != 0:
            A[idx].setd(d)  # sets thickness of layer if and only if not substrate layer

        # move onto the next layer
        m_j = m_i
        idx = idx + 1


    new_d = []
    new_epsilon = []
    rough = []
    cool_d = [0]
    dto = 0
    for i in range(len(copy_slabs)):
        if i == 0:
            new_d.append(0)
        else:
            new_d.append(A[i].d())
            dto = dto + A[i].d()
            cool_d.append(dto)

        new_epsilon.append(A[i].epsxx())
        rough.append(A[i].sigma())


    new_epsilon = np.array(new_epsilon).imag
    copy_slabs = [int(i) for i in copy_slabs]
    old_eps = np.array(epsilon[copy_slabs]).imag


    Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / np.pi  # initial angle
    Rtemp = pr.Reflectivity(A, Theta, wavelength, MagneticCutoff=1e-20)  # Computes the reflectivity

    plt.figure(5)
    plt.plot(qz, np.log10(Rtemp[0]))
    plt.plot(qz, np.log10(R[0]))
    plt.legend(['GO-RXR','ReMagX'])
    plt.show()


    for i, m_i in enumerate(my_slabs):
        if m_i != 0:
            m_i = int(m_i)
            m_j = int(m_j)


            my_thickness.append(thickness[m_i-1])
            my_thickness.append(thickness[m_i-1])


            eps = delta[m_j]  # computes the dielectric constant value to use

            my_delta.append(eps)
            my_beta.append(beta[m_j])
            my_delta.append(delta[m_i])
            my_beta.append(beta[m_i])

            m_j=m_i


    plt.figure(2)
    #plt.plot(data[:,0], data[:,1])
    plt.stem(my_thickness, my_delta, bottom=0, markerfmt='b-', linefmt='b-')
    plt.stem(new_thickness, new_delta, bottom=0, markerfmt='r-', linefmt='r-')

    plt.legend(['GO-RXR', 'ReMagX'])
    #plt.plot(thickness, delta,'r')
    plt.suptitle('GO-RXR Segmentation')
    #plt.legend(['Segmentation', 'Original'])
    plt.show()

    

    plt.figure(3)
    plt.plot(remagx[:,0],remagx[:,1])
    #plt.plot(thickness, delta)
    #plt.stem(my_thickness, my_delta, linefmt='r-',bottom=0, markerfmt='')
    plt.suptitle('ReMagX Segmentation')
    plt.show()

    plt.figure(4)
    plt.plot(data[:,0], data[:,1])
    #plt.legend(['ReMagX', 'Go-RXR'])
    plt.suptitle('GO-RXR Optical Profile')

    """