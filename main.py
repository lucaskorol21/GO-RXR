import numpy as np
import data_structure as ds
import matplotlib.pyplot as plt
import global_optimization as go
from scipy.signal import *
import copy


if __name__ == '__main__':
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v3.h5"

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_magnetic = sample.density_profile()

    keys = ['Ti','O','C1','C2','Sr','La','Mn']

    density['Mn'] = density['Mn2+'] + density['Mn3+']
    density['Sr'] = copy.copy(density['Sr_f'])
    density['La0.7Sr0.3'] = density['Sr'] + density['La']

    plt.figure(1)
    for key in keys:
        plt.plot(thickness,density[key])

    plt.legend(keys)
    plt.ylabel('Reflectivity, R')
    plt.xlabel('Depth (Angstrom)')
    plt.show()