from data_structure import *
from material_structure import *
from GlobalOptimization import *
import matplotlib.pyplot as plt
import numpy as np
from numba import *
from scipy import interpolate
from scipy.fft import fft, fftfreq,fftshift, rfft, irfft
from scipy.interpolate import UnivariateSpline
from scipy import signal
import time




if __name__ == '__main__':
    # Define the sample model
    sample1 = slab(2)
    sample1.addlayer(0, 'SrTiO3', 50,roughness=4)
    sample1.addlayer(1, 'LaMnO3', [30,20,30],roughness=1.5, linked_roughness=1)

    sample2 = slab(2)
    sample2.addlayer(0, 'SrTiO3', 50, density=4.8579022, roughness=5.97722238)
    sample2.addlayer(1, 'LaMnO3', 36.64782587,density=6.77792568, roughness=5.06593207, linked_roughness=0.62078224)

    sample1.plot_density_profile(1)
    plt.xlim([-25,50])
    plt.show()

    fname = 'Pim10uc.h5'
    info, data_dict, sim_dict = ReadDataHDF5(fname)
    mydata = data_dict[info[71][2]]
    print(info[71][2])
    E = mydata['Data'][3]
    R = mydata['Data'][2]
    Theta = mydata['Angle']
    E, Rn = sample1.reflectivity(Theta, E)
    Rn = Rn['S']

    plt.figure(2)
    plt.plot(E, R)
    #plt.plot(qz, Rn)
    plt.ylabel('Reflection Intensity (R)')
    plt.xlabel(r'Energy, E (eV)')
    plt.yscale('log')
    plt.show()

    """
    # save new sample model to current file
    fname = 'Pim10uc.h5'
    #WriteSampleHDF5(fname, sample)

    scans, parameters, bounds = getGlobOptParams(fname)

    #parameters = [[1, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
    #              [2, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
    #              [3, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
    #              [4, 'STRUCTURAL', 'COMPOUND', 'THICKNESS']]

    #lw = [3.5, 3.5, 17.3, 8.5]
    #up = [6.5, 6.5, 19.8, 11.5]
    #bounds = list(zip(lw, up))
    #scans = [1, 2, 3, 4, 5, 6]

    # determines the bounds for the scans
    sBounds = [[(0.1, 0.8)],
               [(0.1, 0.3), (0.3, 0.5), (0.6, 0.8)],
               [(0.1, 0.6), (0.7, 0.8)],
               [(0.1, 0.5)],
               [(0.2, 0.6), (0.6, 0.8)],
               [(0.1, 0.8)]]

    # Determines the weights you want to use for each bound
    sWeights = [[1],
                [1, 0.2, 0.5],
                [1, 0.1],
                [0.5],
                [1, 0.8],
                [0.7]]

    scanBounds = createBoundsDatatype(fname, scans, sBounds, sWeights=sWeights)

    start = time.time()
    x, fun = differential_evolution(fname, scans, parameters, bounds, scanBounds, mIter=10, display=True,
                                    tolerance=1e-6)
    end = time.time()
    print(end - start)

    comparisonScanPlots(fname, x, parameters, scans)
    # Show the results for the global optimization
    # compare with previous version
    # allow user to use current sample model for next optimization
    # give user option to save new sample

    """


