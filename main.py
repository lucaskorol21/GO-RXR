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
    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density =[0.027904,0.027904,0.083712], roughness=[7.58207,False,5.77093])
    sample.addlayer(1, 'SrTiO3', 5, density=[0, 0.027904, 0], roughness=[7.58207, 4.03102, 5.77093])

    sample.addlayer(2,'LaMnO3', 5, density=[0.021798,0.0209,0.084], roughness=[3.77764,2,2],linked_roughness=[False, 0.5, False])
    sample.polymorphous(2,'Mn',['Mn2+', 'Mn3+'], [1,0], sf=['Mn', 'Fe'])
    sample.magnetization(2,['Mn2+','Mn3+'], [0.025,0], ['Co','Ni'])

    sample.addlayer(3,'LaMnO3', 18.8437, density=[0.021798,0.0209,0.084], roughness=[3.77764,2,2])
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 10, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 3, density=[0.025, 0.024, 0.05], roughness=[0.25, 0.25, 2])
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 4, density=[0.025, 0.042, 0.04], roughness=[0.25, 0.25, 2])
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.0053, 0], ['Co', 'Ni'])

    sample.addlayer(7,'CCO', 10.1373, density =[0.05,0.05,0.01], roughness=2, linked_roughness=[3,1.5,False])
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
    print('hello'.split())



