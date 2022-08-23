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
    sample = slab(2)

    sample.addlayer(0, 'SrTiO3', 50, density = [0.27, 0.29,0.84],roughness=4)
    sample.addlayer(1, 'LaMnO3', 30, density=[0.26,0.28,0.82],roughness=1.5, linked_roughness=0.5)

    sample.plot_density_profile(1)
    plt.xlim([-20,40])
    plt.show()

    thickness, density, density_magnetic = sample.density_profile(step=0.1)  # Computes the density profile
    E = 450
    sf = {'Sr': array([24.97849388, 17.13689607]), 'Ti': array([-11.6888725 ,   2.07084738]), 'O': array([5.49082694, 0.30418793]), 'La': array([24.56415562,  9.75841986]), 'Mn': array([14.23069285,  2.83051044])}


    delta, beta = index_of_refraction(density, sf, E)  # calculates dielectric constant for structural component
    plt.figure(2)
    plt.plot(thickness, delta)
    plt.plot(thickness, beta)
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel(r'$\delta, \;\; \beta$')
    plt.legend([r'$\delta$',r'$\beta$'])
    plt.xlim([-20,40])
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


