import numpy as np
import data_structure as ds
import matplotlib.pyplot as plt
import global_optimization as go
from scipy.signal import *
import copy
from scipy import interpolate


from scipy.fft import fft, fftfreq, fftshift, ifft




def noise_removal(qz, R, s=1, k=3):
    tck = interpolate.splrep(qz,R, s=s, k=k)
    return interpolate.splev(qz, tck, der=0)





if __name__ == '__main__':
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v5.h5"

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_magnetic = sample.density_profile()

    density['Mn'] = density['Mn2+'] + density['Mn3+']
    density['Sr/La'] = density['La'] + density['Sr']

    keys = ['Ti', 'O', 'C1','C2','Sr/La','Mn']

    plt.figure()
    for name in keys:
        plt.plot(thickness,density[name])

    plt.legend(keys)
    # plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('Density (mol/cm^3)')
    plt.show()

    print(density.keys())



