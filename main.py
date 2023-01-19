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
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_magnetic = sample.density_profile()

    density['Mn'] = density['Mn2+'] + density['Mn3+']
    density['Sr/La'] = density['La'] + density['Sr']

    #keys = ['Ti', 'O','Sr','La','Mn2+','Mn3+']
    keys = ['Ti', 'O','Sr/La','Mn']

    plt.figure()
    for name in keys:
        plt.plot(thickness,density[name])

    plt.legend(keys)
    #plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('Density (mol/cm^3)')
    plt.show()

    electron = {'Sr':2,'Ti':4,'La':3,'Mn2+':2,'Mn3+':3.3,'O':-2}
    d = 40

    electric_conservation = np.zeros(len(density['Sr']))
    idx = np.argwhere(thickness<d)

    for key in list(electron.keys()):
        electric_conservation[idx] = electric_conservation[idx] - electron[key]*density[key][idx]

    electric_conservation[idx] = electric_conservation[idx]/density['O'][idx]
    plt.figure()
    plt.plot(thickness, electric_conservation)
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('# electrons per oxygen atom')
    plt.show()


