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

def hamming(v, No, tua):
    omegaC = 2*np.pi*v/No # cutoff frequency
    k = np.arange(-tua,tua,1)
    sink = np.sinc(omegaC*k)*omegaC/np.pi
    window = 0.54+0.46*np.cos(k/tua)
    filter = sink*window

    return filter

def Fourier_noise_removal(qz, R):
    from scipy.interpolate import UnivariateSpline
    from scipy.fft import fft, fftfreq, fftshift, ifft, ifftshift

    # convert to log(R)
    y_old = np.log10(R)
    y = np.log10(R)
    x = qz
    spl = UnivariateSpline(x,y)
    plt.figure()
    plt.plot(x,y-spl(x))

    # perform spline subtration
    spl = UnivariateSpline(x,y)
    y_sub = spl(x)  # spline for subtraction

    y = y - y_sub  # spline subtraction

    # Fourier
    N = len(y)  # number of sample points
    T = np.average(np.diff(x))  # sample spacing

    # I need to make sure that each data point is equally spaced (use some kind of interpolation for this!)

    yf = fft(y)
    xf = fftfreq(N,T)
    xf = fftshift(xf)
    yplot = fftshift(yf)
    plt.figure()
    plt.plot(xf, 1.0/N*np.abs(yplot))
    plt.show()

    filter = np.zeros(len(xf))
    val = 0.2
    for idx in range(len(filter)):
        if xf[idx] > -val and xf[idx] < val:
            filter[idx] = 1

    yfiltered = yplot*filter
    plt.figure()
    plt.plot(xf, 1.0/N*np.abs(yfiltered))
    plt.show()

    yfiltered = ifftshift(yfiltered)

    ytemp = ifft(yfiltered)
    print(len(ytemp))
    plt.figure()
    plt.plot(x,ytemp)
    plt.show()

    ynew = ytemp+y_sub
    plt.figure()
    plt.plot(x, y_old)
    plt.plot(x, ynew)
    plt.legend(['Data','Filtered'])
    plt.show()



    # transform back to original data set

    return x,y



if __name__ == '__main__':
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-300K_v3.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)

    thickness, density, density_magnetic = sample.density_profile()
    print(density_magnetic)
    sample.plot_density_profile()
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
    d = 50

    electric_conservation = np.zeros(len(density['Sr']))
    idx = np.argwhere(thickness<d)

    for key in list(electron.keys()):
        electric_conservation[idx] = electric_conservation[idx] - electron[key]*density[key][idx]

    electric_conservation[idx] = electric_conservation[idx]/density['O'][idx]
    plt.figure()
    plt.plot(thickness, electric_conservation)
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('# electrons per oxygen atom')
    #plt.show()


    f1 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/150.txt"
    f2 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/300.txt"

    sample_150 = np.loadtxt(f1)
    sample_300 = np.loadtxt(f2)

    plt.figure()
    plt.suptitle('Magnetic Profile for E2 Sample at 150 and 300K')
    plt.plot(sample_150[:,0], sample_150[:,1]*(-1))
    plt.plot(sample_300[:,0], sample_300[:,1]*(-1))
    plt.xlabel('Thickness (angstrom)')
    plt.ylabel('Density')
    plt.legend(['150K', '300K'])
    plt.show()

    d1 = -6
    d2 = 47.6
    m = len(sample_300[:,1])
    my_division = np.zeros(m)

    for i in range(m):
        if sample_150[i,0] > d1 and sample_150[i,0] < d2:
            my_division[i] = sample_300[i,1]/sample_150[i,1]

    plt.figure()
    plt.suptitle('Ratio of E2 Sample at 150 and 300K')
    plt.plot(sample_300[:,0], my_division)
    plt.xlabel('Thickness (angstroms)')
    plt.ylabel('300K/150K')
    plt.show()

    print(sum(sample_300[:,1])/sum(sample_150[:,1]))