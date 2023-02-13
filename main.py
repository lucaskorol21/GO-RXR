import numpy as np
import data_structure as ds
import matplotlib.pyplot as plt
import global_optimization as go
from scipy.signal import *
import copy
from scipy import interpolate
import material_model as mm
import os


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

def fourier_transform(x,y):
    from scipy.fft import fft, fftfreq, fftshift, ifft, ifftshift

    # Fourier
    print(np.diff(x))
    N = len(y)  # number of sample points
    T = np.average(np.diff(x))  # sample spacing

    # I need to make sure that each data point is equally spaced (use some kind of interpolation for this!)

    yf = fft(y)
    xf = fftfreq(N, T)
    xf = fftshift(xf)
    yf = fftshift(yf)
    return xf,yf, N, T



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

def adaptive_running_average(R, window):
    pass

def index_correction(theta, theta_c, E, R):
    idx = [i for i in range(len(theta)) if theta[i]>=theta_c]
    theta = theta[idx]
    R = R[idx]
    qz = E*(0.001013546247)*np.sqrt(np.power(np.sin(theta*np.pi/180),2)-np.power(np.sin(theta_c*np.pi/180),2))
    return qz, R

def find_cga(theta,R):
    from scipy.interpolate import interp1d

    f = interp1d(R,theta)
    return f(0.5)


if __name__ == '__main__':
    from scipy.interpolate import UnivariateSpline
    import material_structure as ms
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9-test.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    struct_names, mag_names = mm._use_given_ff(os.getcwd())  # look for form factors in directory

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


