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
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    #sample = ms.slab(2)
    #sample.addlayer(0,'SrTiO3',50,roughness=0)
    #sample.addlayer(1,'LaMnO3',246,roughness=0)
    #sample.energy_shift()


    struct_names, mag_names = mm._use_given_ff(os.getcwd())  # look for form factors in directory
    #self._sampleWidget.struct_ff = struct_names
    #self._sampleWidget.mag_ff = mag_names

    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    name = "33_834.59_S"
    qz_temp = data_dict[name]['Data'][0]
    E = data_dict[name]['Energy']

    qz_min = qz_temp[0]
    qz_max = qz_temp[-1]

    qz = np.linspace(qz_min, qz_max, num=15000)

    qz, R = sample.reflectivity(E, qz)

    R = R['S']

    dat = np.array([qz,R])
    dat = np.transpose(dat)

    np.savetxt('834.59_15000.txt',dat)
    """
    thickness, density, density_magnetic = sample.density_profile()

    E = 834 #eV
    order = 5
    b = 0
    theta_i = 0.01
    theta_f = 60.01

    qz_i = np.sin(theta_i*np.pi/180)*E*(0.001013546247)
    qz_f = np.sin(theta_f*np.pi/180)*E*(0.001013546247)

    qz = np.linspace(qz_i,qz_f, num=10001)
    theta = np.arcsin(qz/(E*(0.001013546247)))*180/np.pi

    qz, R = sample.reflectivity(E, qz)
    R = R['S']



    theta_c = find_cga(theta,R)
    qz, R = index_correction(theta, theta_c, E,R)
    #R = np.log10(R)

    spl = UnivariateSpline(qz,np.log10(R), k=order)


    R4 = np.multiply(R, np.power(qz, 4))
    Rs = np.log10(R)-spl(qz)

    x,y,N,T = fourier_transform(qz,Rs)


    from scipy.signal import find_peaks

    peaks = find_peaks(np.abs(y), height=0)
    peaks = peaks[0]
    print(x[peaks])

    plt.figure(1)
    plt.plot(qz, Rs)
    #plt.plot(qz,R4*1e6)

    plt.figure(2)
    plt.plot(qz,np.log10(R))
    plt.plot(qz,spl(qz))
    #plt.plot(qz,spl(qz))


    plt.figure(3)
    plt.plot(x,1.0/N*np.abs(y))
    plt.xlabel('Position (angstroms)')
    plt.show()
    """

