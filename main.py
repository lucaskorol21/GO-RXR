from data_structure import *
from material_structure import *
import matplotlib.pyplot as plt
import numpy as np
from numba import *


if __name__ == '__main__':

    # Example 2: Simple sample creation
    sample = slab(6)  # Initializing four layers
    s = 0.1
    mag_dense = 0.1
    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, density=5.120891853, roughness=2, link=[False, True, True])  # substrate layer
    sample.addlayer(1, 'SrTiO3', 4, density=5.120891853, roughness=2, link=[False, True, True])  # substrate layer

    sample.addlayer(2, 'LaMnO3', 4, density=6.8195658, roughness=2)
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0, 0], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 30, density=6.8195658, roughness=2)
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [mag_dense, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 4, density=6.8195658, roughness=2)
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 4, density=6.8195658, roughness=2)
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0, 0], ['Co', 'Ni'])

    # sample.addlayer(4, 'CCC', 4, density = 0, roughness = 2)

    #sample.plot_density_profile()
    #start = time()
    thickness, density, density_magnetic = sample.density_profile(step=0.1)
    E = 642.2  # eV
    
    elements = density.keys()
    sf = dict()
    elements_mag = density_magnetic.keys()
    sfm = dict()
    # Non-Magnetic Scattering Factor
    for e in sample.find_sf[0].keys():
        sf[e] = find_form_factor(sample.find_sf[0][e], E, False)
    # Magnetic Scattering Factor
    for em in sample.find_sf[1].keys():
        sfm[em] = find_form_factor(sample.find_sf[1][em], E, True)

    delta, beta = index_of_refraction(density, sf, E)

    plt.figure(2)
    plt.plot(thickness, delta, thickness, beta)
    plt.suptitle('Optical Profile')
    plt.xlabel('Thickness')
    plt.ylabel('Profile')
    plt.legend(['delta', 'beta'])

    h = 4.1257e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s
    y = h * c / (E * 1e-10)  # wavelength m

    F = np.loadtxt('test_example.txt')
    qz = F[:, 0]
    p1 = 1e-6
    p2 = 1e-7

    qz, Rtemp = sample.reflectivity(E, qz, precision=1e-15, s_min=0.1)  # baseline

    qz1, R1temp = sample.reflectivity(E, qz, precision=p1, s_min=0.1)


    qz2, R2temp = sample.reflectivity(E, qz, precision=p2, s_min=0.1)


    # R = np.log10(Rtemp['LC'])
    # R1 = np.log10(R1temp['LC'])
    # R2 = np.log10(R2temp['LC'])

    R = Rtemp['AC']
    R1 = R1temp['AC']
    R2 = R2temp['AC']

    plt.figure(4)
    plt.plot(qz, R, 'k-')
    plt.plot(qz1, R1, 'b--')
    plt.plot(qz2, R2, 'r--')
    plt.legend(['baseline', str(p1), str(p2)])
    # plt.yscale("log")
    plt.xlabel('qz')
    plt.ylabel('Reflectivity ' + "$(log_{10}(R))$")
    plt.title('ReMagX vs. Python Script (800 eV)')

    diff_1 = abs(R - R1) / abs(R + R1)
    figure(5)
    plt.plot(qz1, diff_1)
    plt.suptitle("Difference in Spectra: " + str(p1))
    plt.xlabel("Thickness (Angstroms)")
    plt.ylabel("Percent Difference")

    diff_2 = abs(R - R2) / abs(R + R2)
    figure(6)
    plt.plot(qz1, diff_2)
    plt.suptitle("Difference in Spectra: " + str(p2))
    plt.xlabel("Thickness (Angstrom)")
    plt.ylabel("Percent Difference")

    max1 = max(abs(R - R1))
    max2 = max(abs(R - R2))
    A1 = simpson(abs(R - R1), qz)
    A2 = simpson(abs(R - R2), qz)

    #print()
    #print()
    #print(tabulate([[p1, max1, A1], [p2, max2, A2]], headers=['Precision', 'Maximum', 'Total Area']))

    q = F[:, 0]
    I = F[:, 1]
    # I = np.log10(I)
    plt.figure(55)
    plt.plot(q, I, 'k')
    plt.plot(qz, R1, 'r--')
    plt.suptitle('Zak Formalism: Asymmetry 642.2 eV (rho=0.5) ')
    plt.xlabel('qz')
    plt.ylabel('Reflectivity ' + "$(log_{10}(R))$")
    plt.legend(['ReMagX', 'Lucas'])

    diff_3 = abs(I - R1) / abs(I + R1)
    plt.figure(56)
    plt.plot(q, diff_3, 'k')
    plt.suptitle('Percent Difference: precision = ' + str(0))
    plt.xlabel('qz')
    plt.ylabel('Percent Difference')
    plt.legend(['ReMagX', 'Lucas'])
    #plt.show()

    hello = np.loadtxt('energy_test.txt')
    E = hello[:,0]

    R = np.array(R)

    start = time()
    E, R = sample.energy_scan(5.0, E)
    end = time()
    print(end-start)
    R = R['S']

    #R = R['P']
    plt.figure(89)
    plt.plot(E,R,'.')
    plt.plot(E, hello[:,1])
    plt.legend(['Simulation', 'Data'])
    plt.show()
    #fname = "FGT-1L.all"
    #Sscan, Sinfo, sample1 = ReadData(fname)
    #selectScan(Sinfo, Sscan, sample)

