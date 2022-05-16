from math import sin
import numpy as np
from scipy import interpolate
from scipy import optimize
import Pythonreflectivity as pr
from material_model import *






if __name__ == "__main__":

    E = 640.2  # eV
    mag_optical_profile = np.loadtxt('mag_optical_profile')
    optical_profile = np.loadtxt('optical_constant')
    spectra = np.loadtxt('test_example.txt')

    # optical profile
    t = optical_profile[:,0]  # thickness
    delta = optical_profile[:, 1]  # delta
    beta = optical_profile[:, 3]  # beta

    # magnetic optical profile
    tm = mag_optical_profile[:, 0]  # thickness
    delta_m = mag_optical_profile[:, 1]  # delta
    beta_m = mag_optical_profile[:, 3]  # beta

    n = 1 + np.vectorize(complex)(-delta, beta)
    #epsilon = 1 + np.vectorize(complex)(-2*delta,2*beta)
    epsilon = n**2
    Q = np.vectorize(complex)(beta_m, delta_m)
    epsilon_mag = Q*epsilon*2

    m = len(t)-1

    A = pr.Generate_structure(m)  # initializes structure

    for idx in range(m):
        d = t[idx+1] - t[idx]
        eps = (epsilon[idx+1]+epsilon[idx])/2
        #eps = epsilon[idx]
        eps_mag = (epsilon_mag[idx+1]+epsilon_mag[idx])/2
        #eps_mag = epsilon_mag[idx]

        A[idx].setmag('y')
        A[idx].seteps([eps,eps,eps,eps_mag])


        A[idx].setd(d)


    qz = spectra[:,0]
    qi = qz[0]
    qf = qz[-1]
    length_qz = len(qz)



    h = 4.1257e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s
    wavelength = h * c / (E * 1e-10)  # wavelength m

    theta_i = arcsin(qi / E / (0.001013546247)) * 180 / pi  # initial angle
    theta_f = arcsin(qf / E / (0.001013546247)) * 180 / pi  # final angle in interval
    delta_theta = (theta_f - theta_i) / (length_qz)  # sets step size
    Theta = np.arange(theta_i, theta_f + delta_theta, delta_theta)  # Angle array



    R = pr.Reflectivity(A, Theta, wavelength,MagneticCutoff=1e-10)
    Rs = R[0]
    Rp = R[1]
    Rl = R[2]
    Rr = R[3]
    Ra = (Rl-Rr)/(Rl+Rr)  # asymmetry

    Rs = np.log10(Rs)
    Rp = np.log10(Rp)
    Rl = np.log10(Rl)
    Rr = np.log10(Rr)

    Re = np.log10(spectra[:, 1])  # ReMagX Spectra
    Rea = spectra[:,1]
    qz_temp = (0.001013546247) * E * sin(Theta * pi / 180)  # transforms angle back into momentum transfer

    # Interpolation
    itr = interpolate.splrep(qz_temp, Ra)
    R_fit = interpolate.splev(qz, itr)
    R_diff = abs(R_fit-Rea)

    fig, (ax1, ax2) = plt.subplots(2,1, sharex=True)
    fig.suptitle(r"$No\;\;Magnetic\;\;Cutoff:\;\;Asymmetry\;\;for\;\;E = 640.2eV\;\;and\;\;\rho=0.001\;\;mol/cm^3$")
    ax1.plot(qz_temp, Ra,'k')
    ax1.plot(qz, Rea,'r--')
    ax1.legend(['Lucas','ReMagX'])
    ax1.set_title("Reflectivity Spectra")
    ax1.set_ylabel(r"$log_{10}(R)y$")

    ax2.plot(qz, R_diff)
    ax2.set_title("Difference between ReMagX and Lucas Code")
    ax2.set_xlabel(r"$Momentum\;\;Transfer\;\;(q_z)$")
    ax2.set_ylabel(r"$Difference log_{10}(R_1/R_2)$")
    plt.subplots_adjust(hspace=0.3,top=0.8)


    # Testing Optical profile Ratios

    # Example 2: Simple sample creation
    sample = slab(2)  # Initializing four layers
    s = 0.1
    mag_dense = 0.0001
    # Substrate Layer
    # Link: Ti-->Mn and O-->O
    sample.addlayer(0, 'SrTiO3', 50, density=0, roughness=4, link=[False, True, True])  # substrate layer

    sample.addlayer(1, 'LaMnO3', 30, density=6.8195658, roughness=4)
    sample.polymorphous(1, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [mag_dense, 0], ['Co', 'Ni'])


    thickness, density, density_magnetic = sample.density_profile(step=0.1)

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

    d, b = index_of_refraction(density, sf,E)  # calculates dielectric constant for structural component

    dm, bm = magnetic_optical_constant(density_magnetic, sfm,E)  # calculates dielectric constant for magnetic component
    dm = dm*(-1)
    bm = bm*(-1)

    itr_d = interpolate.splrep(thickness, d)
    d_fit = interpolate.splev(t, itr_d)
    d_diff = abs(d_fit-delta)
    d_rat = d_fit/delta

    itr_b = interpolate.splrep(thickness, b)
    b_fit = interpolate.splev(t, itr_b)
    b_diff = abs(b_fit - beta)
    b_rat = b_fit/beta

    itr_dm = interpolate.splrep(thickness, dm)
    dm_fit = interpolate.splev(t, itr_dm)
    dm_diff = abs(dm_fit - delta_m)/abs(dm_fit+delta_m)
    dm_rat = dm_fit/delta_m

    itr_bm = interpolate.splrep(thickness, bm)
    bm_fit = interpolate.splev(t, itr_bm)
    bm_diff = abs(bm_fit - beta_m)
    bm_rat = bm_fit/beta_m

    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
    fig.suptitle(r"$Optical\;\;Profile\;\;Comparison\;\;(\delta_m)\;\;for\;\;\rho=0.01\;\;mol/cm^3$")
    ax1.plot(t, d_rat)
    ax1.set_title("Ratio")
    ax1.set_ylabel(r"$\delta_m^{L}/\delta_m^{R}$", fontsize=14)

    ax2.plot(t, dm_diff)
    ax2.set_title("Difference")
    ax2.set_xlabel("Thickness (A)")
    ax2.set_ylabel(r"$|\delta_m^{L}=\delta_m^{R}|$", fontsize=14)
    plt.subplots_adjust(hspace=0.3, top=0.85)
    plt.show()

    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)
    fig.suptitle('Optical Profile Comparison')
    ax1.plot(t, d_fit,'k')
    ax1.plot(t,delta,'r--')
    ax1.set_title(r"$\delta$")
    ax1.legend(['Lucas','ReMagX'],loc='upper right', fontsize=8)

    ax2.plot(t,b_fit, 'k')
    ax2.plot(t, beta, 'r--')
    ax2.set_title(r"$\beta$")

    ax3.plot(t, dm_fit,'k')
    ax3.plot(t, delta_m,'r--')
    ax3.set_title(r"$\delta_m$")

    ax4.plot(t, bm_fit,'k')
    ax4.plot(t, beta_m,'r--')
    ax4.set_title(r"$\beta_m$")
    plt.subplots_adjust(hspace=0.3, wspace=0.3)


    plt.show()
