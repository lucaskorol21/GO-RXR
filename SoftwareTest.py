from material_structure import *
import matplotlib.pyplot as plt
import numpy as np


def plot_density_profile(sample, fig=1):
    thickness, density, density_magnetic = sample.density_profile()
    val = list(density.values())
    mag_val = list(density_magnetic.values())
    check = []
    for key in list(density.keys()):
        if key[-1].isdigit():
            check.append(True)
        else:
            check.append(False)

    plt.figure(fig)
    for idx in range(len(val)):
        if check[idx]:
            plt.plot(thickness, val[idx], ':')
        else:
            plt.plot(thickness, val[idx])

    for idx in range(len(mag_val)):
        plt.plot(thickness, -mag_val[idx], '--')

    center = np.zeros(len(thickness))
    plt.plot(thickness, center, 'k-.', linewidth=2)
    my_legend = list(density.keys())

    for key in list(density_magnetic.keys()):
        my_legend.append('Mag: ' + key)

    plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02, 0.5))
    plt.xlabel('Thickness (Angstrom)')
    plt.ylabel('Density (mol/cm^3)')

if __name__ == "__main__":
    # In this script I will provide some examples on how to setup the sample model,
    # and how to calculate the reflectivity. Your goal is to re-create the provided model
    # (enclosed in the PDF). I will provide further instructions within the script.
    # Good Luck!!!

    # My code relies on using the same momentum transfer as the input data.
    # This will be preset for you when you are testing out the software
    E = 640
    theta_i = 0.1
    theta_f = 89.9
    delta_theta = (theta_f - theta_i) / 301  # sets step size
    Theta = np.arange(theta_i, theta_f + delta_theta, delta_theta)  # Angle array
    qz = (0.001013546247) * E * sin(Theta * pi / 180)  # transforms angle back into momentum transfer
    #------------------------ Sample Setup Examples-------------------------------#

    # ------- Strutural Model
    N = 3
    example1 = slab(N)  # defines the number of slabs in the layer

    # To add a layer call the class method sample.addlayer. This method must be used for all layers defined (N times).
    #  - Required Parameters (num_layer, formula, thickness)
    #       - num_layer: The layer number you are setting. Substrate or bottom most layer has the index 0, and
    #                    the top layer has index N-1. Hint: I find it easiest to start at substrate (bottom most) layer.
    #       - formula: Formula of the elements contained in the specified layer. Use the element symbol and don't
    #                  forget to include the stoichiometry (e.g. LaMnO3).
    #       - thickness: Thickness of the slab in Angstroms.
    #   - Optional Parameters (density, roughness, link)
    #       - density: The density of all the elements combined in units of g/cm^3. If user decides to not
    #                  input their own value, a pre-determined value will be used. Note that this can only be used
    #                  for perovskite materials at the moment.
    #       - roughness: The roughness of the material of the order Angstroms. A roughness of zero corresponds
    #                    to a smooth surface, while a large value corresponds to a rough surface.
    #       - link: Do not worry about this parameter at the moment as it currently has no significance

    example1.addlayer(0,'SrTiO3',50, roughness=2)  # Substrate made up of a strontium titanium oxide
    example1.addlayer(1, 'SrTiO3',2, roughness=2)  # Substrate extends slightly into the thin film

    # Sometimes we want to be able to distinguish between different  polymorphs (forms) of the same element.
    # For example: it can be important to distinguish between the different oxidation states. To incorporate a
    # a polymorph into the your structure you need to use the class method polymorphous. Once you set one element
    # to be a polymorph, you must maintain consistent and set this for the remainder of the layers that include this
    # element.
    #   - Required Inputs (lay, ele, polymorph, poly_ratio)
    #       - lay: The layer the polymorph element is found
    #       - ele: The symbol of the polymorph element in the selected layer
    #       - polymorph: List of the polymorph 'names' (e.g. ['Mn2+','Mn3+'])
    #       - poly_ratio: A list of ratios that determine how much of the polymorph makes up the total elements density
    #           * polymorph = ['Mn2+','Mn3+'], poly_ratio = [0.2,0.8]
    #               ~ Density(Mn2+) = Density(Mn)*0.2
    #               ~ Density(Mn3+) = Density(Mn)*0.8
    #           * note that the sum of all elements in poly_ratio must equal 1
    #   - Optional Inputs (sf)
    #       - sf: List of the scattering factors you want to use for the polymorphs. If no element symbol is chosen
    #             then the polymorph symplot inputted will be used.

    example1.addlayer(2, 'LaMnO3',40, roughness=2)
    example1.polymorphous(2,'Mn', ['Mn2+','Mn3+'], [0.1,0.9], sf=['Mn','Fe'])

    plot_density_profile(example1,fig=1)  # plot the density profile
    # ------- Magnetic Model

    # Now let's implement the same structure, but include a magnetic component.
    # The method used to implement a magnetic layer is called 'magnetization'.
    # The current implementation required you to first define the element as a polymorph before defining
    # the element as magnetic.
    #   - Required Parameters
    #       - lay: Layer the magnetic material is found
    #       - identifier: The symbol of the magnetic element
    #       - density: The density of the magnetism (mol/cm^3).
    #       - sf: Scattering factors of the magnetic components
    #   - Optional Parameters
    #       - phi: Azimuthal angle to the normal of surface (degrees)
    #       - theta: Polar angle to the normal of surface (degrees)

    # Important: Density is in mol/cm^3 unlike the structural definition of g/cm^3.
    # ** Note that the magnetic layer will have the exact same roughness as defined in the structural component.

    N = 3
    example2 = slab(N)

    # Convention is to build one layer at a time
    example2.addlayer(0, 'SrTiO3', 50, roughness=2)  # Substrate

    example2.addlayer(1, 'SrTiO3', 2, roughness=2)  # Substrate extends slightly into the thin film

    example2.addlayer(2, 'LaMnO3', 40, roughness=2)
    example2.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [0.1, 0.9], sf=['Mn', 'Fe'])
    example2.magnetization(2, ['Mn2+','Mn3+'],[0.1,0.1],['Co','Ni'])

    # plot_density_profile(example2,fig=2)  # Plot the density profile

    #-------- Complicated magnetic structure and polymorph structure
    # Sometimes yoy may want a more complicated structure. This can easily be done by dividing a single slab
    # into multiple slabs.

    N = 5
    example3 = slab(N)

    # Convention is to build one layer at a time
    example3.addlayer(0, 'SrTiO3', 50, roughness=2)  # Substrate

    example3.addlayer(1, 'SrTiO3', 2, roughness=1)  # Substrate extends slightly into the thin film

    # Divide slab that is 40 A into a slab of 10, 25, 5 A
    example3.addlayer(2, 'LaMnO3', 10, roughness=2)
    example3.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [0.01, 0.99], sf=['Mn', 'Fe'])
    example3.magnetization(2, ['Mn2+', 'Mn3+'], [0.01, 0.02], ['Co', 'Ni'])

    example3.addlayer(3, 'LaMnO3', 25,roughness=3)
    example3.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    example3.magnetization(3, ['Mn2+', 'Mn3+'], [0.02, 0.05], ['Co', 'Ni'])

    example3.addlayer(4, 'LaMnO3', 5,roughness=2)
    example3.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [0.99, 0.01], sf=['Mn', 'Fe'])
    example3.magnetization(4, ['Mn2+', 'Mn3+'], [0.04, 0.001], ['Co', 'Ni'])

    plot_density_profile(example3, fig=3)  # plot density profile


    # ----------------------------- Reflectivity Computation ---------------------------------#

    # To compute the reflectivity the class method relfectivity must be used.
    #   - Required input parameters (E)
    #       - E: Energy of incoming photon in units of eV
    #       - qz: Momentum Transfer (pre-defined)
    #   - Optional input parameters (precision, s_min)
    #       - precision: precision value that determines maximum slab thickness based on maximum variation
    #       - s_min: value that sets the minimum slab thickness that can be used (A)
    #   - OUTPUT
    #       - qz: Momentum Transfer
    #       - R: Reflectivity
    #           ~ R[0]: s-polarized light
    #           ~ R[1]: p-polarized light
    #           ~ R[2]: left circular polarized light
    #           ~ R[3]: right circular polarized light
    #       - [thickness, plot_t]
    #           ~ thickness: thickness in Angstrom
    #           ~ plot_t: thickness to show how sample is sliced for computation
    #       - [real(epsilon), plot_e]
    #           ~ e: delta (dispersive component) of index of refraction
    #           ~ plot_e: delta to show how sample is sliced for computation
    #       - IGNORE LAST TWO OUTPUT VALUES

    # ignore t and e
    qz, R, t, e = example3.reflectivity(E, qz)  # incoming photon energy of 360 eV

    # Note that when plotting the reflectivity you must take the logarithm to see anything significant
    Rs = np.log10(R[0])  # s-polarized
    Rp = np.log10(R[1])  # p-polarized
    Rl = np.log10(R[2])  # left circular
    Rr = np.log10(R[3])  # right circular
    # Note that if you did not specifiy any magnetism or set the denity to zero for all layer, the function
    # will not return the reflectivity for right or left circular light.

    plt.figure(4)
    plt.plot(qz, Rs)
    plt.suptitle('Reflectivity Spectra')
    plt.xlabel(r"$Momentum\;\;Transfer\;\;(q_z)$")
    plt.ylabel(r'$log_{10}(R)$')

    ####################################################################################################################
    ################# Now your turn! Construct the sample as outlines in the provided word document ####################
    ####################################################################################################################

    E = 0  # Change Energy

    # Ignore this portion
    theta_i = 0.1
    theta_f = 89.9
    delta_theta = (theta_f - theta_i) / 301
    Theta = np.arange(theta_i, theta_f + delta_theta, delta_theta)
    qz = (0.001013546247) * E * sin(Theta * pi / 180)

    ################### Construct Sample - Use this Space to construct model ##########################################

    # N = ?
    # sample = slab(N)




    # plot_density_profile(sample,10)  # use this to plot the density profile of your sample

    ################### Compute reflectivity - Use this space to compute reflectivity ##################################

    # ignore t and e
    # qz, R, t, e = sample.reflectivity(E, qz)
    # R_answer = ?  # select the correct polarization



    ################### Model Check - Uncomment the code below to check your model #####################################

    # Check terminal for printed out results
    """
    reference = np.load('Trial1_check.npy')
    qz_r = reference[0]
    Rr  =reference[1]

    total_diff = sum(abs(Rr-R_answer))
    if total_diff == 0:
        print('Your model is correct!')
    elif total_diff < 1e-4:
        print('Your model is close!')
    else:
        print('Either your model or reflectivity computation is incorrect.')
    """
    plt.show()


