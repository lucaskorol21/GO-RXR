from material_model import *
from material_structure import *
import numpy as np

def check_results(R_guess):
    E = 625  # Change Energy

    # Ignore this portion
    theta_i = 0.1
    theta_f = 89.9
    delta_theta = (theta_f - theta_i) / 301
    Theta = np.arange(theta_i, theta_f + delta_theta, delta_theta)
    qz = (0.001013546247) * E * sin(Theta * pi / 180)

    ################### Construct Sample - Use this Space to construct model
    N = 8
    sample = slab(N)

    # Convention is to build one layer at a time
    sample.addlayer(0, 'SrTiO3', 50, roughness=2)  # Substrate

    sample.addlayer(1, 'SrTiO3', 5, density=5, roughness=1)  # Substrate extends slightly into the thin film

    # Divide slab that is 40 A into a slab of 10, 25, 5 A
    sample.addlayer(2, 'LaMnO3', 10, roughness=1.5)
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [0.01, 0.99], sf=['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0.01, 0.02], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 25, roughness=5)
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.001, 0.025], ['Co', 'Ni'])

    sample.addlayer(4, 'LaAlO3', 15, roughness=3)

    sample.addlayer(5, 'LaMnO3', 7, roughness=2)
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [0.99, 0.01], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.04, 0.001], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 9, roughness=0.1)
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.6, 0.4], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.03, 0.01], ['Co', 'Ni'])

    sample.addlayer(7, 'LaMnO3', 5, roughness=2)
    sample.polymorphous(7, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(7, ['Mn2+', 'Mn3+'], [0.04, 0.001], ['Co', 'Ni'])


    ################### Compute reflectivity - Use this space to compute reflectivity

    # ignore t and e
    qz, R, t, e = sample.reflectivity(E, qz)  # incoming photon energy of 360 eV
    Rl = np.log10(R[2])

    ################### Model Check - Use this space to check your model with expected output
    okay = sum(Rl-R_guess)
    return okay