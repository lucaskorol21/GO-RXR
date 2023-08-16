import os
import pickle
import material_structure as ms
import matplotlib.pyplot as plt
import material_model as mm
import numpy as np


def saveTest(data, fname):

    my_dir = os.getcwd()

    filename = my_dir + '/test_data/' + fname
    with open(filename, 'wb') as file:
        pickle.dump(data, file)


def simple_model():
    sample = ms.slab(5)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(1, 'LaMnO3', 10, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])
    sample.addlayer(2, 'CaTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[5, 2, 0.5])
    sample.addlayer(3, 'LaAlO3', 30, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(4, 'AlYO2', 50, density=4.5, roughness=[3, 2, 0.77])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'simple_density.pkl')

def var_model():
    sample = ms.slab(2)
    sample.addlayer(0, 'ABO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(1, 'ABO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])

    sample.polymorphous(0, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(0, 'B', ['Ti', 'Mn2', 'Mn3'], [1, 0, 0], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(1, 'A', ['Sr', 'La'], [0.25, 0.75], sf=['Sr', 'La'])
    sample.polymorphous(1, 'B', ['Ti', 'Mn2', 'Mn3'], [0.2, 0.35, 0.45], sf=['Ti', 'Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'var_density.pkl')

def surfaceImpurity_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 2, 2.5])
    sample.addlayer(1, 'LaMnO3', 10, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3])
    sample.addlayer(2, 'CCO', [15, 10, 4], density=[0, 0.08, 0.05], roughness=[2, 3, 0.75])

    thickness, density, mag_density = sample.density_profile()

    test_case = density
    test_case['Thickness'] = thickness

    saveTest(test_case, 'surfaceImpurity_density.pkl')

def unitCell_model():
    sample = ms.slab(6)
    sample.addlayer(0, 'ABO3', 50, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(1, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(2, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(3, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(4, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)
    sample.addlayer(5, 'ABO3', 3.97, density=[0.028, 0.028, 0.084], roughness=0)

    sample.polymorphous(0, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(0, 'B', ['Ti', 'Mn2', 'Mn3'], [1, 0, 0], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(1, 'A', ['Sr', 'La'], [0.75, 0.25], sf=['Sr', 'La'])
    sample.polymorphous(1, 'B', ['Ti', 'Mn2', 'Mn3'], [0.2, 0.35, 0.45], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(2, 'A', ['Sr', 'La'], [0.5, 0.5], sf=['Sr', 'La'])
    sample.polymorphous(2, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0.25, 0.75], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(3, 'A', ['Sr', 'La'], [0.25, 0.75], sf=['Sr', 'La'])
    sample.polymorphous(3, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0, 1], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(4, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(4, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0, 1], sf=['Ti', 'Fe', 'Mn'])
    sample.polymorphous(5, 'A', ['Sr', 'La'], [0, 1], sf=['Sr', 'La'])
    sample.polymorphous(5, 'B', ['Ti', 'Mn2', 'Mn3'], [0, 0, 1], sf=['Ti', 'Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    test_case = density
    test_case['Thickness'] = thickness

    saveTest(test_case, 'unitCell_density.pkl')

def linkedRough_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],linked_roughness=[False,0.25,False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5], linked_roughness=[4,4.5,False])

    sample.polymorphous(1, 'Mn', ['Mn2','Mn3'], [0.5,0.5], sf=['Fe', 'Mn'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'linked_density.pkl')

def magnetic_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, 4.5, False])

    sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

    sample.magnetization(1, ['Mn2','Mn3'], [0.015,0.01], ['Co','Ni'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness
    for key in mag_density:
        data['Mag:'+key] = mag_density[key]

    saveTest(data, 'mag_density.pkl')

def dummy_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnX0', 15, density=4, roughness=[0, 1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiQ0', 20, density=3.5, roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, 4.5, False])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'dummy_density.pkl')

def negative_model():
    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[-0.028, 0.028, 0.084], roughness=[1.5, 5, -2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, -0.084], roughness=[0, -1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, -0.028, 0.084], roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, -4.5, False])

    sample.magnetization(1, ['Mn'], [-0.015], ['Co'])

    thickness, density, mag_density = sample.density_profile()

    data = density
    data['Thickness'] = thickness

    saveTest(data, 'negative_density.pkl')

if __name__ == "__main__":

    sample = ms.slab(3)
    sample.addlayer(0, 'SrTiO3', 50, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5])
    sample.addlayer(1, 'LaMnO3', 15, density=[0.028, 0.028, 0.084], roughness=[0, 1, 3],
                    linked_roughness=[False, 0.25, False])
    sample.addlayer(2, 'SrTiO3', 20, density=[0.028, 0.028, 0.084], roughness=[1.5, 5, 2.5],
                    linked_roughness=[4, 4.5, False])

    sample.polymorphous(1, 'Mn', ['Mn2', 'Mn3'], [0.5, 0.5], sf=['Fe', 'Mn'])

    sample.magnetization(1, ['Mn2', 'Mn3'], [0.015, 0.01], ['Co', 'Ni'])

    sample.energy_shift()

    energy = np.array([641,642.1,650.5,653])

    thickness, density, density_magnetic = sample.density_profile(step=0.1)  # Computes the density profile
    # Magnetic Scattering Factor
    sfm = dict()
    sf = dict()


    for e in sample.find_sf[0].keys():
        sf[e] = mm.find_form_factor(sample.find_sf[0][e], energy, False)
    # Magnetic Scattering Factor
    for em in sample.find_sf[1].keys():
        sfm[em] = mm.find_form_factor(sample.find_sf[1][em], energy, True)


    d_len = len(thickness)
    delta, beta = mm.IoR(density, sf, energy)  # gets absorptive and dispersive components of refractive index

    delta_m, beta_m = mm.MOC(density_magnetic, sfm, energy,
                          d_len)  # absorptive and dispersive components for magnetic components

    temp = np.zeros((16,751))
    print(len(delta[0,:]), len(temp[0,:]))

    temp[0,:] = delta[0,:]
    temp[1, :] = delta[1, :]
    temp[2, :] = delta[2, :]
    temp[3, :] = delta[3, :]

    temp[4, :] = beta[0, :]
    temp[5, :] = beta[1, :]
    temp[6, :] = beta[2, :]
    temp[7, :] = beta[3, :]

    temp[8, :] = delta_m[0, :]
    temp[9, :] = delta_m[1, :]
    temp[10, :] = delta_m[2, :]
    temp[11, :] = delta_m[3, :]

    temp[12, :] = beta_m[0, :]
    temp[13, :] = beta_m[1, :]
    temp[14, :] = beta_m[2, :]
    temp[15, :] = beta_m[3, :]


    #my_data = np.array([delta, beta, delta_m, beta_m])
    #fname = os.getcwd() + '/test_data/optical_energy.txt'
    #np.savetxt(fname, temp)

    form_factors = ['Co', 'Ni']
    energy = [400, 600, 625, 641.51, 648.25, 800]

    solution = [np.array([0., 0.]), np.array([0., 0.]), np.array([-0.00681891, 0.0]),
                np.array([0.30686376, 2.47434212]), np.array([0.24675064, -0.34005255]), np.array([0., 0.]),
                np.array([0., 0.]), np.array([0., 0.]), np.array([0.1187884, 0.00154338]),
                np.array([-4.77536493, -3.46825118]), np.array([-0.92991243, -0.23784557]), np.array([0., 0.])]
    test_case = []
    for ff in form_factors:
        for E in energy:
            test_case.append(mm.find_form_factor(ff, E, True))

    my_sum = 0
    for idx in range(len(test_case)):
        my_sum += abs(sum(solution[idx] - test_case[idx]))

    print(my_sum)


