from data_structure import *
from material_structure import *
import matplotlib.pyplot as plt
import numpy as np
from numba import *
from tqdm import tqdm
from time import sleep


if __name__ == '__main__':

    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density =[0.027904,0.027904,0.083712], roughness=[7.58207,False,5.77093])
    sample.addlayer(1, 'SrTiO3', 5, density=[0, 0.027904, 0], roughness=[7.58207, 4.03102, 5.77093])

    sample.addlayer(2,'LaMnO3', 5, density=[0.021798,0.0209,0.084], roughness=[3.77764,2,2],linked_roughness=[False, 0.5, False])
    sample.polymorphous(2,'Mn',['Mn2+', 'Mn3+'], [1,0], sf=['Mn', 'Fe'])
    sample.magnetization(2,['Mn2+','Mn3+'], [0.025,0], ['Co','Ni'])

    sample.addlayer(3,'LaMnO3', 18.8437, density=[0.021798,0.0209,0.084], roughness=[3.77764,2,2])
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 10, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 3, density=[0.025, 0.024, 0.05], roughness=[0.25, 0.25, 2])
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 4, density=[0.025, 0.042, 0.04], roughness=[0.25, 0.25, 2])
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.0053, 0], ['Co', 'Ni'])

    sample.addlayer(7,'CCO', 10.1373, density =[0.05,0.05,0.01], roughness=2, linked_roughness=[3,1.5,False])


    fname = 'Pim10uc.h5'



    f = h5py.File(fname, 'a')
    experiment = f['Experimental_data']
    simulated = f['Simulated_data']

    RS = experiment['Reflectivity_Scan']
    SimR = simulated['Reflectivity_Scan']

    ES = experiment['Energy_Scan']
    SimE = simulated['Energy_Scan']

    # Collects data information to print to terminal
    data = list()
    data_dict = dict()
    sim_dict = dict()

    for Rkey in RS.keys():
        sim = np.array(list(SimR[Rkey])[2])
        new = sim + sim*np.random.rand(len(sim))
        RS[Rkey][2] = new


    for Ekey in ES.keys():
        sim = np.array(list(SimE[Ekey])[2])
        new = sim + sim * np.random.rand(len(sim))
        ES[Ekey][2] = new

    f.close()





