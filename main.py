import numpy as np
import data_structure as ds
import matplotlib.pyplot as plt


def createFF(directory, filename, my_dict, source, date):
    import os


    f = os.path.join(directory, filename)
    if not (filename.startswith('README')):
        name = f.split('.')[0]
        data = np.loadtxt(f)
        my_dict[name]['Data'] = data  # stores the form factor data
        my_dict[name]['Source'] = source  # stores the url of the source
        my_dict[name]['Updated'] = date  # stores the date last updated

    return my_dict

def createFFM(directory, filename, my_dict,source, date):
    import os
    import pickle

    f = os.path.join(directory, filename)
    if not (filename.startswith('README')):
        name = f.split('.')[0]
        data = np.loadtxt(f)
        my_dict[name]['Data'] = data  # stores the form factor data
        my_dict[name]['Source'] = source  # stores the url of the source
        my_dict[name]['Updated'] = date  # stores the date last updated

    return my_dict

if __name__ == '__main__':

    from scipy.interpolate import UnivariateSpline
    import material_structure as ms
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim7uc_unitCell_v3.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()
    Theta = 15
    energy = 640.2
    qz = np.sin(Theta * np.pi / 180) * (energy * 0.001013546143)

    """
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
    
    


    thickness, density, mag_density = sample.density_profile()

    my_keys = ['Sr', 'Ti', 'La', 'Mn','O']
    my_keys = ['Mn2','Mn3']


    d = 37.550
    idx = [i for i in range(len(thickness)) if thickness[i] < d]


    electron = {'Sr':2,'La':3,'Ti':4,'Mn2':2,'Mn3':3, 'O':-2}
    density['Mn'] = density['Mn2'] + density['Mn3']

    from scipy import integrate

    total = integrate.trapezoid(mag_density['Mn3'], x=thickness)
    print(total/10)
    #print(total/10/4)

    x = [4,7,10]
    y = [0.009153213732481595/4, 0.06715542943115324/7, 0.09526111771582044/10]
    plt.figure(1)
    plt.bar(x,y,)
    plt.xlabel('Sample (uc)')
    plt.ylabel('Magnetic Density per units cell')
    plt.show()

    plt.figure(2)
    for key in my_keys:
        plt.plot(thickness, density[key], c=np.random.rand(3,))
    
    plt.ylabel('Density (mol/cm^3)')
    plt.xlabel('Thickness (angstroms)')
    plt.legend(my_keys)
    plt.show()

    oxidation = np.zeros(len(thickness))
    total = np.zeros(len(thickness))
    for key in electron.keys():
        oxidation = oxidation + electron[key]*density[key]
        total = total+density[key]

    oxidation = oxidation/total

    plt.figure()
    plt.plot(thickness[idx], oxidation[idx])
    plt.show()
    
    x = [4,7,10]
    y = [0.0007044,0.00220607,0.00279876]

    plt.figure(3)
    plt.plot(x,y)
    plt.ylabel('Average Magnetic Density (mol/cm^2)')
    plt.xlabel('Unit Cells (uc)')


    
    
    #KK = variationKK(E_prime,E0,E2,E4)
    

    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"

    struct_names, mag_names = mm._use_given_ff("//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas")  # look for form factors in directory

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()


    
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)

    keys = ['26_452.77_S', '35_460.76_S', '19_500.71_S', '31_635.99_S', '22_640.99_S', '24_644.02_S', '33_834.59_S',
            '9_642.12_LC', '10_642.12_RC', '9-10_642.12_AC_Asymm', '13_644.03_LC', '14_644.03_RC',
            '13-14_644.03_AC_Asymm',
            '16_653.06_LC', '17_653.06_RC', '16-17_653.06_AC_Asymm']


    
    destination = '\\cabinet\work$\lsk601\My Documents\Data_for_Jesus\RXR-Twente-E1-150K-simulation-0.1eV'
    for i in range(0,300+1,1):
        E = 450 + i*0.1
        Theta = np.linspace(0.1, 60, num=2000)
        qz = np.sin(Theta * np.pi / 180) * E * (0.001013546247)

        qz, R = sample.reflectivity(E, qz)

        filename = 'E1_' + str(E)
        dat = np.transpose(np.array([qz, R['S']]))
        file = r"\\cabinet\work$\lsk601\My Documents\Data_for_Jesus\RXR-Twente-E1-150K-simulation-0.1eV"
        file = file + '\E1_' + str(E)
        np.savetxt(file, dat)



    #np.savetxt('E1_' + str(E), dat)
    

    #E = np.arange(10,3000,1)
    #imaginary = np.zeros(len(E))
    #real = np.zeros(len(E))

    #data = np.transpose(np.array([E, imaginary, real]))

    #np.savetxt('Z.txt', data,fmt='%.4s')

    
    dir= "C:/Users/lsk601/PycharmProjects/MaterialReflection/Magnetic_Scattering_Factor"
    import pickle
    import os

    #with open('form_factor.pkl', 'rb') as f:
    #    my_dict = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    #f.close()

    my_dict = dict()

    source = 'The Center for X-ray Optics'
    date = 'April 6, 2023'
    for filename in os.listdir(dir):
        f = os.path.join(dir, filename)
        if not (filename.startswith('README')):
            name = filename.split('.')[0]
            data = np.loadtxt(f)
            my_dict[name] = dict()
            my_dict[name]['Data'] = data  # stores the form factor data
            my_dict[name]['Source'] = source  # stores the url of the source
            my_dict[name]['Updated'] = date  # stores the date last updated


    import pickle
    with open('form_factor_magnetic.pkl', 'wb') as handle:
        pickle.dump(my_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    

    """
    import pickle
    with open('form_factor.pkl', 'rb') as f:
        ff = pickle.load(f)  # This is made a global variable so we do not have to keep on loading in the file
    f.close()

    print(ff['Mn'])