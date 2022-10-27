import matplotlib.pyplot as plt
from prettytable import PrettyTable
import os
from material_structure import *
from material_model import *
from time import *
import ast
import h5py




def WriteDataHDF5(fname, AScans,AInfo, EScans, EInfo, sample):
    """
    Purpose: Take in data and sample model and save it as an hdf5 file
    :param fname: File name with .hdf5 file extension
    :param AScans: Reflectivity scan data from ProcessRXR
    :param AInfo: Reflectivity scan information from ProcessRXR
    :param EScans: Energy scan data from Process RXR
    :param EInfo: Energy scan information from Process RXR
    :param sample: sample object
    :return:
    """

    # Checking if fname already exists
    cwd = os.getcwd()
    path = cwd + '/' + fname
    if os.path.exists(path):
        raise OSError("HDF5 file already exists. To write new file remove the old file from the current working directory.")

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # creating group that will contain the sample information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)


    scattering_factor = []
    mag_scattering_factor = []

    dsLayer = 0
    for my_layer in sample.structure:

        # Layer information
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            # Element information
            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # May be changed in the future as layermagnetized also contains this information
            # Original implemented to avoid problem of trying to load in the magnetic data that does not exist
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

                if type(my_layer[ele].mag_scattering_factor) is list or type(my_layer[ele].mag_scattering_factor) is np.ndarray:
                    for msf in my_layer[ele].mag_scattering_factor:
                        if msf not in mag_scattering_factor:
                            mag_scattering_factor.append([msf, 0])
                else:
                    if my_layer[ele].mag_scattering_factor not in mag_scattering_factor:
                        mag_scattering_factor.append([my_layer[ele].mag_scattering_factor, 0])


            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            if type(my_layer[ele].scattering_factor) is list or type(my_layer[ele].scattering_factor) is np.ndarray:
                for my_sf in my_layer[ele].scattering_factor:
                    if my_sf not in scattering_factor:
                        scattering_factor.append([my_sf, 0])
            else:
                if my_layer[ele].scattering_factor not in scattering_factor:
                    scattering_factor.append([my_layer[ele].scattering_factor, 0])
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    # setting the eShift parameters depending on the user input
    if len(list(sample.eShift.keys())) == 0:
        grp1.attrs['FormFactors'] = str(scattering_factor)  # scattering factors
    else:
        temp = []
        for key in list(sample.eShift.key()):
            temp.append([key, sample.eShift[key]])
        grp1.attrs['FormFactors'] = str(temp)

    # setting the mag eShift depending on the user input
    if len(list(sample.mag_eShift.keys())) == 0:
        grp1.attrs['MagFormFactors'] = str(mag_scattering_factor)  # magnetic scattering factors
    else:
        temp = []
        for key in list(sample.mag_eShift.key()):
            temp.append([key, sample.mag_eShift[key]])
        grp1.attrs['FormFactors'] = str(temp)

    # Loading in the experimental data and simulated data
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    #Loading reflectivity scans
    dsNum = 1
    for i in range(len(AScans)):
        qz = AScans[i][:,2]  # momentum transfer
        R0 = AScans[i][:,3]  # reflectivity
        energy = float(AInfo[i][3])  # energy of scan
        Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / pi  # initial angle
        dat = np.array([qz, Theta, R0])

        energy = float(AInfo[i][3])
        polarization = 'S'
        if (AInfo[i][1] == "S"):
            polarization = 'S'
        elif (AInfo[i][1] == "P"):
            polarization = 'P'
        elif (AInfo[i][1] == "L"):
            polarization = 'LC'
        elif (AInfo[i][1] == "R"):
            polarization = 'RC'
        datasetpoints = len(qz)

        name = str(AInfo[i][0]) + "_" + str(np.round(energy,2)) + "_" + polarization
        qz, R = sample.reflectivity(energy, qz)
        sim = np.array([qz, Theta, R[polarization]])
        dset = grpR.create_dataset(name, data=dat)
        dset1 = subR.create_dataset(name, data=sim)

        dset.attrs['Energy'] = float(energy)
        dset1.attrs['Energy'] = float(energy)

        # newly added
        dset.attrs['Background Shift'] = float(0)
        dset1.attrs['Background Shift'] = float(0)

        # newly added
        dset.attrs['Scaling Factor'] = float(1)
        dset1.attrs['Scaling Factor'] = float(1)

        dset.attrs['Polarization'] = str(polarization)
        dset1.attrs['Polarization'] = str(polarization)

        dset.attrs['DataPoints'] = int(datasetpoints)
        dset1.attrs['DataPoints'] = int(datasetpoints)

        dset.attrs['DatasetNumber'] = int(dsNum)
        dset1.attrs['DatasetNumber'] = int(dsNum)

        dsNum = dsNum + 1
        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):

                qz = AScans[i][:,2]
                A = (AScans[i - 1][:,3] - AScans[i][:,3]) / (AScans[i - 1][:,3]) + AScans[i][:,3]
                energy = float(AInfo[i][3])
                Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / pi  # initial angle
                dat = np.array([qz, Theta, A])

                datasetpoints = len(qz)

                if (AInfo[i - 1][1] == "S" or AInfo[i - 1][1] == "P"):
                    polarization = "AL"
                elif (AInfo[i - 1][1] == "L" or AInfo[i - 1][1] == "R"):
                    polarization = "AC"
                name = str(AInfo[i-1][0])+ "-" + str(AInfo[i][0]) + "_" + str(np.round(energy,2)) +"_"+ polarization + "_Asymm"
                dset = grpR.create_dataset(name, data=dat)
                qz, R = sample.reflectivity(energy, qz)
                sim = np.array([qz,Theta,R[polarization]])
                dset1 = subR.create_dataset(name,data = sim)

                dset.attrs['Energy'] = float(energy)
                dset1.attrs['Energy'] = float(energy)

                # newly added
                dset.attrs['Background Shift'] = float(0)
                dset1.attrs['Background Shift'] = float(0)

                # newly added
                dset.attrs['Scaling Factor'] = float(1)
                dset1.attrs['Scaling Factor'] = float(1)

                dset.attrs['Polarization'] = polarization
                dset1.attrs['Polarization'] = polarization

                dset.attrs['DataPoints'] = int(datasetpoints)
                dset1.attrs['DataPoints'] = int(datasetpoints)

                dset.attrs['DatasetNumber'] = int(dsNum)
                dset1.attrs['DatasetNumber'] = int(dsNum)
                dsNum = dsNum + 1


    for i in range(len(EScans)):

        qz = EScans[i][:,2]
        R0 = EScans[i][:,3]
        E = EScans[i][:,0]

        Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / pi  # initial angle
        datasetpoints = len(qz)
        dat = np.array([qz,Theta,R0,E])

        energy = float(EInfo[i][3])
        polarization = 'S'
        if (EInfo[i][1] == "S"):
            polarization = 'S'
        elif (EInfo[i][1] == "P"):
            polarization = 'P'
        elif (EInfo[i][1] == "L"):
            polarization = 'LC'
        elif (EInfo[i][1] == "R"):
            polarization = 'RC'
        angle = float(EInfo[i][4])

        name = str(EInfo[i][0]) + "_E" + str(np.round(energy,2)) + "_Th" + str(np.round(angle,2)) + "_" + polarization
        E, R = sample.energy_scan(angle,E)
        sim = np.array([qz, Theta, R[polarization], E])

        dset = grpE.create_dataset(name, data=dat)
        dset1 = subE.create_dataset(name, data=sim)

        dset.attrs['Energy'] = float(energy)
        dset1.attrs['Energy'] = float(energy)

        # newly added
        dset.attrs['Background Shift'] = float(0)
        dset1.attrs['Background Shift'] = float(0)

        # newly added
        dset.attrs['Scaling Factor'] = float(1)
        dset1.attrs['Scaling Factor'] = float(1)

        dset.attrs['Polarization'] = polarization
        dset1.attrs['Polarization'] = polarization

        dset.attrs['DataPoints'] = int(datasetpoints)
        dset1.attrs['DataPoints'] = int(datasetpoints)

        dset.attrs['Angle'] = float(angle)
        dset1.attrs['Angle'] = float(angle)

        dset.attrs['DatasetNumber'] = int(dsNum)
        dset1.attrs['DatasetNumber'] = int(dsNum)
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
                qz = EScans[i][:,2]
                E = EScans[i][:,0]
                A = (EScans[i-1][:,3]-EScans[i][:,3])/(EScans[i-1][:,3]+EScans[i][:,3])
                Theta = np.arcsin(qz / E / (0.001013546247)) * 180 / pi  # initial angle
                dat = np.array([qz, Theta, A,E])
                energy = float(EInfo[i][3])

                if (EInfo[i - 1][1] == "S" or EInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                elif (EInfo[i - 1][1] == "L" or EInfo[i - 1][1] == "R"):
                    polarization = 'AC'

                angle = float(EInfo[i][4])
                name = str(EInfo[i-1][0])+"-"+ str(EInfo[i][0]) + "_E" + str(np.round(energy,2)) + "_Th" + str(np.round(angle, 2)) + "_" + polarization + "_Asymm"
                E, R = sample.energy_scan(angle, E)
                sim = np.array([qz, Theta, R[polarization], E])

                dset = grpE.create_dataset(name, data=dat)
                dset1 = subE.create_dataset(name, data=sim)

                dset.attrs['Energy'] = float(energy)
                dset1.attrs['Energy'] = float(energy)

                # newly added
                dset.attrs['Background Shift'] = float(0)
                dset1.attrs['Background Shift'] = float(0)

                # newly added
                dset.attrs['Scaling Factor'] = float(1)
                dset1.attrs['Scaling Factor'] = float(1)

                dset.attrs['Polarization'] = polarization
                dset1.attrs['Polarization'] = polarization

                dset.attrs['Angle'] = float(angle)
                dset1.attrs['Angle'] = float(angle)

                dset.attrs['DataPoints'] = int(datasetpoints)
                dset1.attrs['DataPoints'] = int(datasetpoints)

                dset.attrs['DatasetNumber'] = int(dsNum)
                dset1.attrs['DatasetNumber'] = int(dsNum)

                dsNum = dsNum + 1

    f.close()

def ReadSampleHDF5(fname):

    """
    Purpose: Read the sample info from hdf5 file and recreate the sample object
    :param fname: File name
    :return: sample - the recreated sample object
    """
    f = h5py.File(fname, 'r')

    S = f['Sample']

    # Retieves the general information of the sample
    m = int(S.attrs['NumberLayers'])
    sample = slab(m)
    sample.poly_elements = ast.literal_eval(S.attrs['PolyElements'])
    sample.mag_elements = ast.literal_eval(S.attrs['MagElements'])
    sample.layer_magnetized = S.attrs['LayerMagnetized']
    scattering_factor = ast.literal_eval(S.attrs['FormFactors'])
    mag_scattering_factors = ast.literal_eval(S.attrs['MagFormFactors'])

    for ff in scattering_factor:
        sample.eShift[ff[0]] = ff[1]

    for ffm in mag_scattering_factors:
        sample.mag_eShift[ffm[0]] = ff[1]

    # Retrieves the general layer information
    for lay_key in S.keys():
        layer = S[lay_key]
        formula = layer.attrs['Formula']
        lay_num = int(layer.attrs['LayerNumber'])
        sample.addlayer(lay_num, formula,20, density=1)  #pre-initialize parameters to random numbers for each layer

        # retrieves the information for each element
        for ele_key in layer.keys():
            element = layer[ele_key]
            sample.structure[lay_num][ele_key].molar_mass = element.attrs['MolarMass']
            sample.structure[lay_num][ele_key].density = element.attrs['Density']
            sample.structure[lay_num][ele_key].thickness = element.attrs['Thickness']
            sample.structure[lay_num][ele_key].roughness = element.attrs['Roughness']
            sample.structure[lay_num][ele_key].linked_roughness = element.attrs['LinkedRoughness']
            sample.structure[lay_num][ele_key].poly_ratio = element.attrs['PolyRatio']
            sample.structure[lay_num][ele_key].polymorph = element.attrs['Polymorph']
            sample.structure[lay_num][ele_key].gamma = element.attrs['Gamma']
            sample.structure[lay_num][ele_key].phi = element.attrs['Phi']
            if element.attrs['Magnetic']:
                sample.structure[lay_num][ele_key].mag_density = element.attrs['MagDensity']
                sample.structure[lay_num][ele_key].mag_scattering_factor = element.attrs['MagScatteringFactor']

            sample.structure[lay_num][ele_key].scattering_factor = element.attrs['ScatteringFactor']
            sample.structure[lay_num][ele_key].position = element.attrs['Position']

    f.close()
    return sample

def ReadDataHDF5(fname):
    """
    Purpose: Reads in the experimental and simulated data from hdf5 file and then plots spectrum chosen by user
    :param fname: File name
    :return:
    """

    f = h5py.File(fname, 'r')
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
        mydata = RS[Rkey]
        data_dict[Rkey] = RS[Rkey]
        sim_dict[Rkey] = SimR[Rkey]
        Rdat = [int(mydata.attrs['DatasetNumber']),'Reflectivity', Rkey]
        data.append(Rdat)

    for Ekey in ES.keys():
        mydata = ES[Ekey]
        data_dict[Ekey] = ES[Ekey]
        sim_dict[Ekey] = SimE[Ekey]
        Edat = [int(mydata.attrs['DatasetNumber']),'Energy', Ekey]
        data.append(Edat)

    # Sorts data in appropriate order
    data = np.array(data)
    sort_idx = np.argsort(data[:,0].astype(int))
    data = data[sort_idx]



    data_dict = hdf5ToDict(data_dict)
    sim_dict = hdf5ToDict(sim_dict)
    f.close()

    return data, data_dict, sim_dict


def WriteSampleHDF5(fname, sample):
    """
    Purpose: Write a new sample to the hdf5 file fname
    :param fname: File name
    :param sample: new sample information as sample type
    :return:
    """

    # deletes the previous sample information
    with h5py.File(fname, "a") as f:
        if 'Sample' in f:
            del f['Sample']
    f.close()

    f = h5py.File(fname, "a")
    grp1 = f.create_group('Sample')

    # Sets the general sample information
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)

    scattering_factor = []
    mag_scattering_factor = sample.mag_eshift

    for key in list(sample.eShift.keys()):
        scattering_factor.append([key,sample.eShift[key]])

    for key in list(sample.mag_eShift.keys()):
        mag_scattering_factor.append([key,sample.mag_eShift[key]])

    grp1.attrs['FormFactors'] = np.array(scattering_factor)
    grp1.attrs['MagFormFactors'] = np.array(mag_scattering_factor)

    # Sets the information for each layer
    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)
        ### Need to change back to previous version

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        # Sets the information for each element
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1



def WriteHDF5Simulation():
    """
        Purpose: Write a new sample to the hdf5 file fname
        :param fname: File name
        :param sample: new sample information as sample type
        :return:
        """

    # deletes the previous sample information
    with h5py.File(fname, "a") as f:
        if 'Sample' in f:
            del f['Sample']
    f.close()

    f = h5py.File(fname, "a")
    grp1 = f.create_group('Sample')

    # Sets the general sample information
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)

    scattering_factor = []
    mag_scattering_factor = sample.mag_eshift

    for key in list(sample.eShift.keys()):
        scattering_factor.append([key, sample.eShift[key]])

    for key in list(sample.mag_eShift.keys()):
        mag_scattering_factor.append([key, sample.mag_eShift[key]])

    grp1.attrs['FormFactors'] = np.array(scattering_factor)
    grp1.attrs['MagFormFactors'] = np.array(mag_scattering_factor)

    # Sets the information for each layer
    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)
        ### Need to change back to previous version

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        # Sets the information for each element
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f["Simulated_data"]

    grpR = grp2["Reflectivity_Scan"]
    grpE = grp2["Energy_Scan"]

    # Recalculate the simulated reflectivity scan data
    for key in list(grpR.keys()):
        dset = grpR[key]  # retrieves the old data set
        Rdata = list(dset)  # gets the dataset information
        qz = np.array(Rdata[0])  # retrieves the momentum transfer
        E = float(dset.attrs['Energy'])  # retrieves the energy
        polarization = dset.attrs['Polarization']  # retrieves the polarization

        qz, R = sample.reflectivity(E, qz)  # computes reflecticity
        Rdata[2] = R[polarization]  # retrieves appropriate reflectivity for the correct polarization
        dset[...] = Rdata  # overwrites previous dataset with new data

    # Recalculates the simulated energy scan data
    for key in list(grpE.keys()):
        dset = grpE[key]  # retrieves the old data set
        Edata = list(dset)  # gets the dataset information
        E = Edata[3]  # retrieves the numpy energy array
        theta = dset.attrs['Angle']  # retrieves the angle
        polarization = dset.attrs['Polarization']  # retrieves the polarization
        E, R = sample.energy_scan(theta, E)  # recomputes the energy scan
        Edata[2] = R[polarization]  # retrieves appropriate reflectivity with the correct polarization
        dset[...] = Edata  # overwrites previous data with new data

    f.close()

def updateHDF5Data(fname, scalingFactor, backgroundShift):

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    f = h5py.File(fname, 'r')
    experiment = f['Experimental_data']

    RS = experiment['Reflectivity_Scan']

    ES = experiment['Energy_Scan']


    # Recalculate the simulated reflectivity scan data
    idx = 0
    for key in list(RS.keys()):
        RS[key].attrs['ScalingFactor'] = scalingFactor[idx]
        RS[key].attrs['BackgroundShift'] = backgroundShift[idx]
        idx = idx + 1

    # Recalculates the simulated energy scan data
    for key in list(ES.keys()):
        ES[key].attrs['ScalingFactor'] = scalingFactor[idx]
        ES[key].attrs['BackgroundShift'] = backgroundShift[idx]
        idx = idx + 1

    f.close()


def WriteSampleASCII(file,sample):
    """
    Purpose: Write the sample as an ASCII file
    :param file: File name
    :param sample: sample information as a slab object
    :return:
    """

    file.write("# Structure \n")  # header defining that the sample information is starting
    n = len(sample.structure)

    # General information for the sample model
    file.write("numberlayers = %s \n" % str(n))
    file.write("polyelements = %s \n" % str(sample.poly_elements))
    file.write("magelements = %s \n" % str(sample.mag_elements))
    file.write("layermagnetized = %s \n" % sample.layer_magnetized)
    file.write("scalingfactor = %s \n" % str(sample.scaling_factor))
    file.write("backgroundshift = %s \n\n" % str(sample.background_shift))

    # writing the layer and element information
    num_lay = 0
    for layer in sample.structure:

        # General layer information
        file.write("layer = %s \n" % str(num_lay))

        # Reconstructing the chemical formula
        formula = ''
        for ele in list(layer.keys()):
            stoich = layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        file.write("formula = %s \n\n" % formula)

        # writing the element information
        for ele in layer.keys():
            file.write("element = %s \n" % ele)
            file.write("molarmass = %f \n" % layer[ele].molar_mass)
            file.write("density = %f \n" % layer[ele].density)
            file.write("thickness = %f \n" % layer[ele].thickness)
            file.write("roughness = %f \n" % layer[ele].roughness)
            file.write("linkedroughness = %f \n" % layer[ele].linked_roughness)
            file.write("scatteringfactor = %s \n" % layer[ele].scattering_factor)
            file.write("polymorph = %s \n" % layer[ele].polymorph)

            poly_ratio = layer[ele].poly_ratio
            if type(poly_ratio) != int:
                poly_ratio = [str(poly) for poly in poly_ratio]
            file.write("polyratio = %s \n" % poly_ratio)


            file.write("gamma = %f \n" % layer[ele].gamma)
            file.write("phi = %f \n" % layer[ele].phi)

            mag_density = layer[ele].mag_density
            mag_density = [str(x) for x in mag_density]
            file.write("magdensity = %s \n" % mag_density)
            sfm = layer[ele].mag_scattering_factor
            file.write("magscatteringfactor = %s \n" % sfm)
            file.write("position = %s \n" % layer[ele].position)
            file.write("\n")

        num_lay = num_lay + 1

def WriteExperimentalDataASCII(file, AScans,AInfo,EScans,EInfo):
    """
    Purpose: Write experimental data to .all file for program use
    :param fname: Name of data file
    :param AScans: Reflection scan experimental data
    :param AInfo: Info about reflection data scans
    :param EScans: Energy scan experimental data
    :param EInfo: Energy scan info
    :param header: Additional header
    :return:
    """


    file.write("# Experimental_Data \n")  # header that states writing experimental data

    # Writing the reflectivity scan data
    dsNum = 1
    for i in range(len(AScans)):

        # Writing the general information of each scan
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % AInfo[i][3])
        if (AInfo[i][1] == "S"):
            file.write("polarization = S \n")
        elif (AInfo[i][1] == "P"):
            file.write("polarization = P \n")
        elif (AInfo[i][1] == "L"):
            file.write("polarization = LC \n")
        elif (AInfo[i][1] == "R"):
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

        # writing the data of each scan
        for j in range(len(AScans[i][:, 0])):
            file.write("dataset_qz = %f \n" % AScans[i][j][2])  # momentum transfer
            file.write("dataset_R0 = %e \n" % AScans[i][j][3])  # reflectivity
            # file.write("dataset_eng = %f \n" % AScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):
                # General info of each scan
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(AInfo[i - 1][0]) + "-" + str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i - 1][
                    1] + "-" + AInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % AInfo[i][3])
                if (AInfo[i-1][1] == "S" or AInfo[i-1][1] =="P"):
                    file.write("polarization = AL \n")
                elif (AInfo[i-1][1] == "L" or AInfo[i-1][1] == "R"):
                    file.write("polarization = AC \n")
                file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

                # writing data
                for j in range(len(AScans[i][:, 0])):
                    file.write("dataset_qz = %f \n" % AScans[i][j][2])  # momentum transfer
                    file.write("dataset_A = %e \n" % (
                                (AScans[i - 1][j][3] - AScans[i][j][3]) / (AScans[i - 1][j][3] + AScans[i][j][3])))

                file.write("\n\n")
                dsNum = dsNum + 1

    # Writing energy scan data
    for i in range(len(EScans)):
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]), 2)) + "_Th" + str(round(float(EInfo[i][4]), 2)) + "_" + EInfo[i - 1][
                   1] + "-" + EInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % EInfo[i][3])
        if (EInfo[i][1] == "S"):
            file.write("polarization = S \n")
        elif (EInfo[i][1] == "P"):
            file.write("polarization = P \n")
        elif (EInfo[i][1] == "L"):
            file.write("polarization = LC \n")
        elif (EInfo[i][1] == "R"):
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(EScans[i][:, 0]))

        # Writing the data for each energy scan
        for j in range(len(EScans[i][:, 0])):
            file.write("dataset_qz = %f \n" % EScans[i][j][2])
            file.write("dataset_R0 = %e \n" % EScans[i][j][3])
            file.write("dataset_eng = %f \n" % EScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):

                # General information
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(EInfo[i - 1][0]) + "-" + str(EInfo[i][0]) + "_E" + str(
                    round(float(EInfo[i][3]), 2)) + "_Th" + str(round(float(EInfo[i][4]), 2)) + "_" + EInfo[i - 1][
                           1] + "-" + EInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % EInfo[i][3])
                if (EInfo[i-1][1] == "S" or EInfo[i-1][1] =="P"):
                    file.write("polarization = AL \n")
                elif (EInfo[i-1][1] == "L" or EInfo[i-1][1] == "R"):
                    file.write("polarization = AC \n")

                # Writing the data for each energy scan
                for j in range(len(EScans[i][:, 0])):
                    file.write("dataset_qz = %f \n" % EScans[i][j][2])  # momentum transfer
                    file.write("dataset_A = %e \n" % (
                                (EScans[i - 1][j][3] - EScans[i][j][3]) / (EScans[i - 1][j][3] + EScans[i][j][3])))
                    file.write("dataset_eng = %f \n" % EScans[i][j][0])

                file.write("\n\n")
                dsNum = dsNum + 1



def WriteSimulationASCII(file, AScans,AInfo,EScans,EInfo, sample):
    """
        Purpose: Write experimental data to .all file for program use
        :param fname: Name of data file
        :param AScans: Reflection scan experimental data
        :param AInfo: Info about reflection data scans
        :param EScans: Energy scan experimental data
        :param EInfo: Energy scan info
        :param header: Additional header
        :return:
        """

    file.write("# Simulation \n")  # header that states start of simulation data
    qz = list()

    # Writing reflectivity data
    dsNum = 1
    polarization = 'S'
    for i in range(len(AScans)):
        # General information
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % AInfo[i][3])
        if (AInfo[i][1] == "S"):
            polarization = 'S'
            file.write("polarization = S \n")
        elif (AInfo[i][1] == "P"):
            polarization = 'P'
            file.write("polarization = P \n")
        elif (AInfo[i][1] == "L"):
            polarization = 'LC'
            file.write("polarization = LC \n")
        elif (AInfo[i][1] == "R"):
            polarization = 'RC'
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

        qz = AScans[i][:,2]  # retrieves momentum transfer
        qz, R = sample.reflectivity(float(AInfo[i][3]), qz)  # computes the reflectivity
        R = R[polarization]  # retrieves the reflectivity of the correct polarization

        # Writes the simulated data
        for j in range(len(qz)):
            file.write("dataset_qz = %f \n" % qz[j])  # momentum transfer
            file.write("dataset_R0 = %e \n" % R[j])  # reflectivity

        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):
                # General information
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(AInfo[i - 1][0]) + "-" + str(AInfo[i][0]) + "_A_" + AInfo[i][3] + "_" + AInfo[i - 1][
                    1] + "-" + AInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % AInfo[i][3])
                if (AInfo[i - 1][1] == "S" or AInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                    file.write("polarization = AL \n")
                elif (AInfo[i - 1][1] == "L" or AInfo[i - 1][1] == "R"):
                    polarization = 'AC'
                    file.write("polarization = AC \n")
                file.write("datasetpoints = %d \n" % len(AScans[i][:, 0]))

                qz = AScans[i][:,2]  # retrieves the momentum transfer
                qz, R = sample.reflectivity(float(AInfo[i][3]), qz)  # computes reflectivity
                R = R[polarization]  # retrieves the reflectivity of the correct polarization

                # Writes the simulated data
                for j in range(len(qz)):
                    file.write("dataset_qz = %f \n" % qz[j])
                    # print(AScans[i-1][j][3]+AScans[i][j][3])
                    file.write("dataset_A = %e \n" % R[j])
                    # file.write("dataset_eng = %f \n" % AScans[i][j][0])
                file.write("\n\n")
                dsNum = dsNum + 1

    # Writing the energy scan simulation data
    for i in range(len(EScans)):
        # General information
        file.write("datasetnumber = %d \n" % dsNum)
        name = str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]), 2)) + "_Th" + str(
            round(float(EInfo[i][4]), 2)) + "_" + EInfo[i][1]
        file.write("datasettitle = %s \n" % name)
        file.write("datasetenergy = %s \n" % EInfo[i][3])
        if (EInfo[i][1] == "S"):
            polarization = 'S'
            file.write("polarization = S \n")
        elif (EInfo[i][1] == "P"):
            polarization = 'P'
            file.write("polarization = P \n")
        elif (EInfo[i][1] == "L"):
            polarization = 'S'
            file.write("polarization = LC \n")
        elif (EInfo[i][1] == "R"):
            polarization = 'S'
            file.write("polarization = RC \n")
        file.write("datasetpoints = %d \n" % len(EScans[i][:, 0]))

        E = EScans[i][:,0]  # numpy energy array
        E, R = sample.energy_scan(float(EInfo[i][4]), E)  # computes the energy scan
        R = R[polarization]  # retrieves the energy scan for the correct polarization

        # Writes the energy scan data
        for j in range(len(E)):
            file.write("dataset_R0 = %e \n" % R[j])  # reflectivity
            file.write("dataset_eng = %f \n" % E[j])  # energy
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
                #General information
                file.write("datasetnumber = %d \n" % dsNum)
                name = str(EInfo[i - 1][0]) + "-" + str(EInfo[i][0]) + "_E" + str(
                    round(float(EInfo[i][3]), 2)) + "_Th" + str(round(float(EInfo[i][4]), 2)) + "_" + EInfo[i - 1][
                           1] + "-" + EInfo[i][1] + "_Asymm"
                file.write("datasettitle = %s \n" % name)
                file.write("datasetenergy = %s \n" % EInfo[i][3])
                if (EInfo[i - 1][1] == "S" or EInfo[i - 1][1] == "P"):
                    polarization = 'AL'
                    file.write("polarization = AL \n")
                elif (EInfo[i - 1][1] == "L" or EInfo[i - 1][1] == "R"):
                    polarization = 'AC'
                    file.write("polarization = AC \n")

                E = EScans[i][:,0]  # numpy energy array
                E, R = sample.energy_scan(float(EInfo[i][4]), E)  # compute energy scan
                R = R[polarization]  # retrieve energy scan with correct polarization

                # write energy scan data
                for j in range(len(E)):
                    file.write("dataset_A = %e \n" % R[j])  # reflectivity
                    file.write("dataset_eng = %f \n" % E[j])  # energy

                file.write("\n\n")
                dsNum = dsNum + 1


def WriteDataASCII(fname,AScans,AInfo,EScans,EInfo, sample):
    """
    Purpose: Write experimental, simulated data, and sample info to an ASCII file with .all extension
    :param fname: File name
    :param AScans: Array that contains the reflectivity scan experimental data
    :param AInfo: Array that contains the information related to each reflecitivty scan
    :param EScans: Array that contains the energy scan experimental data
    :param EInfo: Array that contains the information related to each energy scan
    :param sample: Slab object containing sample information
    :return:
    """
    # If you decide to use your own functions to create an ASCII file you must include the following headers before
    # you write the data
    #  - # Structure
    #  - # Experimental_Data
    #  - # Simulation

    # Information for the data fitting and neural network parameters will be added later on

    file = open(fname, "w")  # creating a new file with name fname
    WriteSampleASCII(file, sample)  # writes the sample information
    WriteExperimentalDataASCII(file, AScans,AInfo,EScans,EInfo)  # writes the experimental data
    WriteSimulationASCII(file, AScans, AInfo, EScans, EInfo, sample)  # writes the simulation data
    f.close()

def getScanInfo(title):
    """
    Purpose: Retrieves important information in the scan title
    :param title: title/label of the scan
    :return:
    """
    title = title.split('_')  # split title into an array

    angle = None
    scan_number = title[0]

    # Determines the scan type
    if title[1] == 'A':
        scanType = 'Reflectivity'
    else:
        scanType = 'Energy'
        angle = title[2].replace('Th','')

    return scan_number, scanType, angle

def createNewDict():
    """
    Purpose: Initializes a dicitonary with the required keys and values
    :return:
    """
    my_dict = dict()
    my_dict['scanNumber'] = None
    my_dict['dataNumber'] = None
    my_dict['scanType'] = None
    my_dict['angle'] = None
    my_dict['energy'] = None
    my_dict['polarization'] = None
    my_dict['numberPoints'] = None
    return my_dict

def ConvertASCIItoHDF5(fascii, fhdf5):
    """
    Purpose: Converts and ASCII file to an hdf5 file with the proper format
    :param fname: File name of the ASCII file
    :return:
    Current implementation does not allow user to change the hdf5 file name from the ASCII file name
    """

    # checks to make sure file type of ASCII and HDF5 files are correct
    if not(fascii.endswith('.all')):
        raise NameError('File name must be a .all file type')
    if not(fhdf5.endswith('.h5')):
        raise NameError('File name must be a .all file type')

    Sinfo, Sscan, SimInfo, SimScan, sample = ReadDataASCII(fascii)  # retrieves information from ASCII file

    # checking to make sure that the HDF5 file asked for does not already exist
    cwd = os.getcwd()
    path = cwd + '/' + fhdf5
    if os.path.exists(path):
        raise OSError(
            "HDF5 file already exists. To write new HDF5 file remove the old file from the current working directory.")

    f = h5py.File(fhdf5, 'a')  # creates a new HDF5 file

    # Sample conversion of general information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)
    grp1.attrs['ScalingFactor'] = float(sample.scaling_factor)
    grp1.attrs['BackgroundShift'] = float(sample.background_shift)


    # general layer information
    dsLayer = 0
    for my_layer in sample.structure:

        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        # reconstructing the chemical formula input by user
        formula = ''
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)
            stoich = int(my_layer[ele].stoichiometry)
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # Magnetic element check
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        layer.attrs['Formula'] = formula
        dsLayer = dsLayer + 1

    # Experimental and simulated data conversion

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    name = ''

    dsNum = 1
    for i in range(len(Sinfo)):
        info = Sinfo[i]
        data = Sscan[i]
        simData = SimScan[i]

        scanType = info['scanType']

        # Reflectivity scans
        if scanType == "Reflectivity":
            qz = np.array(data[0])
            R0 = np.array(data[1])
            R0sim = np.array(simData[1])
            energy = info['energy']
            Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / pi  # initial angle
            dat = [qz, Theta, R0]
            sim = [qz, Theta, R0sim]

            polarization = info['polarization']

            if polarization == 'S' or polarization == 'P' or polarization == 'LC' or polarization == 'RC':
                name = str(info['dataNumber']) + "_" + str(np.round(energy,2)) + "_" + info['polarization']
            elif polarization == 'AL' or polarization == 'AC':
                name = str(info['dataNumber']) + "_" +str(np.round(energy,2)) + "_" + info['polarization'] +"_Asymm"

            datasetpoints = len(qz)

            dset = grpR.create_dataset(name, data=dat)
            dset1 = subR.create_dataset(name, data=sim)


            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = polarization
            dset1.attrs['Polarization'] = polarization

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)

            dsNum = dsNum + 1

        elif scanType == 'Energy':
            E = np.array(data[0])
            R = np.array(data[1])
            Rsim = np.array(simData[1])
            qz = np.array(data[2])
            energy = info['energy']
            angle = info['angle']
            Theta = np.arcsin(qz / energy / (0.001013546247)) * 180 / pi  # initial angle

            dat = [qz, Theta, R, E]
            sim = [qz, Theta, Rsim, E]

            polarization = info['polarization']

            if polarization == 'S' or polarization == 'P' or polarization == 'LC' or polarization == 'RC':
                name = str(info['dataNumber']) + "_E" +str(np.round(energy,2)) + "_Th" + str(np.round(angle,2)) + "_" + info['polarization']
            elif polarization == 'AL' or polarization == 'AC':
                name = str(info['dataNumber']) + "_E" +str(np.round(energy,2)) + "_Th" + str(np.round(angle,2)) + "_" + info['polarization'] + "_Asymm"

            datasetpoints = len(E)
            dset = grpE.create_dataset(name, data=dat)
            dset1 = subE.create_dataset(name, data=sim)

            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = polarization
            dset1.attrs['Polarization'] = polarization

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['Angle'] = float(angle)
            dset1.attrs['Angle'] = float(angle)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)
            dsNum = dsNum + 1


    f.close()


def ReadDataASCII(fname):
    """
    Purpose: Read .all file that contains sample and experimental information
    :param fname: Name of file contatining related information
    :return: Sscan - list of energy/reflectivity data
             Sinfo - contains a dictionary of relevant scan information
             sample - pre-initialized slab class
    """

    file = open(fname)
    Structure = False
    experimental_data = False
    simulation = False
    for line in file:
        line = line.split()
        if len(line) != 0:
            if line[0] == '#':
                if line[1] == 'Structure':
                    Structure = True
                    experimental_data = False
                    simulation = False
                    layer = 0

                    new_element = False
                    element = ''

                elif line[1] == 'Experimental_Data':

                    if simulation:
                        SimInfo[idx]['scanNumber'] = scan_number
                        SimInfo[idx]['scanType'] = scanType
                        SimInfo[idx]['angle'] = angle
                        SimInfo[idx]['energy'] = energy
                        SimInfo[idx]['polarization'] = polarization
                        SimInfo[idx]['numberPoints'] = numberPoints
                        if scanType == 'Energy':
                            SimScan.append([x_axis, y_axis, energy_qz])
                        elif scanType == 'Reflectivity':
                            SimScan.append([x_axis, y_axis])

                    Structure = False
                    experimental_data = True
                    simulation = False

                    idx = 0  # index of scan
                    Sinfo = []  # will contain information of sample
                    Sinfo.append(createNewDict())
                    Sscan = []

                    # initialization of parameters
                    x_axis = list()
                    y_axis = list()
                    energy_qz = list()
                    scan_number = 0
                    scanType = 0
                    angle = 0
                    energy = 0
                    polarization = 0
                    numberPoints = 0

                    #NewScan = False  # determines if we have a new scan
                    first = True
                    # Read in each line one at a time
                elif line[1] == 'Simulation':
                    # Attaches last experiment element
                    if experimental_data:
                        Sinfo[idx]['scanNumber'] = scan_number
                        Sinfo[idx]['scanType'] = scanType
                        Sinfo[idx]['angle'] = angle
                        Sinfo[idx]['energy'] = energy
                        Sinfo[idx]['polarization'] = polarization
                        Sinfo[idx]['numberPoints'] = numberPoints
                        if scanType == 'Energy':
                            Sscan.append([x_axis, y_axis, energy_qz])
                        elif scanType == 'Reflectivity':
                            Sscan.append([x_axis, y_axis])

                    Structure = False
                    experimental_data = False
                    simulation = True

                    idx = 0  # index of scan
                    SimInfo = []  # will contain information of sample
                    SimInfo.append(createNewDict())
                    SimScan = []

                    # initialization of parameters
                    x_axis = list()
                    y_axis = list()
                    energy_qz = list()
                    scan_number = 0
                    scanType = 0
                    angle = 0
                    energy = 0
                    polarization = 0
                    numberPoints = 0

                    first = True
                    # Read in each line one at a time
            else:
                if Structure:
                    if '=' not in line:
                        raise SyntaxError('Data file is improperly initialized.')
                    line.remove('=')  # removes the equal sign

                    # initializing the slab with correct number of layers
                    if line[0] == 'numberlayers':
                        sample = slab(int(line[1]))
                    elif line[0] == 'polyelements':
                        line.pop(0)
                        line = ''.join(line)
                        polyelements = ast.literal_eval(line)
                        sample.poly_elements = polyelements
                    elif line[0] == 'magelements':
                        line.pop(0)
                        line = ''.join(line)
                        magelements = ast.literal_eval(line)
                        sample.mag_elements = magelements

                    elif line[0] == 'layermagnetized':
                        line.pop(0)
                        line = ''.join(line)
                        layermagnetized = ast.literal_eval(line)
                        sample.layer_magnetized = layermagnetized
                    elif line[0] == 'scalingfactor':
                        scaling_factor = float(line[1])
                        sample.scaling_factor = scaling_factor
                    elif line[0] == 'backgroundshift':
                        background_shift = float(line[1])
                        sample.background_shift = background_shift
                    elif line[0] == 'layer':
                        layer = int(line[1])
                    elif line[0] == 'formula':
                        formula = line[1]

                        thickness = 20  # temporary assignment of the layer thickness
                        sample.addlayer(layer, formula, thickness)  # add layer

                    elif line[0] == 'element':
                        element = line[1]
                        new_element = True
                        sample.structure[layer][element].name = element

                    if new_element:
                        if line[0] == 'molarmass':
                            molarmass= float(line[1])
                            sample.structure[layer][element].molar_mass = molarmass

                        elif line[0] == 'density':
                            density = float(line[1])
                            sample.structure[layer][element].density = density
                        elif line[0] == 'thickness':
                            thickness = float(line[1])
                            sample.structure[layer][element].thickness = thickness
                        elif line[0] == 'roughness':
                            roughness = float(line[1])
                            sample.structure[layer][element].roughness = roughness
                        elif line[0] == 'linkedroughness':
                            linked_roughness = float(line[1])
                            sample.structure[layer][element].linked_roughness = linked_roughness
                        elif line[0] == 'scatteringfactor':
                            line.pop(0)
                            if len(line) == 1:
                                scatteringfactor = line[0]
                            else:
                                scatteringfactor = ast.literal_eval(''.join(line))
                            sample.structure[layer][element].scattering_factor = scatteringfactor

                        elif line[0] == 'polymorph':
                            line.pop(0)
                            polymorph = ast.literal_eval(''.join(line))
                            sample.structure[layer][element].polymorph = polymorph
                        elif line[0] == 'polyratio':
                            line.pop(0)
                            if len(line) == 1:
                                polyratio = 1
                            else:
                                polyratio = ast.literal_eval(''.join(line))
                                polyratio = np.array([float(poly) for poly in polyratio])
                            sample.structure[layer][element].poly_ratio = np.array(polyratio)
                        elif line[0] == 'gamma':
                            gamma = float(line[1])
                            sample.structure[layer][element].gamma = gamma
                        elif line[0] == 'phi':
                            phi = float(line[1])
                            sample.structure[layer][element].phi = phi
                        elif line[0] == 'magdensity':
                            line.pop(0)
                            magdensity = ast.literal_eval(''.join(line))
                            if len(magdensity) != 0:
                                magdensity = np.array([float(mag) for mag in magdensity])
                            sample.structure[layer][element].mag_density = magdensity

                        elif line[0] == 'magscatteringfactor':
                            line.pop(0)
                            if len(line)==0:
                                magscatteringfactor = None
                            else:
                                magscatteringfactor = ast.literal_eval(''.join(line))
                            sample.structure[layer][element].mag_scattering_factor = magscatteringfactor
                        elif line[0] == 'position':
                            position = int(line[1])
                            sample.structure[layer][element].position = position
                            new_element = False  # make sure we do not enter this if statement

                elif experimental_data:

                    if '=' not in line:
                        raise SyntaxError('Data file is improperly initialized.')
                    line.remove('=')  # removes the equal sign

                    info = line[0]  # data identifier
                    data = line[1]  # data value

                    # retrieves the datasetnumber
                    if info == 'datasetnumber':
                        if first:
                            data = int(data)
                            Sinfo[idx]['dataNumber'] = data
                            first = False
                        else:
                            Sinfo[idx]['scanNumber'] = scan_number
                            Sinfo[idx]['scanType'] = scanType
                            Sinfo[idx]['angle'] = angle
                            Sinfo[idx]['energy'] = energy
                            Sinfo[idx]['polarization'] = polarization
                            Sinfo[idx]['numberPoints'] = numberPoints

                            # Attaches appropriate dataset for different scan types
                            if scanType == 'Energy':
                                Sscan.append([x_axis, y_axis, energy_qz])
                            elif scanType == 'Reflectivity':
                                Sscan.append([x_axis, y_axis])


                            # resets all parameters when new scan
                            Sinfo.append(createNewDict())
                            x_axis = list()
                            y_axis = list()
                            energy_qz = list()
                            scan_number = None
                            scanType = None
                            angle = None
                            energy = None
                            numberPoints = None
                            idx = idx + 1
                            data = int(data)
                            Sinfo[idx]['dataNumber'] = data

                    # retrieves the data set title
                    if info == 'datasettitle':
                        scan_number, scanType, angle = getScanInfo(data)
                        if angle != None:
                            angle = float(angle)

                    # sets parameters based on scan type
                    if scanType == 'Energy':
                        if info == 'datasetenergy':
                            energy = float(data)
                        if info == 'polarization':
                            polarization = data
                        if info == 'datasetpoints':
                            numberPoints = int(data)
                        if info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_eng':
                            data = float(data)
                            x_axis.append(data)
                        if info == 'dataset_qz':
                            data = float(data)
                            energy_qz.append(data)
                    elif scanType == 'Reflectivity':
                        if info == 'datasetenergy':
                            energy = float(data)
                        elif info == 'polarization':
                            polarization = data
                        elif info == 'datasetpoints':
                            numberPoints = int(data)
                        elif info == 'dataset_qz':
                            data = float(data)
                            x_axis.append(data)
                        elif info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        elif info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)

                elif simulation:

                    if '=' not in line:
                        raise SyntaxError('Data file is improperly initialized.')
                    line.remove('=')  # removes the equal sign

                    info = line[0]  # data identifier
                    data = line[1]  # data value

                    # retrieves the datasetnumber
                    if info == 'datasetnumber':

                        if first:
                            data = int(data)
                            SimInfo[idx]['dataNumber'] = data
                            first = False
                        else:
                            SimInfo[idx]['scanNumber'] = scan_number
                            SimInfo[idx]['scanType'] = scanType
                            SimInfo[idx]['angle'] = angle
                            SimInfo[idx]['energy'] = energy
                            SimInfo[idx]['polarization'] = polarization
                            SimInfo[idx]['numberPoints'] = numberPoints

                            # Attaches different elements for different scan types
                            if scanType == 'Energy':
                                SimScan.append([x_axis, y_axis, energy_qz])
                            elif scanType == 'Reflectivity':
                                SimScan.append([x_axis, y_axis])


                            # resets all parameters when new scan
                            SimInfo.append(createNewDict())
                            x_axis = list()
                            y_axis = list()
                            energy_qz = list()
                            scan_number = None
                            scanType = None
                            angle = None
                            energy = None
                            numberPoints = None
                            idx = idx + 1
                            data = int(data)
                            SimInfo[idx]['dataNumber'] = data

                    # retrieves the data set title
                    if info == 'datasettitle':
                        scan_number, scanType, angle = getScanInfo(data)
                        if angle != None:
                            angle = float(angle)

                    # sets parameters based on scan type
                    if scanType == 'Energy':
                        if info == 'datasetenergy':
                            energy = float(data)
                        if info == 'polarization':
                            polarization = data
                        if info == 'datasetpoints':
                            numberPoints = int(data)
                        if info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)
                        if info == 'dataset_eng':
                            data = float(data)
                            x_axis.append(data)
                        if info == 'dataset_qz':
                            data = float(data)
                            energy_qz.append(data)
                    elif scanType == 'Reflectivity':
                        if info == 'datasetenergy':
                            energy = float(data)
                        elif info == 'polarization':
                            polarization = data
                        elif info == 'datasetpoints':
                            numberPoints = int(data)
                        elif info == 'dataset_qz':
                            data = float(data)
                            x_axis.append(data)
                        elif info == 'dataset_R0':
                            data = float(data)
                            y_axis.append(data)
                        elif info == 'dataset_A':
                            data = float(data)
                            y_axis.append(data)


    if experimental_data:
        Sinfo[idx]['scanNumber'] = scan_number
        Sinfo[idx]['scanType'] = scanType
        Sinfo[idx]['angle'] = angle
        Sinfo[idx]['energy'] = energy
        Sinfo[idx]['polarization'] = polarization
        Sinfo[idx]['numberPoints'] = numberPoints
        if scanType == 'Energy':
            Sscan.append([x_axis, y_axis, energy_qz])
        elif scanType == 'Reflectivity':
            Sscan.append([x_axis, y_axis])

    elif simulation:
        SimInfo[idx]['scanNumber'] = scan_number
        SimInfo[idx]['scanType'] = scanType
        SimInfo[idx]['angle'] = angle
        SimInfo[idx]['energy'] = energy
        SimInfo[idx]['polarization'] = polarization
        SimInfo[idx]['numberPoints'] = numberPoints
        if scanType == 'Energy':
            SimScan.append([x_axis, y_axis, energy_qz])
        elif scanType == 'Reflectivity':
            SimScan.append([x_axis, y_axis])

    f.close()

    return Sinfo, Sscan, SimInfo, SimScan, sample




def selectScan(fname):
    """
    Purpose: Takes in the read in data and plots the data and the simulated data
    :param Sinfo: Scan info
    :param Sscan: Scan data
    :param sample: Data sample for simulation
    :return:
    """
    Sinfo, Sscan, SimInfo, SimScan, sample = ReadDataASCII(fname)

    # Prints out the scans and their information
    header = ['#', 'Scan Type', 'Energy', 'Angle', 'Polarization']
    tab = PrettyTable(header)
    dataNumber = list()
    for scan in Sinfo:
        data = [scan['dataNumber'], scan['scanType'], scan['energy'], scan['angle'], scan['polarization']]
        dataNumber.append(str(scan['dataNumber']))
        tab.add_row(data)
    print(tab)
    val = input('Select scan # you would like to use: ')
    while val in dataNumber:
        # Determines the scan to use based on #
        val = int(val)
        info = Sinfo[val-1]
        simInfo = SimInfo[val-1]
        data = Sscan[val-1]
        simData = SimScan[val-1]

        scan_type = info['scanType']  # determine the scan type
        pol = info['polarization']  # determines the polarization of the scan
        Rdata = data[1]  # retrieves the reflectivity information


        if scan_type == 'Reflectivity':
            E = info['energy']  # retrieves the energy
            qz = np.array(data[0])  # gets momentum transfer of data

            R = np.array(simData[1])
            # Determines if the reflectivity of the data should be calculated
            plt.figure()
            plt.plot(qz, Rdata)
            plt.plot(qz, R)
            plt.legend(['Data', 'Simulation'])
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')



        elif scan_type == 'Energy':
            Theta = info['angle']  # angle of energy scan
            energy = np.array(data[0])  # energy array
            R = np.array(simData[1])

            plt.figure(2)
            plt.plot(energy, Rdata)
            plt.plot(energy, R)
            plt.legend(['Data', 'Simulation'])
            # Again, determines if natural logarithm needs to be calculated
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')



        plt.show()

        val = input('Select scan # you would like to use: ')


def hdf5ToDict(hform):

    mydict = dict()
    for key in hform.keys():
        mydict[key] = dict()
        mydict[key]['Data'] = list(hform[key])

        for attrskey,val in hform[key].attrs.items():
            mydict[key][attrskey] = val

    return mydict

def saveNewFile(fname, info, data_dict, sample):
    """
        Purpose: Take in data and sample model and save it as an hdf5 file
        :param fname: File name with .hdf5 file extension
        :param AScans: Reflectivity scan data from ProcessRXR
        :param AInfo: Reflectivity scan information from ProcessRXR
        :param EScans: Energy scan data from Process RXR
        :param EInfo: Energy scan information from Process RXR
        :param sample: sample object
        :return:
        """

    # Checking if fname already exists
    cwd = os.getcwd()
    path = cwd + '/' + fname
    if os.path.exists(path):
        raise OSError(
            "HDF5 file already exists. To write new file remove the old file from the current working directory.")

    f = h5py.File(fname, 'a')  # create fname hdf5 file

    # creating group that will contain the sample information
    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)
    grp1.attrs['ScalingFactor'] = float(sample.scaling_factor)
    grp1.attrs['BackgroundShift'] = float(sample.background_shift)

    dsLayer = 0
    for my_layer in sample.structure:

        # Layer information
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)

        formula = ''
        for ele in list(my_layer.keys()):
            stoich = my_layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)
        layer.attrs['Formula'] = formula

        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            # Element information
            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
            element.attrs['LinkedRoughness'] = my_layer[ele].linked_roughness
            element.attrs['PolyRatio'] = my_layer[ele].poly_ratio
            element.attrs['Polymorph'] = my_layer[ele].polymorph
            element.attrs['Gamma'] = my_layer[ele].gamma
            element.attrs['Phi'] = my_layer[ele].phi

            # May be changed in the future as layermagnetized also contains this information
            # Original implemented to avoid problem of trying to load in the magnetic data that does not exist
            if len(my_layer[ele].mag_density) == 0:
                element.attrs['Magnetic'] = False
            else:
                element.attrs['Magnetic'] = True
                element.attrs['MagDensity'] = my_layer[ele].mag_density
                element.attrs['MagScatteringFactor'] = my_layer[ele].mag_scattering_factor

            element.attrs['ScatteringFactor'] = my_layer[ele].scattering_factor
            element.attrs['Position'] = my_layer[ele].position

        dsLayer = dsLayer + 1

    # Loading in the experimental data and simulated data
    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    for inf in info:
        name = inf[2]

        if inf[1] == 'Energy':
            energy = data_dict[name]['Energy']
            angle = data_dict[name]['Angle']
            polarization = data_dict[name]['Polarization']
            datasetpoints = data_dict[name]['DataPoints']
            dsNum = data_dict[name]['DatasetNumber']
            data = data_dict[name]['Data']
            theta = data[1]
            E = data[3]

            E, R = sample.energy_scan(angle, E)
            R = R[polarization]

            sim = np.array([data[0], data[1], R, data[3]])
            dat = np.array([data[0], data[1], data[2], data[3]])

            dset = grpE.create_dataset(name, data=dat)
            dset1 = subE.create_dataset(name, data=sim)

            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = polarization
            dset1.attrs['Polarization'] = polarization

            dset.attrs['Angle'] = float(angle)
            dset1.attrs['Angle'] = float(angle)

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)

        elif inf[1] == 'Reflectivity':

            energy = data_dict[name]['Energy']
            polarization = data_dict[name]['Polarization']
            datasetpoints = data_dict[name]['DataPoints']
            dsNum = data_dict[name]['DatasetNumber']
            data = data_dict[name]['Data']
            dat = np.array([data[0], data[1], data[2]])
            qz, R = sample.reflectivity(energy, data[0])
            R = R[polarization]

            sim = np.array([data[0],data[1],R])

            dset = grpR.create_dataset(name, data=dat)
            dset1 = subR.create_dataset(name, data=sim)

            dset.attrs['Energy'] = float(energy)
            dset1.attrs['Energy'] = float(energy)

            dset.attrs['Polarization'] = str(polarization)
            dset1.attrs['Polarization'] = str(polarization)

            dset.attrs['DataPoints'] = int(datasetpoints)
            dset1.attrs['DataPoints'] = int(datasetpoints)

            dset.attrs['DatasetNumber'] = int(dsNum)
            dset1.attrs['DatasetNumber'] = int(dsNum)

    f.close()
    return

if __name__ == "__main__":

    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density=5.12, roughness=4)
    sample.addlayer(1, 'SrTiO3', 5.28602, density=5.12, roughness=2)

    sample.addlayer(2, 'LaMnO3', 18.84,density = 5, roughness=3.77)
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [0.9,0.1], sf =['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0.00023,0], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 7.73, roughness= 4.54)
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [0.9, 0.1], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.001, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 3.51, roughness=0.788804)
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [0.9, 0.1], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.024, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 2.43, roughness=1.3591)
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 3.86, roughness=0.8215)
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.5, 0.5], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.005, 0], ['Co', 'Ni'])

    sample.addlayer(7, 'CCO', 4, density = 2.5, roughness = 2)

    fname = "Pim10uc.h5"
    fnew = 'test.h5'
    info, data_dict, sim_dict=ReadDataHDF5(fname)
    #print(len(data_dict['59_E429.58_Th5.0_S']['Data'][3]))

    #WriteSampleHDF5(fname, sample)




