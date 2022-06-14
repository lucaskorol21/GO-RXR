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
    grp1.attrs['Links'] = np.array(sample.link)

    # Retrieve the sample information of each layer
    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)
        formula = ''

        # retrieve the sample information for each element in each layer
        for ele in list(my_layer.keys()):
            element = layer.create_group(ele)

            # reconstructs the formula
            stoich = int(my_layer[ele].stoichiometry)
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)

            element.attrs['MolarMass'] = my_layer[ele].molar_mass
            element.attrs['Density'] = my_layer[ele].density
            element.attrs['Thickness'] = my_layer[ele].thickness
            element.attrs['Roughness'] = my_layer[ele].roughness
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

        layer.attrs['Formula'] = formula
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
        #dset1 = subR.create_dataset(name,dset.shape, dtype=dset.dtype)
        dset1 = subR.create_dataset(name, data=sim)

        dset.attrs['Energy'] = float(energy)
        dset1.attrs['Energy'] = float(energy)

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
                #dset1 = subR.create_dataset(name, dset.shape, dtype=dset.dtype)

                dset.attrs['Energy'] = float(energy)
                dset1.attrs['Energy'] = float(energy)

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
        #dset1 = subE.create_dataset(name, dset.shape, dtype=dset.dtype)
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

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
                qz = EScans[i][:,2]
                E = EScans[i][:,0]
                A = (EScans[i-1][:,2]-EScans[i][:,3])/(EScans[i-1][:,2]+EScans[i][:,3])
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
                #dset1 = subE.create_dataset(name, dset.shape, dtype=dset.dtype)
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

                dsNum = dsNum + 1

    f.close()

def ReadSampleHDF5(fname):
    f = h5py.File(fname, 'r')

    S = f['Sample']

    # Sample Information
    m = int(S.attrs['NumberLayers'])
    sample = slab(m)
    sample.poly_elements = ast.literal_eval(S.attrs['PolyElements'])
    sample.mag_elements = ast.literal_eval(S.attrs['MagElements'])
    sample.layer_magnetized = S.attrs['LayerMagnetized']


    for lay_key in S.keys():
        layer = S[lay_key]
        formula = layer.attrs['Formula']
        lay_num = int(layer.attrs['LayerNumber'])
        sample.addlayer(lay_num, formula,20, density=1)  #pre-initialize parameters to random numbers for each layer
        for ele_key in layer.keys():
            element = layer[ele_key]
            sample.structure[lay_num][ele_key].molar_mass = element.attrs['MolarMass']
            sample.structure[lay_num][ele_key].density = element.attrs['Density']
            sample.structure[lay_num][ele_key].thickness = element.attrs['Thickness']
            sample.structure[lay_num][ele_key].roughness = element.attrs['Roughness']
            sample.structure[lay_num][ele_key].poly_ratio = element.attrs['PolyRatio']
            sample.structure[lay_num][ele_key].polymorph = element.attrs['Polymorph']
            sample.structure[lay_num][ele_key].gamma = element.attrs['Gamma']
            sample.structure[lay_num][ele_key].phi = element.attrs['Phi']
            if element.attrs['Magnetic']:
                sample.structure[lay_num][ele_key].mag_density = element.attrs['MagDensity']
                sample.structure[lay_num][ele_key].mag_scattering_factor = element.attrs['MagScatteringFactor']

            sample.structure[lay_num][ele_key].scattering_factor = element.attrs['ScatteringFactor']
            sample.structure[lay_num][ele_key].position = element.attrs['Position']

    sample.link = S.attrs['Links']
    f.close()
    return sample

def ReadDataHDF5(fname):
    f = h5py.File(fname, 'r')
    experiment = f['Experimental_data']
    simulated = f['Simulated_data']

    RS = experiment['Reflectivity_Scan']
    SimR = simulated['Reflectivity_Scan']

    ES = experiment['Energy_Scan']
    SimE = simulated['Energy_Scan']

    data = list()
    for Rkey in RS.keys():
        mydata = RS[Rkey]
        Rdat = [int(mydata.attrs['DatasetNumber']),'Reflectivity', Rkey]
        data.append(Rdat)
    for Ekey in ES.keys():
        mydata = ES[Ekey]
        Edat = [int(mydata.attrs['DatasetNumber']),'Energy', Ekey]
        data.append(Edat)

    data = np.array(data)
    sort_idx = np.argsort(data[:,0].astype(int))
    data = data[sort_idx]
    header = ['#', 'Scan Type', 'Scan']
    tab = PrettyTable(header)
    tab.add_rows(data)
    print(tab)

    val = input('Choose scan # you want to use: ')
    while val in data[:,0]:
        idx = int(val)-1
        myScan = data[idx]
        scanType = myScan[1]
        scanName = myScan[2]

        if scanType == 'Reflectivity':
            Rdata = list(RS[scanName])
            Rsim = list(SimR[scanName])
            qz = Rdata[0]
            R = Rdata[2]
            Rs = Rsim[2]
            plt.figure()
            plt.plot(qz, R)
            plt.plot(qz, Rs)
            plt.yscale('log')
            if not(RS[scanName].attrs['Polarization'] == 'AL') and not(RS[scanName].attrs['Polarization'] == 'AC'):
                plt.yscale('log')

            plt.legend(['Data', 'Simulation'])

        elif scanType == 'Energy':
            Edata = list(ES[scanName])
            ESim = list(SimE[scanName])
            E = Edata[3]
            R = Edata[2]
            Rs = ESim[2]
            plt.figure()
            plt.plot(E,R)
            plt.plot(E, Rs)
            plt.legend(['Data','Simulation'])
        plt.show()
        val = input('Choose scan # you want to use: ')


    f.close()

def WriteSampleHDF5(fname, sample):
    with h5py.File(fname, "a") as f:
        del f['Sample']
    f.close()

    f = h5py.File(fname, "a")
    grp1 = f.create_group('Sample')
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)
    grp1.attrs['Links'] = np.array(sample.link)

    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)
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

        layer.attrs['Formula'] = formula
        dsLayer = dsLayer + 1

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f["Simulated_data"]

    grpR = grp2["Reflectivity_Scan"]
    grpE = grp2["Energy_Scan"]

    # start of loading data
    dsNum = 1
    for key in list(grpR.keys()):
        dset = grpR[key]
        Rdata = list(dset)
        qz = np.array(Rdata[0])
        E = float(dset.attrs['Energy'])
        polarization = dset.attrs['Polarization']

        qz, R = sample.reflectivity(E, qz)
        Rdata[2] = R[polarization]
        dset[...] = Rdata

    for key in list(grpE.keys()):
        dset = grpE[key]
        Edata = list(dset)
        E = Edata[3]
        theta = dset.attrs['Angle']
        polarization = dset.attrs['Polarization']
        E, R = sample.energy_scan(theta, E)
        Edata[2] = R[polarization]
        dset[...] = Edata

    f.close()

def WriteSampleASCII(file,sample):

    file.write("# Structure \n")
    n = len(sample.structure)
    link = sample.link

    file.write("numberlayers = %s \n" % str(n))
    file.write("polyelements = %s \n" % str(sample.poly_elements))
    file.write("magelements = %s \n" % str(sample.mag_elements))

    layermagnetized = sample.layer_magnetized
    file.write("layermagnetized = %s \n" % layermagnetized)

    link = sample.link
    file.write("link =  %s \n\n" % link)

    num_lay = 0
    molarmass = 0
    density = 0
    thickness = 0
    roughness = 0
    polymorph = 0
    poly_ratio = 0
    gamma = 0
    phi = 0
    mag_density = 0
    scattering_factor = 0
    position = 0
    for layer in sample.structure:

        file.write("layer = %s \n" % str(num_lay))
        formula = ''

        # retieves the chemical formula
        for ele in layer.keys():

            # retrieve the chemical formula
            stoich = layer[ele].stoichiometry
            if stoich == 1:
                formula = formula + ele
            else:
                formula = formula + ele + str(stoich)


        file.write("formula = %s \n\n" % formula)

        for ele in layer.keys():
            file.write("element = %s \n" % ele)
            file.write("molarmass = %f \n" % layer[ele].molar_mass)
            file.write("density = %f \n" % layer[ele].density)
            file.write("thickness = %f \n" % layer[ele].thickness)
            file.write("roughness = %f \n" % layer[ele].roughness)

            poly_names = ''
            poly_ratio = ''
            sf = ''

            scatfact = layer[ele].scattering_factor
            sf = layer[ele].scattering_factor
            file.write("scatteringfactor = %s \n" % sf)

            poly_names = layer[ele].polymorph
            poly_ratio = layer[ele].poly_ratio
            if type(poly_ratio) != int:
                poly_ratio = [str(poly) for poly in poly_ratio]


            file.write("polymorph = %s \n" % poly_names)
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


    file.write("# Experimental_Data \n")

    dsNum = 1
    for i in range(len(AScans)):
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
        for j in range(len(AScans[i][:, 0])):
            file.write("dataset_qz = %f \n" % AScans[i][j][2])
            file.write("dataset_R0 = %e \n" % AScans[i][j][3])
            # file.write("dataset_eng = %f \n" % AScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):
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
                for j in range(len(AScans[i][:, 0])):
                    file.write("dataset_qz = %f \n" % AScans[i][j][2])
                    # print(AScans[i-1][j][3]+AScans[i][j][3])
                    file.write("dataset_A = %e \n" % (
                                (AScans[i - 1][j][3] - AScans[i][j][3]) / (AScans[i - 1][j][3] + AScans[i][j][3])))
                    # file.write("dataset_eng = %f \n" % AScans[i][j][0])
                file.write("\n\n")
                dsNum = dsNum + 1

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


                for j in range(len(EScans[i][:, 0])):
                    file.write("dataset_qz = %f \n" % EScans[i][j][2])
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

    file.write("# Simulation \n")
    qz = list()
    dsNum = 1
    polarization = 'S'
    for i in range(len(AScans)):
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
        qz = AScans[i][:,2]
        qz, R = sample.reflectivity(float(AInfo[i][3]), qz)
        R = R[polarization]
        for j in range(len(qz)):
            file.write("dataset_qz = %f \n" % qz[j])
            file.write("dataset_R0 = %e \n" % R[j])
            # file.write("dataset_eng = %f \n" % AScans[i][j][0])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (AInfo[i - 1][3] == AInfo[i][3]):
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
                qz = AScans[i][:,2]
                qz, R = sample.reflectivity(float(AInfo[i][3]), qz)
                R = R[polarization]
                for j in range(len(qz)):
                    file.write("dataset_qz = %f \n" % qz[j])
                    # print(AScans[i-1][j][3]+AScans[i][j][3])
                    file.write("dataset_A = %e \n" % R[j])
                    # file.write("dataset_eng = %f \n" % AScans[i][j][0])
                file.write("\n\n")
                dsNum = dsNum + 1

    for i in range(len(EScans)):
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
        E = EScans[i][:,0]
        E, R = sample.energy_scan(float(EInfo[i][4]), E)
        R = R[polarization]
        for j in range(len(E)):
            file.write("dataset_R0 = %e \n" % R[j])
            file.write("dataset_eng = %f \n" % E[j])
        file.write("\n\n")
        dsNum = dsNum + 1

        # write asymmetry if possible
        if i > 0:
            if (abs(float(EInfo[i - 1][3]) - float(EInfo[i][3])) < 0.015 and abs(
                    float(EInfo[i - 1][4]) - float(EInfo[i][4])) < 0.1):
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

                E = EScans[i][:,0]
                E, R = sample.energy_scan(float(EInfo[i][4]), E)
                R = R[polarization]
                for j in range(len(E)):
                    file.write("dataset_A = %e \n" % R[j])
                    file.write("dataset_eng = %f \n" % E[j])

                file.write("\n\n")
                dsNum = dsNum + 1


def WriteDataASCII(fname,AScans,AInfo,EScans,EInfo, sample):

    file = open(fname, "w")
    WriteSampleASCII(file, sample)
    WriteExperimentalDataASCII(file, AScans,AInfo,EScans,EInfo)
    WriteSimulationASCII(file, AScans, AInfo, EScans, EInfo, sample)
    f.close()

def getScanInfo(title):
    """
    Purpose: Retrieves important information in the scan title
    :param title: title/label of the scan
    :return:
    """
    title = title.split('_')

    scanType = None
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

def ConvertASCIItoHDF5(fname):

    Sinfo, Sscan, SimInfo, SimScan, sample = ReadDataASCII(fname)

    if fname.endswith('.all'):
        fname = fname[:-4]
        fname = fname + '.hdf5'
    else:
        raise NameError('File name must be a .all file type')

    cwd = os.getcwd()
    path = cwd + '/' + fname

    if os.path.exists(path):
        raise OSError(
            "HDF5 file already exists. To write new HDF5 file remove the old file from the current working directory.")

    f = h5py.File(fname, 'a')


    grp1 = f.create_group("Sample")
    m = len(sample.structure)
    grp1.attrs['NumberLayers'] = int(m)
    grp1.attrs['PolyElements'] = str(sample.poly_elements)
    grp1.attrs['MagElements'] = str(sample.mag_elements)
    grp1.attrs['LayerMagnetized'] = np.array(sample.layer_magnetized)
    grp1.attrs['Links'] = np.array(sample.link)

    dsLayer = 0
    for my_layer in sample.structure:
        name = "Layer_" + str(dsLayer)
        layer = grp1.create_group(name)
        layer.attrs['LayerNumber'] = int(dsLayer)
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

        layer.attrs['Formula'] = formula
        dsLayer = dsLayer + 1

    h = 4.135667696e-15  # Plank's constant eV*s
    c = 2.99792458e8  # speed of light m/s

    grp2 = f.create_group("Experimental_data")
    grp3 = f.create_group("Simulated_data")

    grpR = grp2.create_group("Reflectivity_Scan")
    subR = grp3.create_group("Reflectivity_Scan")

    grpE = grp2.create_group("Energy_Scan")
    subE = grp3.create_group("Energy_Scan")

    name = ''
    # start of loading data
    dsNum = 1

    for i in range(len(Sinfo)):
        info = Sinfo[i]
        data = Sscan[i]
        simData = SimScan[i]

        scanType = info['scanType']
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
            plt.figure()
            plt.plot(E,Rsim)
            plt.show()
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
                    layermagnetized = []
                    formula = ''
                    my_link = []
                    polyelements = dict()
                    magelements = dict()
                    new_element = False
                    element = ''
                    thickness = 20
                    roughness = 0
                    scatteringfactor = []
                    polymorph = []
                    polyratio = []
                    gamma = 90
                    phi = 90
                    magdensity = []
                    magscatteringfactor = []
                    position = 0

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
                    theta = list()
                    energy_qz = list()
                    scan_number = 0
                    scanType = 0
                    angle = 0
                    energy = 0
                    polarization = 0
                    numberPoints = 0

                    NewScan = False  # determines if we have a new scan
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

                    NewScan = False  # determines if we have a new scan
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
                    elif line[0] == 'link':
                        line.pop(0)
                        line = ''.join(line)
                        link = ast.literal_eval(line)
                        sample.link = link
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

                    #line = line.split()
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
                        #Sinfo[idx]['scanNumber'] = scan_number
                        #Sinfo[idx]['scanType'] = scanType
                        if angle != None:
                            angle = float(angle)
                            #Sinfo[idx]['angle'] = angle

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

                    # line = line.split()
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
                        # Sinfo[idx]['scanNumber'] = scan_number
                        # Sinfo[idx]['scanType'] = scanType
                        if angle != None:
                            angle = float(angle)
                            # Sinfo[idx]['angle'] = angle

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


    #fname = "Pim10uc.hdf5"
    #WriteSampleHDF5(fname, sample)
    #ReadDataHDF5(fname)

