from SpecFileProcessing import *
from FGT_RXR_Processing import *
import matplotlib.pyplot as plt
from prettytable import PrettyTable
import os
from material_structure import *
from material_model import *
from time import *
import ast
import h5py


def WriteDataASCII(fname,AScans,AInfo,EScans,EInfo, sample):

    file = open(fname, "w")
    sampleFormat(file, sample)
    WriteExperimentalData(file, AScans,AInfo,EScans,EInfo)

def WriteDataHDF5(fname, AScans,AInfo, EScans, EInfo, sample):

    cwd = os.getcwd()
    path = cwd + '/' + fname

    if os.path.exists(path):
        raise OSError("HDF5 file already exists. To write new file remove the old file from the current working directory.")

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

    #start of loading data
    dsNum = 1
    for i in range(len(AScans)):
        qz = AScans[i][:,2]
        R0 = AScans[i][:,3]
        energy = float(AInfo[i][3])
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

def sampleFormat(file,sample):

    file.write("# Structure \n")
    n = len(sample.structure)
    link = sample.link

    file.write("numberlayers = %s \n" % str(n))
    file.write("polyelements = %s \n" % str(sample.poly_elements))
    file.write("magelements = %s \n" % str(sample.mag_elements))

    laymag = sample.layer_magnetized
    layermagnetized = ''
    for lm in laymag:
        layermagnetized = layermagnetized + str(lm) + " "
    file.write("layermagnetized = %s \n\n" % layermagnetized)

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


        file.write("formula = %s \n" % formula)



        my_link = ''

        for li in link[num_lay]:
            my_link = my_link + str(li) + " "

        file.write("link = %s \n\n" % my_link)

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

            if type(scatfact) == list:
                for s in scatfact:
                    sf =sf + s + " "
            else:
                sf = scatfact

            file.write("scatteringfactor = %s \n" % sf)

            for poly in layer[ele].polymorph:
                poly_names = poly_names + poly + " "

            if type(layer[ele].poly_ratio) != int:
                for rat in layer[ele].poly_ratio:
                    poly_ratio = poly_ratio + str(rat)+ " "


            file.write("polymorph = %s \n" % poly_names)
            file.write("polyratio = %s \n" % poly_ratio)


            file.write("gamma = %f \n" % layer[ele].gamma)
            file.write("phi = %f \n" % layer[ele].phi)

            mag_density = ''
            if len(layer[ele].mag_density) != 0:
                for md in layer[ele].mag_density:
                    mag_density = mag_density + str(md) + " "

            file.write("magdensity = %s \n" % mag_density)

            sfm = ''
            if layer[ele].mag_scattering_factor != None:
                for magscat in layer[ele].mag_scattering_factor:
                    sfm = sfm + magscat + " "
            file.write("magscatteringfactor = %s \n" % sfm)
            file.write("position = %s \n" % layer[ele].position)
            file.write("\n")

        num_lay = num_lay + 1

def WriteExperimentalData(file, AScans,AInfo,EScans,EInfo):
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
        name = str(EInfo[i][0]) + "_E" + str(round(float(EInfo[i][3]), 2)) + "_Th" + str(
            round(float(EInfo[i][4]), 2)) + "_" + EInfo[i][1]
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

    file.close()

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


def ReadData(fname):
    """
    Purpose: Read .all file that contains sample and experimental information
    :param fname: Name of file contatining related information
    :return: Sscan - list of energy/reflectivity data
             Sinfo - contains a dictionary of relevant scan information
             sample - pre-initialized slab class
    """

    file = open(fname)
    structure = False
    experimental_data = False
    for line in file:
        line = line.split()
        if len(line) != 0:
            if line[0] == '#':
                if line[1] == 'Structure':
                    Structure = True
                    experimental_data = False
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
                    Structure = False
                    experimental_data = True

                    idx = 0  # index of scan
                    Sinfo = []  # will contain information of sample
                    Sinfo.append(createNewDict())
                    Sscan = []

                    # initialization of parameters
                    x_axis = list()
                    y_axis = list()
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
                    elif line[0] == 'layer':
                        layer = int(line[1])
                    elif line[0] == 'formula':
                        formula = line[1]
                    elif line[0] == 'layermagnetized':
                        line.pop(0)
                        for laymag in line:
                            if laymag == 'False':
                                layermagnetized.append(False)
                            elif laymag == 'True':
                                layermagnetized.append(True)
                    elif line[0] == 'link':
                        line.pop(0)
                        for val in line:
                            if val == 'True':
                                my_link.append(True)
                            elif val == 'False':
                                my_link.append(False)
                        thickness = 20
                        sample.addlayer(layer,formula,thickness,link=my_link)
                        sample.layer_magnetized = layermagnetized
                        my_link = []
                        layermagnetized = []


                    elif line[0] == 'element':
                        element = line[1]
                        new_element = True

                    if new_element:
                        if line[0] == 'molarmass':
                            molarmass= float(line[1])
                        elif line[0] == 'density':
                            density = float(line[1])
                        elif line[0] == 'thickness':
                            thickness = float(line[1])
                        elif line[0] == 'roughness':
                            roughness = float(line[1])
                        elif line[0] == 'scatteringfactor':
                            line.pop(0)
                            if len(line) == 1:
                                scatteringfactor = line[0]
                            else:
                                scatteringfactor = line

                        elif line[0] == 'polymorph':
                            line.pop(0)
                            if len(line) == 0:
                                polymorph = []
                            else:
                                polymorph = line

                        elif line[0] == 'polyratio':
                            line.pop(0)
                            if len(line) == 0:
                                polyratio = 1
                            else:
                                for rat in line:
                                    polyratio.append(float(rat))
                        elif line[0] == 'gamma':
                            gamma = line[1]
                        elif line[0] == 'phi':
                            phi = line[1]
                        elif line[0] == 'magdensity':
                            line.pop(0)
                            if len(line) == 0:
                                magdensity = []
                            if len(line) == 1:
                                magdensity = float(line[0])
                            else:
                                for mag in line:
                                    magdensity.append(float(mag))

                        elif line[0] == 'magscatteringfactor':
                            line.pop(0)
                            if len(line)==0:
                                magscatteringfactor = None
                            else:
                                magscatteringfactor = line
                        elif line[0] == 'position':
                            position = int(line[1])

                            # End of arguments to load in
                            sample.structure[layer][element].name = element
                            sample.structure[layer][element].molar_mass = molarmass
                            sample.structure[layer][element].density = density
                            sample.structure[layer][element].thickness = thickness
                            sample.structure[layer][element].roughness = roughness
                            sample.structure[layer][element].poly_ratio = np.array(polyratio)
                            sample.structure[layer][element].polymorph = polymorph
                            sample.structure[layer][element].gamma = gamma
                            sample.structure[layer][element].phi = phi
                            sample.structure[layer][element].mag_density = magdensity
                            sample.structure[layer][element].scattering_factor = scatteringfactor
                            sample.structure[layer][element].mag_scattering_factor = magscatteringfactor
                            sample.structure[layer][element].position = position
                            polyratio = []
                            new_element=False  # make sure we do not enter this if statement


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
                            Sscan.append([x_axis, y_axis])

                            # resets all parameters when new scan
                            Sinfo.append(createNewDict())
                            x_axis = list()
                            y_axis = list()
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
        Sscan.append([x_axis, y_axis])

    return Sscan, Sinfo, sample


def ReadLucasFile(fname):
    """
    Purpose: Read in datafile and transform into a usable format for dta analysis
    :param fname: The sample file name
    :return:
    """
    idx = 0  # index of scan
    Sinfo = []  # will contain information of sample
    Sinfo.append(createNewDict())
    Sscan = []
    file = open(fname)

    # initialization of parameters
    x_axis = list()
    y_axis = list()
    scan_number = 0
    scanType = 0
    angle = 0
    energy = 0
    polarization = 0
    numberPoints = 0

    NewScan = True  # determines if we have a new scan
    # Read in each line one at a time
    for line in file:
        if line == "\n":
            if NewScan:
                # resets all parameters when new scan
                Sscan.append([x_axis, y_axis])
                x_axis = list()
                y_axis = list()
                scan_number = None
                scanType = None
                angle=None
                energy = None
                numberPoints = None
                Sinfo.append(createNewDict())
                NewScan = False
                idx = idx + 1
        else:

            NewScan = True
            line = line.split()
            if '=' not in line:
                raise SyntaxError('Data file is improperly initialized.')


            line.remove('=')  # removes the equal sign
            info = line[0]  # data identifier
            data = line[1]  # data value
            # retrieves the datasetnumber
            if info == 'datasetnumber':
                data = int(data)
                Sinfo[idx]['dataNumber'] = data

            # retrieves the data set title
            if info == 'datasettitle':
                scan_number, scanType, angle = getScanInfo(data)
                Sinfo[idx]['scanNumber'] = scan_number
                Sinfo[idx]['scanType'] = scanType
                if angle != None:
                    angle = float(angle)
                    Sinfo[idx]['angle'] = angle

            # sets parameters based on scan type
            if scanType == 'Energy':
                if info == 'datasetenergy':
                    data = float(data)
                    Sinfo[idx]['energy'] = data
                if info == 'polarization':
                    Sinfo[idx]['polarization'] = data
                if info == 'datasetpoints':
                    Sinfo[idx]['numberPoints'] = polarization
                if info == 'dataset_R0':
                    data = float(data)
                    y_axis.append(data)
                if info == 'dataset_A':
                    data = float(data)
                    y_axis.append(data)
                if info == 'dataset_eng':
                    data = float(data)
                    x_axis.append(data)
            elif scanType == 'Reflectivity':
                if info == 'datasetenergy':
                    data = float(data)
                    Sinfo[idx]['energy'] = data
                if info == 'polarization':
                    Sinfo[idx]['polarization'] = data
                if info == 'datasetpoints':
                    Sinfo[idx]['numberPoints'] = polarization
                if info == 'dataset_qz':
                    data = float(data)
                    x_axis.append(data)
                if info == 'dataset_R0':
                    data = float(data)
                    y_axis.append(data)
                if info == 'dataset_A':
                    data = float(data)
                    y_axis.append(data)

    if len(Sscan) != len(Sinfo):
        Sinfo.pop()

    return Sscan, Sinfo



def selectScan(Sinfo, Sscan, sample):
    """
    Purpose: Takes in the read in data and plots the data and the simulated data
    :param Sinfo: Scan info
    :param Sscan: Scan data
    :param sample: Data sample for simulation
    :return:
    """
    sample.plot_density_profile(fig=1)  # plots the sample density profile
    # Prints out the scans and their information
    header = ['#', 'Scan Type', 'Energy', 'Angle', 'Polarization']
    tab = PrettyTable(header)
    for scan in Sinfo:
        data = [scan['dataNumber'], scan['scanType'], scan['energy'], scan['angle'], scan['polarization']]
        tab.add_row(data)
    print(tab)
    val = input('Select scan # you would like to use: ')
    val = int(val)
    while val != 0:
        # Determines the scan to use based on #

        info = Sinfo[val-1]
        data = Sscan[val-1]

        scan_type = info['scanType']  # determine the scan type
        pol = info['polarization']  # determines the polarization of the scan
        Rdata = data[1]  # retrieves the reflectivity information


        if scan_type == 'Reflectivity':
            E = info['energy']  # retrieves the energy
            qz = np.array(data[0])  # gets momentum transfer of data

            start = time()
            qz, R = sample.reflectivity(E,qz)  # performs reflectivity simulation calculation
            end = time()

            print(end-start)
            # Determines if the reflectivity of the data should be calculated
            plt.figure()
            plt.plot(qz, Rdata)
            plt.plot(qz, R[pol])
            plt.legend(['Data', 'Simulation'])
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')



        elif scan_type == 'Energy':
            Theta = info['angle']  # angle of energy scan
            energy = np.array(data[0])  # energy array

            start = time()
            energy, R = sample.energy_scan(Theta, energy)  # simulated energy scan
            #energy, R = sample.energy_scan_multi(Theta, energy)  # simulated energy scan
            end = time()
            print(end-start)


            plt.figure(2)
            plt.plot(energy, Rdata)
            plt.plot(energy, R[pol])
            plt.legend(['Data', 'Simulation'])
            # Again, determines if natural logarithm needs to be calculated
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')



        plt.show()

        val = input('Select scan # you would like to use: ')
        val = int(val)





if __name__ == "__main__":

    """

    fname = "FGT-1L.all"


    sample = slab(2)

    sample.addlayer(0, 'SrTiO3', 50, density=5.12, roughness=2)

    sample.addlayer(1, 'LaMnO3', 10, density=6.5, roughness=2)
    sample.polymorphous(1, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [0.5, 0], ['Co', 'Ni'])

    fnameCorr = "FGT-2L"
    samples = ["FGT-2L", "FGT-1L"]

    fnameCorr = "FGT-2L"
    samples = ["FGT-2L", "FGT-1L"]
    for sam in range(len(samples)):
        EScan, AScan, ECal, Geo, Corr = GetSampleInfo(datadir + fnameCorr, datadir + samples[sam])
        AsData, AsInfo = ProcessRXR(datadir + samples[sam], AScan, ECal, Geo, Corr, "A")
        EsData, EsInfo = ProcessRXR(datadir + samples[sam], EScan, ECal, Geo, Corr, "E")

        #remagxHeader = GetReMagXHeader(names[sam], densities[sam], thicknesses[sam])

        #WriteDataASCII(samples[sam] + ".all", AsData, AsInfo, EsData, EsInfo, sample)
        start = time()
        WriteDataHDF5(samples[sam] + ".hdf5", AsData, AsInfo, EsData, EsInfo, sample)
        end = time()
        print(end-start)

        print()
    print('########################### DONE!!!! #######################################')
    """

    fname = "FGT-1L.hdf5"
    #f = h5py.File(fname, 'r')
    #print(f['Simulated_data/Energy_Scan'].keys())
    #ReadDataHDF5(fname)

    sample = slab(2)

    sample.addlayer(0, 'SrTiO3', 50, density=5.12, roughness=2)

    sample.addlayer(1, 'LaMnO3', 20, density=5.5, roughness=2)
    sample.polymorphous(1, 'Mn', ['Mn2+', 'Mn3+'], [0.5, 0.5], sf=['Mn', 'Fe'])
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [0.5, 0], ['Co', 'Ni'])
    WriteSampleHDF5(fname, sample)

    ReadDataHDF5(fname)
    #Sscan, Sinfo, sample1 = ReadData(fname)
    #selectScan(Sinfo, Sscan, sample)
    #sample.plot_density_profile(fig=1)
    #sample1.plot_density_profile(fig=2)
    #plt.show()
