from SpecFileProcessing import *
from FGT_RXR_Processing import *
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from material_structure import *
from material_model import *
from time import time
import ast
import h5py


def WriteData(fname,AScans,AInfo,EScans,EInfo, sample):

    file = open(fname, "w")
    sampleFormat(file, sample)
    WriteExperimentalData(file, AScans,AInfo,EScans,EInfo)

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

                    if NewScan:
                        # resets all parameters when new scan
                        Sscan.append([x_axis, y_axis])
                        x_axis = list()
                        y_axis = list()
                        scan_number = None
                        scanType = None
                        angle = None
                        energy = None
                        numberPoints = None
                        Sinfo.append(createNewDict())
                        NewScan = False

                        if not(first):
                            idx = idx + 1
                    else:

                        line = line.split()
                        if '=' not in line:
                            raise SyntaxError('Data file is improperly initialized.')

                        line.remove('=')  # removes the equal sign
                        info = line[0]  # data identifier
                        data = line[1]  # data value

                        # retrieves the datasetnumber
                        if info == 'datasetnumber':
                            if first:
                                first = False
                            else:
                                NewScan = True
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
                            elif info == 'polarization':
                                Sinfo[idx]['polarization'] = data
                            elif info == 'datasetpoints':
                                Sinfo[idx]['numberPoints'] = polarization
                            elif info == 'dataset_qz':
                                data = float(data)
                                x_axis.append(data)
                            elif info == 'dataset_R0':
                                data = float(data)
                                y_axis.append(data)
                            elif info == 'dataset_A':
                                data = float(data)
                                y_axis.append(data)

    # sometimes the data file has a new line at the end of the file and creates too long of a list
    if len(Sscan) != len(Sinfo):
        Sinfo.pop()

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

    # sometimes the data file has a new line at the end of the file and creates too long of a list
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
            qz, R, t, e = sample.reflectivity(E,qz)  # performs reflectivity simulation calculation
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


            plt.figure()
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
    fname = "FGT-1L.all"


    sample = slab(2)

    sample.addlayer(0, 'SrTiO3', 50, density=5.12, roughness=2)

    sample.addlayer(1, 'LaMnO3', 10, density=6.5, roughness=2)
    sample.polymorphous(1, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(1, ['Mn2+', 'Mn3+'], [0.5, 0], ['Co', 'Ni'])



    """
    fnameCorr = "FGT-2L"
    samples = ["FGT-2L", "FGT-1L"]

    fnameCorr = "FGT-2L"
    samples = ["FGT-2L", "FGT-1L"]
    for sam in range(len(samples)):
        EScan, AScan, ECal, Geo, Corr = GetSampleInfo(datadir + fnameCorr, datadir + samples[sam])
        AsData, AsInfo = ProcessRXR(datadir + samples[sam], AScan, ECal, Geo, Corr, "A")
        EsData, EsInfo = ProcessRXR(datadir + samples[sam], EScan, ECal, Geo, Corr, "E")

        #remagxHeader = GetReMagXHeader(names[sam], densities[sam], thicknesses[sam])

        WriteData(samples[sam] + ".all", AsData, AsInfo, EsData, EsInfo, sample)
        print()
    print('########################### DONE!!!! #######################################')


    """
    """
    energy = np.linspace(600,602, num=20)
    num_E = len(energy)


    energy, R = sample.energy_scan(15, energy)


    avg_opt = mean(opt_comp)
    avg_adapt = mean(adapt_comp)
    avg_init = mean(init_comp)
    avg_R = mean(R_compt)
    total = avg_opt + avg_adapt + avg_init + avg_R
    pie_chart = np.array([avg_opt,avg_adapt, avg_init, avg_R])*100/total

    my_trunc = np.trunc(np.array(pie_chart)*100)/100
    my_label = []
    for idx in range(len(my_trunc)):
        label = my_trunc[idx]
        label = str(label) + "%"
        my_label.append(label)

    my_legend = ['Optical','Adaptive','Structure','Reflectivity']
    plt.figure()
    plt.suptitle("Two Layer Sample - Timing Chart")
    plt.pie(pie_chart, labels=my_label)
    plt.legend(my_legend, loc="lower right", bbox_to_anchor=(1.3,-0.15))
    plt.show()

    """
    #Sscan, Sinfo, sample1 = ReadData(fname)

    #sample.plot_density_profile(fig=1)
    #sample1.plot_density_profile(fig=2)
    #plt.show()
    #Sscan, Sinfo = ReadLucasFile(fname)

    #selectScan(Sinfo, Sscan, sample)
    #sampleFormat('testascii.all',sample)

    data = np.loadtxt('energy_test.txt')
    energy = data[:, 0]

    start = time()
    sample.energy_scan(15, energy)
    end = time()
    print('Energy Scan: ', end - start)

    E = 642.2
    data = np.loadtxt('test_example.txt')
    qz = data[:,0]
    start_r = time()
    sample.reflectivity(E, qz)
    end_r = time()
    print('Reflectivity Scan: ', end_r-start_r)