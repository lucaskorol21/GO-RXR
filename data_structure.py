from SpecFileProcessing import *
from FGT_RXR_Processing import *
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from material_structure import *
from material_model import *
import h5py

def WriteLucas(fname,AScans,AInfo,EScans,EInfo,header):
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
    file = open(fname, "w")

    startstr = """ # My header"""

    # Current version does not use header
    if (header == ""):
        file.write(startstr)
    elif (header == 'None'):
        print('No header')
    else:
        file.write(header)

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

            qz, R, t, e = sample.reflectivity(E,qz)  # performs reflectivity simulation calculation

            # Determines if the reflectivity of the data should be calculated
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                is_all_zero = np.all((R[pol] == 0))
                Rtest = R[pol]
                Rdata = np.log10(Rdata)
                if not(is_all_zero):
                    Rtest = np.log10(R[pol])

            else:
                Rtest = R[pol]

            plt.figure()
            plt.plot(qz, Rdata)
            plt.plot(qz, Rtest)
            plt.legend(['Data','Simulation'])

        elif scan_type == 'Energy':
            Theta = info['angle']  # angle of energy scan
            energy = np.array(data[0])  # energy array
            energy, R = sample.energy_scan(Theta, energy)  # simulated energy scan

            # Again, determines if natural logarithm needs to be calculated
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                is_all_zero = np.all((R[pol] == 0))
                Rtest = R[pol]
                Rdata = np.log10(Rdata)
                if not (is_all_zero):
                    Rtest = np.log10(R[pol])
            else:
                Rtest = R[pol]

            plt.figure()
            plt.plot(energy,Rdata)
            plt.plot(energy, Rtest)
            plt.legend(['Data', 'Simulation'])

        plt.show()
        val = input('Select scan # you would like to use: ')
        val = int(val)


if __name__ == "__main__":
    """
    fnameCorr = "FGT-2L"
    samples = ["FGT-2L", "FGT-1L"]

    names = [["Ge", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Ge5O", "Ge5O3C"],
             ["Ge", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Fe300Ge100Te200Co", "Ge5O", "Ge5O3C"]]

    densities = [[5.323, 7.3, 7.3, 7.3, 5.323, 5.323],
                 [5.323, 7.3, 7.3, 7.3, 5.323, 5.323]]

    thicknesses = [[0, 6, 6, 6, 25, 25],
                   [0, 3, 3, 3, 25, 25]]
    
    for sam in range(len(samples)):
        EScan, AScan, ECal, Geo, Corr = GetSampleInfo(datadir + fnameCorr, datadir + samples[sam])
        AsData, AsInfo = ProcessRXR(datadir + samples[sam], AScan, ECal, Geo, Corr, "A")
        EsData, EsInfo = ProcessRXR(datadir + samples[sam], EScan, ECal, Geo, Corr, "E")

        #remagxHeader = GetReMagXHeader(names[sam], densities[sam], thicknesses[sam])
        header = 'None'
        WriteLucas(samples[sam] + ".all", AsData, AsInfo, EsData, EsInfo, header)
        print()
    print('########################### DONE!!!! #######################################')
    """
    """
    sample =slab(2)

    sample.addlayer(0,'SrTiO3',50, density = 5.12,roughness=2)

    sample.addlayer(1,'LaMnO3', 40, density= 6.5, roughness=2)
    sample.polymorphous(1,'Mn',['Mn2+','Mn3+'],[1,0], sf=['Mn','Fe'])
    sample.magnetization(1,['Mn2+','Mn3+'],[2,0],['Co', 'Ni'])

    fname = "FGT-1L.all"
    Sscan, Sinfo = ReadLucasFile(fname)

    selectScan(Sinfo, Sscan, sample)

    """
    f = h5py.File("mytestfile.hdf5","w")
    dset = f.create_dataset("mydtaset",(100,), dtype='i')

