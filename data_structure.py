from SpecFileProcessing import *
from FGT_RXR_Processing import *
import matplotlib.pyplot as plt
from prettytable import PrettyTable

def WriteLucas(fname,AScans,AInfo,EScans,EInfo,header):
    file = open(fname, "w")

    startstr = """ # My header"""


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
    title = title.split('_')

    scanType = None
    angle = None

    scan_number = title[0]
    if title[1] == 'A':
        scanType = 'Reflectivity'
    else:
        scanType = 'Energy'
        angle = title[2].replace('Th','')

    return scan_number, scanType, angle
def createNewDict():
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
    idx = 0  # index of scan
    Sinfo = []
    Sinfo.append(createNewDict())
    Sscan = []
    file = open(fname)
    x_axis = list()
    y_axis = list()
    scan_number = 0
    scanType = 0
    angle = 0
    energy = 0
    polarization = 0
    numberPoints = 0

    NewScan = True
    # Read in each line one at a time
    for line in file:
        if line == "\n":
            if NewScan:
                x_axis = list()
                y_axis = list()
                scan_number = None
                scanType = None
                angle = None
                energy = None
                polarization = None
                numberPoints = None
                Sinfo.append(createNewDict())
                NewScan = False
                Sscan.append([x_axis, y_axis])
                idx = idx + 1
        else:
            NewScan = True
            line = line.split()
            if '=' not in line:
                raise SyntaxError('Data file is improperly initialized.')


            line.remove('=')  # removes the equal sign
            info = line[0]  # data identifier
            data = line[1]  # data value

            if info == 'datasetnumber':
                data = int(data)
                Sinfo[idx]['dataNumber'] = data


            if info == 'datasettitle':
                scan_number, scanType, angle = getScanInfo(data)
                Sinfo[idx]['scanNumber'] = scan_number
                Sinfo[idx]['scanType'] = scanType
                if angle != None:
                    angle = float(angle)
                    Sinfo[idx]['angle'] = angle


            if info == 'datasetenergy':
                data = float(data)
                Sinfo[idx]['energy'] = data
            if info == 'polarization':
                Sinfo[idx]['polarization'] = polarization
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
            if info == 'dataset_eng':
                data = float(data)
                y_axis.append(data)
    return Sscan, Sinfo

def selectScan(Sinfo, Sscan):

    header = ['#', 'Scan Type', 'Energy', 'Angle', 'Polarization']
    tab = PrettyTable(header)
    for scan in Sinfo:
        data = [scan['dataNumber'], scan['scanType'], scan['energy'], scan['angle'], scan['polarization']]
        tab.add_row(data)

    val = input('Select scan # you would like to you')


# The purpose of this set of python functions is to create the data structure that will be in my python
# reflectivity software.

# Initial conversations revealed that HDF5 or ASCII are good options
# Personally, I think that HDF5 is the best option

if __name__ == "__main__":
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
    fname = "FGT-1L.all"
    Sscan, Sinfo = ReadLucasFile(fname)

    selectScan(Sinfo, Sscan)

