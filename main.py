import numpy as np
import data_structure as ds
import matplotlib.pyplot as plt
import global_optimization as go
from scipy.signal import *



if __name__ == '__main__':
    fname = "//cabinet/work$/lsk601/My Documents/Masters_Project/SrTiO3-LaMnO3/Pim10uc.all"
    f = "//cabinet/work$/lsk601/My Documents/Masters_Project/SrTiO3-LaMnO3/test_2.txt"
    data_info, data_dict = ds.Read_ReMagX(fname)

    idx = 2
    name = '51_640.2_S'
    #name = data_info[idx][2]
    scan_type = data_info[idx][1]

    data = data_dict[name]['Data']

    temp_data = np.loadtxt(f)

    if scan_type == 'Reflectivity':
        qz = data[0]
        angle = data[1]
        R = data[2]

        #Rf = savgol_filter(np.log10(R), 31, 3)  # window size 51, polynomial order 3

        #Rdiff = savgol_filter(np.diff(np.log(R)),51,5)

        n = len(qz)

        Rf = resample(np.log10(R),n,t=qz)
        Rf = Rf[0]
        print(Rf)

        plt.figure(1)
        plt.plot(qz,Rf)
        plt.plot(qz,np.log10(R))
        #plt.plot(temp_data[:,0], np.log10(temp_data[:,1]))

        plt.show()

        plt.figure(2)
        plt.plot(np.diff(np.log10(temp_data[:,1])))
        #plt.plot(Rdiff)
        plt.show()

        #tv = go.total_variation(np.log10(temp_data[:,1][14:-14]),Rroll[14:-14])

        #print(tv)