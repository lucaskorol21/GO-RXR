from data_structure import *
from material_structure import *
from GlobalOptimization import *
import matplotlib.pyplot as plt
import numpy as np
from numba import *
from scipy import interpolate
from scipy.fft import fft, fftfreq,fftshift, rfft, irfft
from scipy.interpolate import UnivariateSpline
from scipy import signal
import time

import asyncio
from tkinter import *
from tkinter import ttk
import multiprocessing as mp

import functools



def getInfo(val):
    data = val[0]
    data_dict = val[1]
    sim_dict = val[2]

    scans = list()
    # select the scan you would like to view
    #   provide option if they want to skip this step and simply just select the scans
    select_scan = 'Y'
    while select_scan.upper() == 'Y' or select_scan.upper() == 'YES':
        scan_number = input("Select which scan you would like to view? ")
        if scan_number.upper() == "EXIT":
            exit()
        while scan_number not in data[:, 0]:
            scan_number = input("Scan number not in dataset. Please select another scan to view: ")
            if scan_number.upper() == 'EXIT':
                exit()

        # Plotting the Scan
        name = data[int(scan_number) - 1][2]
        scanType = data[int(scan_number) - 1][1]

        dat = data_dict[name]
        sim = sim_dict[name]
        print(name, dat)
        pol = dat.attrs['Polarization']
        print(scanType)
        if scanType.upper() == 'REFLECTIVITY':
            temp_data = list(dat)
            qz_data = temp_data[0]
            Rdata = temp_data[2]

            temp_sim = list(sim)
            qz_sim = temp_sim[0]
            Rsim = temp_sim[2]
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                plt.figure()
                plt.plot(qz_data, Rdata)
                plt.plot(qz_sim, Rsim)
                plt.suptitle('Dataset ' + name + " : " + "Reflectivity Scan")
                plt.xlabel('Momentum Transfer, qz (A^{-1})')
                plt.ylabel('Reflectivity')
                plt.yscale('log')
                plt.show()
            else:

                plt.figure()
                plt.plot(qz_data, Rdata)
                plt.plot(qz_sim, Rsim)
                plt.suptitle('Dataset ' + name + " : " + "Reflectivity Scan")
                plt.xlabel('Momentum Transfer, qz (A^{-1})')
                plt.ylabel('Reflectivity')
                plt.show()

        if scanType.upper() == 'ENERGY':
            print('bye')



        # Questions for the user
        my_scan = input('Would you like to add scan ' + str(scan_number) + ' to data optimization (y/n)?')
        if my_scan.upper() == 'EXIT':
            exit()

        while my_scan.upper() != 'Y' and my_scan.upper() != 'YES' and my_scan.upper() != 'N' and my_scan.upper() != 'NO':
            my_scan = input('Please select (y/n). If you want to exit please type \'exit\': ')
            if my_scan.upper() == 'EXIT':
                exit()

        if my_scan.upper() == 'Y' or my_scan.upper == 'YES':
            scans.append(int(scan_number))
        select_scan = input('Would you like to select another scan (y/n)?')
        if select_scan.upper() != 'Y' and select_scan.upper() != 'YES' and select_scan.upper() != 'N' and select_scan.upper() != 'NO':
            select_scan = input('Please select (y/n). If you want to exit please type \'exit\': ')
        if select_scan.upper() == 'EXIT':
            exit()


def createTable(data):
    # View the different scans
    ws = Tk()
    ws.title('PythonGuides')
    ws.geometry('700x250')
    ws['bg'] = '#AC99F2'

    data_frame = Frame(ws)
    data_frame.pack()

    # scrollbar
    data_scroll = Scrollbar(data_frame)
    data_scroll.pack(side=RIGHT, fill=Y)

    # data_scroll = Scrollbar(data_frame, orient='horizontal')
    # data_scroll.pack(side=BOTTOM, fill=X)

    my_data = ttk.Treeview(data_frame, yscrollcommand=data_scroll.set, xscrollcommand=data_scroll.set)
    my_data.pack()

    # data_scroll.config(command=my_data.yview)
    data_scroll.config(command=my_data.xview)

    # define our column

    my_data['columns'] = ('Scan Number', 'Scan Type', 'Scan Name')

    # format our column
    my_data.column("#0", width=0, stretch=NO)
    my_data.column("Scan Number", anchor=CENTER, width=160)
    my_data.column("Scan Type", anchor=CENTER, width=160)
    my_data.column("Scan Name", anchor=CENTER, width=320)

    # Create Headings
    my_data.heading("#0", text="", anchor=CENTER)
    my_data.heading("Scan Number", text="Scan Number", anchor=CENTER)
    my_data.heading("Scan Type", text="Scan Type", anchor=CENTER)
    my_data.heading("Scan Name", text="Name", anchor=CENTER)

    for idx in range(len(data)):
        my_data.insert(parent='', index='end', iid=idx, text='',
                       values=(data[idx][0], data[idx][1], data[idx][2]))

    my_data.pack()
    ws.mainloop()

if __name__ == '__main__':
    # Define the sample model
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
    
    # save new sample model to current file
    fname = 'Pim10uc.h5'
    #WriteSampleHDF5(fname, sample)


    # Load in the sample information
    #fname = 'Pim10uc.h5'
    sample = ReadSampleHDF5(fname)

    # load in the data and simulation data
    f, data, data_dict, sim_dict = ReadDataHDF5(fname)  # file must remain open as we process the dataset


    # Running the two tasks at once
    f1 = functools.partial(createTable, data)
    f2 = functools.partial(getInfo, [data, data_dict, sim_dict])

    p1 = mp.Process(target=f1)
    p1.start()

    p2 = mp.Process(target=getInfo([data, data_dict, sim_dict]))
    p2.start()

    p2.join()
    p1.join()




    f.close()





