from scipy import optimize
from material_structure import *
import numpy as np
from data_structure import *
import matplotlib.pyplot as plt
import time
from tkinter import *
from material_model import *
from tkinter import ttk
import multiprocessing as mp
import sys
import os
from PIL import Image, ImageTk

import functools

def saveScans(data, data_dict, sim_dict, scans, sample):


    dir = 'Plot_Scans'


    for file in os.scandir(dir):
        os.remove(file.path)

    sample.plot_density_profile(save=True, fig=99)

    my_index = list()  # contains the indices of the appropriate scans
    # Finds the indices of each scan
    for s in scans:
        my_index.append((data[:, 0].astype(int).tolist().index(s)))

    my_scans = data[my_index]  # gets all the appropriate scans

    # plot the sample model
    fig_idx = 1
    for s in my_scans:

        name = s[2]  # name of the scan
        scanType = s[1]  # scan type
        scan_num = s[0]  # scan number
        dat = data_dict[name]
        sim = sim_dict[name]
        pol = dat['Polarization']

        if scanType.upper() == 'REFLECTIVITY':
            qz = dat['Data'][0]  # momentum transfer
            R = dat['Data'][2]  # reflectivity
            Rsim = sim['Data'][2]
            plt.figure(fig_idx)
            plt.suptitle('Reflectivity Scan ' + str(scan_num) + ': ' + name)
            plt.plot(qz, R)  # data
            plt.plot(qz, Rsim)  # simulation

            # Change log scale for non-asymmetry scans
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')
                plt.ylabel('log(R)')
            else:
                plt.ylabel('A')
            plt.xlabel('Momentum Transfer, qz (A^{-1})')
            plt.legend(['Experiment', 'Simulation'])

        elif scanType.upper() == 'ENERGY':
            E = dat['Data'][3]  # Energy
            R = dat['Data'][2]  # Reflectivity
            Rsim = sim['Data'][2]
            plt.figure(fig_idx)
            plt.suptitle('Energy Scan' + str(scan_num) + ': ' + name)
            plt.plot(E, R)  # data
            plt.plot(E, Rsim)  # simulation

            # Changes the yscale depending on if scan is an asymmetry scan
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')
                plt.ylabel('log(R)')
            else:
                plt.ylabel('A')
            plt.xlabel('Energy, E (eV)')
            plt.legend(['Experiment', 'Simulation'])

        fig_idx = fig_idx + 1
        saveto = dir +'/' + name + '.png'
        plt.savefig(saveto)
        plt.close()

    return

def plotScansWidget(sample):

    dir = 'Plot_Scans'
    root = Tk()
    root.geometry('900x900')
    root.title('Show Selected Scans')

    # create a notebook
    tabControl = ttk.Notebook(root)

    parameters = ['Thickness', 'Density', 'Roughness', 'Linked Roughness',
                  'Stoichiometry', 'Polymorph', 'Gamma', 'Phi', 'Magnetic Density',
                  'Scattering Factor', 'Magnetic Scattering Factor']
    #root = tk.Tk()
    #root.title('Sample Information')
    #root.geometry("900x900")
    tree = ttk.Treeview(root)
    ttk.Style().configure('Treeview', rowheight=30)

    structure = sample.structure

    # Loop through all the different elements
    for layer in range(len(structure)):
        if layer == 0:
            layer_name = 'Substrate'
        else:
            layer_name = 'Layer ' + str(layer)
        tree.insert("", layer, layer_name, text=layer_name)  # adding the layer to the tree

        elements = list(structure[layer].keys())

        for e in range(len(elements)):
            element = elements[e]
            element_name = element + " " + str(layer)
            tree.insert(layer_name, e, element_name, text=element)

            for param in parameters:
                if param == 'Thickness':
                    thick_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 1, thick_name, text=param)
                    thick_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(thick_name, 1, thick_data_name, text=str(structure[layer][element].thickness) + " (A)")
                elif param == 'Density':
                    density_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 2, density_name, text=param)
                    density_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(density_name, 1, density_data_name,
                                text=str(structure[layer][element].density) + " (mol/cm^3)")
                elif param == 'Roughness':
                    rough_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 3, rough_name, text=param)
                    rough_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(rough_name, 1, rough_data_name, text=str(structure[layer][element].roughness)) + " (A)"
                elif param == 'Linked Roughness':
                    link_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 4, link_name, text=param)
                    link_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(link_name, 1, link_data_name,
                                text=str(structure[layer][element].linked_roughness)) + " (A)"
                elif param == 'Stoichiometry':
                    stoich_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 5, stoich_name, text=param)
                    stoich_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(stoich_name, 1, stoich_data_name, text=structure[layer][element].stoichiometry)
                elif param == 'Polymorph':
                    polymorphs = structure[layer][element].polymorph
                    if len(polymorphs) != 0:
                        polymorphs_name = param + " " + element + " " + str(layer)
                        tree.insert(element_name, 6, polymorphs_name, text=param)
                        for poly in range(len(polymorphs)):
                            poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                            tree.insert(polymorphs_name, poly, poly_name, text=polymorphs[poly])
                            poly_density_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                            poly_data = 'Ratio = ' + str(structure[layer][element].poly_ratio[poly])
                            tree.insert(poly_name, 1, poly_density_name, text=poly_data)
                elif param == 'Gamma':
                    gamma_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 7, gamma_name, text=param)
                    gamma_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(gamma_name, 1, gamma_data_name, text=str(structure[layer][element].gamma) + " degrees")
                elif param == 'Phi':
                    phi_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 8, phi_name, text=param)
                    phi_data_name = param + " " + element + " " + str(layer) + " data"
                    tree.insert(phi_name, 1, phi_data_name, text=str(structure[layer][element].phi) + " degrees")
                elif param == 'Magnetic Density':
                    polymorphs = structure[layer][element].polymorph
                    if len(polymorphs) != 0:
                        magnetic_name = param + " " + element + " " + str(layer)
                        tree.insert(element_name, 9, magnetic_name, text=param)
                        for poly in range(len(polymorphs)):
                            poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                            tree.insert(magnetic_name, poly, poly_name, text=polymorphs[poly])
                            mag_density_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                            poly_data = str(structure[layer][element].mag_density[poly]) + " mol/cm^3"
                            tree.insert(poly_name, 1, mag_density_name, text=poly_data)
                elif param == 'Scattering Factor':
                    polymorphs = structure[layer][element].polymorph
                    scat_name = param + " " + element + " " + str(layer)
                    tree.insert(element_name, 9, scat_name, text=param)
                    if len(polymorphs) != 0:
                        for poly in range(len(polymorphs)):
                            poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                            tree.insert(scat_name, poly, poly_name, text=polymorphs[poly])
                            scat_fact_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                            poly_data = str(structure[layer][element].scattering_factor[poly])
                            tree.insert(poly_name, 1, scat_fact_name, text=poly_data)
                    else:
                        scat_fact_name = param + " " + element + " " + str(layer) + " " + "data"
                        tree.insert(scat_name, 1, scat_fact_name, text=structure[layer][element].scattering_factor)
                elif param == 'Magnetic Scattering Factor':
                    polymorphs = structure[layer][element].polymorph
                    if len(structure[layer][element].mag_density) != 0:
                        mag_scat_name = param + " " + element + " " + str(layer)
                        tree.insert(element_name, 9, mag_scat_name, text=param)
                        if len(polymorphs) != 0:
                            for poly in range(len(polymorphs)):
                                poly_name = param + " " + polymorphs[poly] + " " + str(layer)
                                tree.insert(mag_scat_name, poly, poly_name, text=polymorphs[poly])
                                mag_scat_fact_name = param + " " + polymorphs[poly] + " " + str(layer) + " " + "data"
                                poly_data = str(structure[layer][element].mag_scattering_factor[poly])
                                tree.insert(poly_name, 1, mag_scat_fact_name, text=poly_data)
                        else:
                            scat_fact_name = param + " " + element + " " + str(layer) + " " + "data"
                            tree.insert(mag_scat_name, 1, scat_fact_name,
                                        text=structure[layer][element].mag_scattering_factor)

    tree.pack(expand=True, fill='both')
    tabControl.add(tree, text='Sample Information')
    # Showing the Scans
    """
    dir = 'Plot_Scans'
    
    root = Tk()
    root.geometry('900x900')
    root.title('Show Selected Scans')

    # create a notebook
    tabControl = ttk.Notebook(root)
    """
    idx = 0
    im = list()
    for filename in os.listdir(dir):

        frame = ttk.Frame(tabControl)
        imageFile = dir + '/' + filename
        im.append(ImageTk.PhotoImage(Image.open(imageFile)))
        label = Label(frame, imag=im[idx])
        label.pack()
        frame.pack()

        tabControl.add(frame, text=filename)

        idx = idx + 1


    tabControl.pack()
    root.mainloop()
    return

def plotScans(data, data_dict, sim_dict, scans):
    """
    Purpose: Plot the all the selected scans
    :param scans: An array containing the scans that you want to plot
    :param data: The two dimensional list that contains scan number and name
    :param data_dict: Data dictionary
    :param sim_dict: Simulation dictionary
    :return:
    """

    my_index = list()  # contains the indices of the appropriate scans
    # Finds the indices of each scan
    for s in scans:
        my_index.append((data[:, 0].astype(int).tolist().index(s)))

    my_scans = data[my_index] # gets all the appropriate scans

    # plot the sample model
    fig_idx = 1
    for s in my_scans:

        name = s[2]  # name of the scan
        scanType = s[1]  # scan type
        scan_num = s[0]  # scan number
        dat = data_dict[name]
        sim = sim_dict[name]
        pol = dat['Polarization']

        if scanType.upper() == 'REFLECTIVITY':
            qz = dat['Data'][0]  # momentum transfer
            R = dat['Data'][2]  # reflectivity
            Rsim = sim['Data'][2]
            plt.figure(fig_idx)
            plt.suptitle('Reflectivity Scan ' + str(scan_num) + ': ' + name)
            plt.plot(qz, R)  # data
            plt.plot(qz, Rsim)  # simulation

            # Change log scale for non-asymmetry scans
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')
                plt.ylabel('log(R)')
            else:
                plt.ylabel('A')
            plt.xlabel('Momentum Transfer, qz (A^{-1})')
            plt.legend(['Experiment', 'Simulation'])

        elif scanType.upper() == 'ENERGY':
            E = dat['Data'][3]  # Energy
            R = dat['Data'][2]  # Reflectivity
            Rsim = sim['Data'][2]
            plt.figure(fig_idx)
            plt.suptitle('Energy Scan' + str(scan_num) + ': ' + name)
            plt.plot(E, R)  # data
            plt.plot(E, Rsim)  # simulation

            # Changes the yscale depending on if scan is an asymmetry scan
            if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':
                plt.yscale('log')
                plt.ylabel('log(R)')
            else:
                plt.ylabel('A')
            plt.xlabel('Energy, E (eV)')
            plt.legend(['Experiment', 'Simulation'])

        fig_idx = fig_idx + 1

    plt.show()
    return

def getScans(data, data_dict, sim_dict, queue):
    """
    Purpose: Asks user what scans they would like to see and use for global optimization.
    :param data: array containing scan number, scan type, and scan name
    :param data_dict: data dictionary
    :param sim_dict: simulation dictionary
    :param queue: multiprocessing queue used to store scans
    :return:
    """
    print('SCAN SELECTION FOR GLOBAL OPTIMIZATION \n')
    print('Select an option:')

    stage1 = {'1': 'Select Scan',
              '2': 'Exit/Finished'}

    stage2 = {'1': 'Use Scan',
              '2': 'Choose Different Scan',
              '3': 'Return',
              '4': 'Exit'}

    stage3 = {'1': 'Select Scan Boundaries',
              '2': 'Use Default Boundary and Weights',
              '3': 'Return',
              '4': 'Exit/Finish'}

    stage4 = {'1': 'Select Scan Boundary Weights',
              '2': 'Use Default Boundary Weights',
              '3': 'Return',
              '4': 'Exit'}


    scanNumberList = list(data[:,0])  # get the scan number list

    scans = list()
    boundaries = list()
    weights = list()

    cont = True

    stg1 = True   # stage 1
    stg2 = False  # stage 2
    stg3 = False  # stage 3
    stg4 = False  # weight selection

    scan = 0
    f = functools.partial(plotScans, data, data_dict, sim_dict)
    while cont:

        # Determines if user wants to select a scan
        if stg1:
            for key in stage1.keys():
                val = stage1[key]
                print('\t'+key+': ' + val)
            toggle = input("\n" + "-> ")
            print()
            if toggle != '1' and toggle != '2':
                while toggle != '1' and toggle != '2':
                    toggle = input('Please select one of the provided options: ')

            if toggle == '2':
                cont = False
            elif toggle == '1':
                stg1 = False
                stg2 = True

        # Determine which scans the user wants to use
        elif stg2:
            print('SCAN SELECTION' + '\n')
            scan = input('Please select scan you would like to view: ')
            while scan not in scanNumberList:
                scan = input('Scan number must be found between '+scanNumberList[0]+' and ' + scanNumberList[-1]+': ')
            while scan in scans:
                scan = input('Scan ' + scan + ' already selected! Choose another scan: ')
            print('\n Select an option: ')

            p = mp.Process(target=f, args=([int(scan)],))
            p.start()
            #p.join()
            for key in stage2.keys():
                val = stage2[key]
                print('\t' + key + ': ' + val)
            toggle = input("\n" + "-> ")
            print()
            if toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                    toggle = input('Please select one of the provided options: ')

            if toggle == '1':
                scans.append(scan)
                stg2 = False
                stg3 = True
            elif toggle == '2':
                pass
            elif toggle == '3':
                pass
                stg1 = True
                stg2 = False
            elif toggle == '4':
                cont = False

        # Selecting bounds
        elif stg3:
            scanNumber = data[int(scan) - 1][0]
            scanType = data[int(scan) - 1][1]
            scanName = data[int(scan) - 1][2]

            print('\n Choose boundary options:')
            for key in stage3.keys():
                val = stage3[key]
                print('\t' + key + ': ' + val)
            toggle = input("\n" + "-> ")
            print()

            if toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                    toggle = input('Please select one of the provided options: ')

            if toggle == '1':
                bound = []
                if scanType == 'Reflectivity':
                    qz_up = str(round(data_dict[scanName]['Data'][0][-1],2))
                    qz_lw = str(round(data_dict[scanName]['Data'][0][0],2))
                    bound = input('\n Enter the boundaries for scan ' + scan + ' in terms of qz (' + qz_lw + ', '+qz_up + '): ')

                elif scanType == 'Energy':
                    E_up = str(round(data_dict[scanName]['Data'][3][-1],2))
                    E_lw = str(round(data_dict[scanName]['Data'][3][0],2))
                    bound = input('\n Enter the boundaries for scan ' + scan + ' in terms of qz (' + E_lw + ', ' + E_up + '): ')

                bound = bound.split()
                bound = [ast.literal_eval(bd) for bd in bound]

                prev = 0
                tuple_eval = True
                limit_eval = True
                while tuple_eval or limit_eval:
                    tuple_eval = False
                    limit_eval = False
                    for bd in bound:

                        if type(bd) != tuple:
                            tuple_eval = True
                        elif bd[0] > bd[1] or prev > bd[0]:
                            limit_eval = True
                        prev = bd[1]
                    if tuple_eval:
                        bound = input('Enter the boundaries as tuples, each separated by spaces:')
                        bound = bound.split()
                        bound = [ast.literal_eval(bd) for bd in bound]
                    if limit_eval:
                        bound = input('Make sure boundary limits are in ascending order: ')
                        bound = bound.split()
                        bound = [ast.literal_eval(bd) for bd in bound]

                boundaries.append(bound)
                # go onto next step
                stg3 = False
                stg4 = True
            elif toggle == '2':
                bound = []
                if scanType == 'Reflectivity':
                    qz_up =data_dict[scanName]['Data'][0][-1]
                    qz_lw = data_dict[scanName]['Data'][0][0]
                    bound = [(qz_lw, qz_up)]


                elif scanType == 'Energy':
                    E_up = data_dict[scanName]['Data'][3][-1]
                    E_lw = data_dict[scanName]['Data'][3][0]
                    bound = [(E_lw, E_up)]

                boundaries.append(bound)
                weights.append([1])
                # set boundary and weights to default values
                stg3 = False
                stg1 = True
                p.terminate()
            elif toggle == '3':
                # remove previous selected scan
                scans.pop()
                stg2 = True
                stg3 = False
            elif toggle == '4':
                # set boundary and weights to default and exit
                cont = False

        elif stg4:
            for key in stage4.keys():
                val = stage4[key]
                print('\t' + key + ': ' + val)

            toggle = input("\n" + "-> ")
            print()

            if toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                    toggle = input('Please select one of the provided options: ')

            if toggle == '1':
                num_weights = len(boundaries[-1])
                weight = input('Input weights for the boundaries selected: ')
                weight = weight.split()
                while len(weight) != num_weights:
                    weight = input('Select '+ str(num_weights)+' weights, each separated by a space: ')
                    weight = weight.split()

                weight = [float(w) for w in weight]

                notFloat = True
                notBig = True
                while notFloat and notBig:
                    notFloat = False
                    notBig = False
                    for w in weight:
                        if type(w) != float:
                            notFloat = True
                        if w<0:
                            notBig = True
                    if notFloat:
                        weight = input('Make sure you enter weight as an integer or float: ')
                    elif notBig:
                        weight = input('Weight must be positive: ')

                weights.append(weight)
                # set weights
                stg1 = True
                stg4 = False


                p.terminate()
            elif toggle == '2':
                # set to default values
                num_weights = len(boundaries[-1])
                weight = [1 for n in range(num_weights)]
                weights.append(weight)
                stg1 = True
                stg4 = False


                p.terminate()
            elif toggle == '3':
                boundaries.pop()  # remove selected boundaries
                # return to previous section
                stg4 = False
                stg3 = True

            elif toggle == '4':
                # set weights to default value and exit
                cont = False

    scanBounds = dict()
    for idx in range(len(boundaries)):
        scanBounds[int(scans[idx])] = (boundaries[idx], weights[idx])

    queue.put([scans, scanBounds])
    return

def createTable(data):
    """
    Purpose: create a table that shows the scan number and it's info in a new window
    :param data: array containing scan number, scan type, and scan name
    :return:
    """
    # Create a tkinter tree
    ws = Tk()  # create a new window
    ws.title('Scan info')
    ws.geometry('700x250')  # geometry of the new window
    ws['bg'] = '#AC99F2'  # colour of the window

    data_frame = Frame(ws)  # creates a new data frame
    data_frame.pack()  # packs the new data frame

    # initializes the scrollbar
    data_scroll = Scrollbar(data_frame)
    data_scroll.pack(side=RIGHT, fill=Y)

    # sets the scrollbar configurations
    my_data = ttk.Treeview(data_frame, yscrollcommand=data_scroll.set, xscrollcommand=data_scroll.set)
    my_data.pack()
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

    # insert the data to the data frame
    for idx in range(len(data)):
        my_data.insert(parent='', index='end', iid=idx, text='',
                       values=(data[idx][0], data[idx][1], data[idx][2]))

    my_data.pack()  # pack all of this into the data frame
    ws.mainloop()  # begin an infinite loop that runs until the window is closed

def changeSampleParams(x, parameters, sample):
    """
    Purpose: Change the parameters for the input sample for x-values
    :param x: A list that contains the new parameters values
    :param parameters: A list that defines which parameter to change in the sample definition
    :param sample: A slab object class
    :return: The new sample
    """
    # Loop through each sample parameter
    for p in range(len(parameters)):
        params = parameters[p]

        if params[0] == "SCALING FACTOR":
            sample.scaling_factor = x[p]
        elif params[0] == "BACKGROUND SHIFT":
            sample.background_shift = x[p]
        elif params[0] == 'SCATTERING FACTOR':
            mode =params[1]
            element = params[2]
            dE = x[p]

            if mode == 'STRUCTURAL':
                FfEnergyShift(element, dE, opt=True)
            elif mode == 'MAGNETIC':
                FfmEnergyShift(element, dE, opt=True)
        else:
            layer = params[0]  # Determine which layer to change
            property = params[1]  # structural/polymorphous/magnetic
            if property == 'STRUCTURAL':
                mode = params[2]
                # Determine if user wants to use compound or element mode
                if mode == 'COMPOUND':
                    characteristic = params[3]  # thickness/density/roughness/linked roughness
                    for ele in list(sample.structure[layer].keys()):
                        if characteristic == 'THICKNESS':
                            sample.structure[layer][ele].thickness = x[p]
                        elif characteristic == 'DENSITY':
                            stoich = sample.structure[layer][ele].stoichiometry  # stoichiometry
                            molar_mass = sample.structure[layer][ele].molar_mass  # molar mass
                            sample.structure[layer][ele].density = x[p]*stoich/molar_mass  # set density
                        elif characteristic == 'ROUGHNESS':
                            sample.structure[layer][ele].roughness = x[p]
                        elif characteristic == 'LINKED ROUGHNESS':
                            sample.structure[layer][ele].linked_roughness = x[p]

                elif mode == 'ELEMENT':
                    element = params[3]  # determines the element to change
                    characteristic = params[4]  # thickness/density/roughness/linked roughness
                    if characteristic == 'THICKNESS':
                        sample.structure[layer][element].thickness = x[p]
                    elif characteristic == 'DENSITY':
                        sample.structure[layer][element].density = x[p]
                    elif characteristic == 'ROUGHNESS':
                        sample.structure[layer][element].roughness = x[p]
                    elif characteristic == 'LINKED ROUGHNESS':
                        sample.structure[layer][element].linked_roughness = x[p]


            elif property == 'POLYMORPHOUS':
                element = params[2]  # determines the element that contains the polymorph
                polymorph = params[3]  # determines the polymorph to change

                ratio = 1 - x[p]  # Assumes only two possible polymorphs for now and computes other polymorph ratio

                poly = np.where(sample.structure[layer][element].polymorph == polymorph)  # determines location of polymorph

                # sets poly_ratio value making sure sum equals to 1
                if poly == 0:
                    sample.structure[layer][element].poly_ratio[0] = x[p]
                    sample.structure[layer][element].poly_ratio[1] = ratio
                elif poly == 1:
                    sample.structure[layer][element].poly_ratio[1] = x[p]
                    sample.structure[layer][element].poly_ratio[0] = ratio

            elif property == 'MAGNETIC':
                element = params[3]  # determines the magnetic element to use

                # determines if magnetic element is polymorphous
                if len(params) == 3:
                    sample.structure[layer][element].mag_density[0] = x[p]  # non-polymorphous case
                else:
                    polymorph = params[4]  # polymorphous case
                    poly = np.where(sample.structure[layer][element].polymorph == polymorph)
                    sample.structure[layer][element].mag_density[poly] = x[p]

    return sample

def scanCompute(x, *args):

    chi2 = 0  # what we will use to determine some values

    sample = args[0]  # slab class
    scans = args[1]  # data info
    data = args[2]  # data dict
    sims = args[3]  # simulation dict
    parameters = args[4]  # defines which parameters to change
    scanBounds = args[5]  # defines the bounds of the scans

    sample = changeSampleParams(x, parameters, sample)


    for scan in scans:
        scan_number = int(scan[0])
        scanType = scan[1]
        name = scan[2]
        scanbounds = scanBounds[scan_number]
        xbound = scanbounds[0]
        weights = scanbounds[1]

        if scanType == 'Reflectivity':
            myDataScan = data[name]
            myData = myDataScan['Data']
            E = myDataScan['Energy']
            pol = myDataScan['Polarization']
            qz = np.array(myData[0])
            Rdat = np.log(np.array(myData[2]))
            qz, Rsim = sample.reflectivity(E, qz)
            Rsim = np.log10(Rsim[pol])*sample.scaling_factor + np.ones(len(Rsim[pol]))*sample.background_shift


            for b in range(len(xbound)):
                lw = xbound[b][0]
                up = xbound[b][1]

                w = weights[b]

                idx = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] <= up]  #

                if len(idx) != 0:
                    chi2 = chi2 + sum((Rdat[idx]-Rsim[idx])**2/abs(Rsim[idx]))*w


            #chi2 = chi2 + sum((Rdat-Rsim)**2/abs(Rsim))

        elif scanType == 'Energy':
            myDataScan = data[name]
            myData = myDataScan['Data']
            Theta = myDataScan['Angle']
            E = np.array(myData[3])
            pol = myDataScan['Polarization']

            Rdat = np.log(np.array(myData[2]))
            Rsim = sample.energy_scan(Theta, E)
            Rsim = np.log(Rsim[pol])

            chi2 = chi2 + sum((Rdat-Rsim)**2/abs(Rsim))

    return chi2

def differential_evolution(fname,scan, parameters, bounds,scanBounds, strat = 'currenttobest1exp', mIter=25, tolerance=0.1, display=True):

    sample = ReadSampleHDF5(fname)  # import the sample information
    data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    # makes sure that scan is a list
    if type(scan) != list and type(scan) != np.ndarray:
        scan = [scan]

    scan = [s - 1 for s in scan]  # makes sure the indices are correct

    scans = data_info[scan]  # gets the appropriate scans

    params = [sample, scans, data, sims, parameters, scanBounds]  # required format for function scanCompute

    # This line will be used to select and use different global optimization algorithms
    ret = optimize.differential_evolution(scanCompute, bounds, args=params, strategy=strat,
                                          maxiter=mIter, tol=tolerance, disp=display)
    x = ret.x
    fun = ret.fun


    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    return x, fun

def shgo(fname, scan,parameters, bounds, scanBounds, N=64, iterations=3):

    sample = ReadSampleHDF5(fname)  # import the sample information
    f, data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    # makes sure that scan is a list
    if type(scan) != list or type(scan) != np.ndarray:
        scan = [scan]

    scan = [s - 1 for s in scan]  # makes sure the indices are correct
    scans = data_info[tuple(scan)]

    params = [sample, scans, data, sims, parameters, scanBounds]  # required format for function scanCompute

    ret = optimize.shgo(scanCompute, bounds, args=tuple(params), n=N, iters=iterations)
    x = ret.x
    fun = ret.fun

    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun

def dual_annealing(fname, scan, parameters, bounds,scanBounds, mIter=300):
    sample = ReadSampleHDF5(fname)  # import the sample information

    f, data_info, data, sims = ReadDataHDF5(fname)  # import the experimental data and simulated data

    # makes sure that scan is a list
    if type(scan) != list or type(scan) != np.ndarray:
        scan = [scan]

    scan = [s - 1 for s in scan]  # makes sure the indices are correct
    scans = data_info[tuple(scan)]

    params = [sample, scans, data, sims, parameters, scanBounds]  # required format for function scanCompute


    ret = optimize.dual_annealing(scanCompute, bounds, args=params, maxiter=mIter)
    x = ret.x
    fun = ret.fun


    print('Chi: ' + str(fun))
    print('Fitting parameters: ', x)

    f.close()
    return x, fun

def parameterSelection(sample, queue):
    """
    Purpose: Ask user to determine the parameters they would like to use
    :param sample: Sample information in slab format
    :param queue: Used for multiprocessing
    :return:
    """

    # Booleans that determine which process to complete

    parameters = list()  # list that keeps track of all the selected parameters
    constraints = list()
    param = list()  # keeps track of the current selection tree

    upperbound = list()  # upper bound of the parameter
    lowerbound = list()  # lower bound of the parameter
    compoundBounds = False  # compound mode bounds selection process
    elementBounds = False  # element mode bounds selection process
    polyBounds = False
    magneticBounds = False  # magnetic bounds selection process
    magBounds = False
    polyBounds = False  # polymorphous bounds selection process
    ffBounds = False  # form factor boundaries
    lastStep = False  # States that program has finished successfully
    val = 0  # keeps track of structure values for default boundaries

    print('PARAMETER SELECTION \n')

    cont = True  # continue parameter selection

    paramType = True   # parameter type

    layerSelect = False
    sampleParam = False  # structural
    bsSelect = False # background shift and scaling factor
    ffSelect = False # form factor energy shift
    constSelect = False  # constraint
    polyRatio = False
    constCompound = False
    constElement = False
    constThick = False
    isFinished = False

    constDict = dict()
    bSelect = False  # background shift
    sSelect = False  # scaling factor

    modeSelect = False  # structural mode selection
    structSelect = False  # parameters
    polySelect = False  # polymorphouse
    magSelect = False  # magnetic

    element_mode = False  # element mode
    compound_mode = False  # compound mode

    # dictionaries to keep track of which parameters have already been selected
    structFf = dict()  # structural  form factor parameters
    magFf = dict()  # magnetic form factor parameters
    polyFf = dict()  # polymorphous form factors parameters
    structDict = dict()  # stuctural parameters
    polyDict = dict()  # polymorphous parameters
    magDict = dict()
    bsTrack = ['Scaling Factor', 'Background Shift']  # keeps track if scaling factor or background shift have been selected

    # Retrieves dictionary info from sample input by user
    for i in range(len(sample.structure)):
        structDict[i] = dict()
        constDict[i] = []
        structDict[i]['compound'] = ['Thickness', 'Density', 'Roughness', 'Linked Roughness']
        structDict[i]['element'] = dict()
        if sample.layer_magnetized[i]:
            magDict[i] = dict()
        # structural scattering factors
        for ele in list(sample.structure[i].keys()):
            constDict[i].append(ele)
            if ele not in list(structFf.keys()):
                if ele not in list(sample.poly_elements.keys()) and ele not in list(sample.mag_elements.keys()):
                    structFf[ele] = sample.structure[i][ele].scattering_factor

            # polymorphous scattering factors
            polymorph = sample.structure[i][ele].polymorph
            if len(polymorph) > 0:
                polyDict[i] = []
                if ele not in list(polyDict.values()):
                    polyDict[i].append(ele)
            for j in range(len(polymorph)):
                if polymorph[j] not in list(polyFf.keys()):
                    polyFf[polymorph[j]] = sample.structure[i][ele].scattering_factor[j]

            # magnetic scattering factors
            mag = sample.structure[i][ele].mag_density
            if len(mag) > 0:

                if len(polymorph) > 0:
                    magDict[i][ele] = polymorph
                    for j in range(len(polymorph)):
                        if polymorph[j] not in list(magFf.keys()):
                            magFf[polymorph[j]] = sample.structure[i][ele].mag_scattering_factor[j]
                else:
                    magDict[i][ele] = ele
                    if ele not in list(magFf.keys()):
                        magFf[ele] = sample.structure[i][ele].mag_scattering_factor

            structDict[i]['element'][ele] = ['Thickness', 'Density', 'Roughness', 'Linked Roughness']



    # start of main process
    while cont:
        # ----------------------------------------------------------------------------------------- #
        # Selecting which kind of parameter user wants to optimize -------------------------------- #
        # ----------------------------------------------------------------------------------------- #
        if paramType:
            param = list()
            print('Select which parameter type to optimize: ')
            print('\t 1: Background Shift/Scaling Factor')
            print('\t 2: Sample Parameters')
            print('\t 3: Form Factor Energy Shift')
            print('\t 4: Constraints')
            print('\t 5: Exit/Finish')

            toggle = input('\n -> ')
            print()

            while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4' and toggle != '5':
                toggle = input('Please input the number of your selection: ')

            if toggle == '1':
                # Background shift and Scaling Factor
                paramType = False
                bsSelect = True
            elif toggle == '2':
                # Sample Parameters
                paramType = False
                layerSelect = True
            elif toggle == '3':
                # Scattering Factor Energy Shift
                paramType = False
                ffSelect = True
            elif toggle == '4':
                # Constraints
                paramType = False
                constSelect = True

            elif toggle == '5':
                # exit the program
                cont = False
                isFinished = True

        # -------------------------------------------------------------------------------------- #
        # Selecting whether to vary background shift or scaling factor ------------------------- #
        # -------------------------------------------------------------------------------------- #
        elif bsSelect:
            temp = dict()
            print('BACKGROUND SHIFT AND SCALING FACTOR \n')
            print('Choose an option: ')
            idx = 1
            for var in bsTrack:
                print('\t '+str(idx)+': ' + var)
                temp[str(idx)] = var
                idx = idx + 1
            if len(bsTrack) == 2:
                print('\t ' + str(idx) + ': Background Shift and Scaling Factor')
                temp[str(idx)] = 'Background Shift and Scaling Factor'
                idx = idx + 1
            print('\t '+str(idx)+': Return')
            temp[str(idx)] = 'Return'
            idx = idx + 1
            print('\t '+str(idx)+': Exit')
            temp[str(idx)] = 'Exit'
            toggle = input('\n ->')
            print()

            while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4' and toggle != '5':
                toggle = input('Please input the number of your selection: ')

            background = False
            scaling = False
            backScale = False
            if temp[toggle] == 'Background Shift':
                # background shift selected
                background = True
            elif temp[toggle] == 'Scaling Factor':
                 # scaling factor selected
                scaling = True
            elif temp[toggle] =='Background Shift and Scaling Factor':
                # both background shift and scaling factor selected
                backScale = True
            elif temp[toggle] == 'Return':
                # return to parameter selection
                bsSelect = False
                paramType = True
            elif temp[toggle] == 'Exit':
                # exit the program where function not completed correctly
                cont = False

            if cont and not(paramType):
                print(temp[toggle].upper() + ' SELECTION \n')
                print('Select an option: ')
                print('\t 1: Select bounds')
                print('\t 2: Use default bounds')
                print('\t 3: Return')
                print('\t 4: Exit')
                toggle = input('\n -> ')
                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                    toggle = input('Select one of the provided options: ')
                print()

                if toggle == '3':
                    background = False
                    scaling = False
                    backScale = False
                elif toggle == '4':
                    cont = False
                else:
                    if background:
                        bsTrack.remove('Background Shift')
                        param.append('BACKGROUND SHIFT')
                        if toggle == '1':
                            bBounds = input('Select the background shift optimization bounds (lower, upper): ')
                            bBounds = ast.literal_eval(bBounds)
                            while type(bBounds) != tuple or float(bBounds[0]) > float(bBounds[1]):
                                bBounds = input('Input bounds as a tuple in ascending order (lower, upper): ')
                                bBounds = ast.literal_eval(bBounds)

                            lowerbound.append(float(bBounds[0]))
                            upperbound.append(float(bBounds[1]))
                            print()

                        elif toggle == '2':
                            lowerbound.append(-5e-7)
                            upperbound.append(5e-7)

                        parameters.append(param.copy())
                    elif scaling:
                        bsTrack.remove('Scaling Factor')
                        param.append('SCALING FACTOR')
                        if toggle == '1':
                            sBounds = input('Select the scaling factor optimization bounds (lower, upper): ')
                            sBounds = ast.literal_eval(sBounds)
                            while type(sBounds) != tuple or float(sBounds[0]) > float(sBounds[1]):
                                sBounds = input('Input bounds as a tuple in ascending order (lower, upper): ')
                                sBounds = ast.literal_eval(sBounds)

                            lowerbound.append(float(sBounds[0]))
                            upperbound.append(float(sBounds[1]))
                            print()

                        elif toggle == '2':
                            lowerbound.append(0.8)
                            upperbound.append(1.2)

                        parameters.append(param.copy())
                    elif backScale:
                        bsTrack.remove('Background Shift')
                        param.append('BACKGROUND SHIFT')
                        if toggle == '1':
                            bBounds = input('Select the background shift optimization bounds (lower, upper): ')
                            bBounds = ast.literal_eval(bBounds)
                            while type(bBounds) != tuple or float(bBounds[0]) > float(bBounds[1]):
                                bBounds = input('Input bounds as a tuple in ascending order (lower, upper): ')
                                bBounds = ast.literal_eval(bBounds)

                            lowerbound.append(float(bBounds[0]))
                            upperbound.append(float(bBounds[1]))

                        elif toggle == '2':
                            lowerbound.append(-5e-7)
                            upperbound.append(5e-7)

                        parameters.append(param.copy())
                        param = list()
                        print(parameters)
                        bsTrack.remove('Scaling Factor')
                        param.append('SCALING FACTOR')
                        if toggle == '1':
                            sBounds = input('Select the scaling factor optimization bounds (lower, upper): ')
                            sBounds = ast.literal_eval(sBounds)
                            while type(sBounds) != tuple or float(sBounds[0]) > float(sBounds[1]):
                                sBounds = input('Input bounds as a tuple in ascending order (lower, upper): ')
                                sBounds = ast.literal_eval(sBounds)

                            lowerbound.append(float(sBounds[0]))
                            upperbound.append(float(sBounds[1]))
                            print()

                        elif toggle == '2':
                            lowerbound.append(0.8)
                            upperbound.append(1.2)

                        parameters.append(param.copy())

                    print('SCALING FACTOR AND BACKGROUND SHIFT FINISH \n')
                    print('Select an option: ')
                    print('\t 1: Select another parameter type')
                    print('\t 2: Return')
                    print('\t 3: Exit/Finish')
                    toggle = input('\n -> ')
                    print()

                    while toggle != '1' and toggle != '2' and toggle != '3':
                        toggle = input('Select one of the provided options: ')

                    if toggle == '1':
                        background = False
                        scaling = False
                        backScale = False
                        bsSelect = False
                        paramType = True
                    elif toggle == '2':
                        if background or scaling:
                            if background:
                                bsTrack.append('Background Shift')
                            elif scaling:
                                bsTrack.append('Scaling Factor')

                            upperbound.pop()
                            lowerbound.pop()
                            parameters.pop()
                            param = list()
                        elif backScale:
                            bsTrack = ['Background Shift', 'Scaling Factor']
                            upperbound.pop()
                            upperbound.pop()
                            lowerbound.pop()
                            lowerbound.pop()
                            parameters.pop()
                            parameters.pop()
                            param = list()
                        background = False
                        scaling = False
                        backScale = False
                        isFinished = True
                    elif toggle == '3':
                        cont = False
        # --------------------------------------------------------------------------------------------------------- #
        # Form factor selection
        # --------------------------------------------------------------------------------------------------------- #
        elif ffSelect:

            print('FORM FACTOR ENERGY SHIFT SELECTION \n')
            print('Select an option: ')
            print('\t 1: Structural form factor')
            print('\t 2: Polymorphous form factor')
            print('\t 3: Magnetic form factor')
            print('\t 4: Return')
            print('\t 5: Exit')
            toggle = input('\n -> ')
            print()

            selected_ele = dict()
            while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4' and toggle != '5':
                toggle = input('Select one of the provided options')

            if toggle == '1':

                print('STRUCTURAL SCATTERING FACTOR \n')
                print('Select which element scattering factor to shift: ')
                temp = dict()
                idx = 1
                for key in list(structFf.keys()):
                    print('\t ' + str(idx) + ': ' + key)
                    temp[str(idx)] = key
                    idx = idx + 1
                print('\t ' + str(idx) + ': Return')
                temp[str(idx)] = 'Return'
                idx = idx + 1
                print('\t ' + str(idx) + ': Exit')
                temp[str(idx)] = 'Exit'
                toggle = input('\n -> ')
                print()

                while toggle not in list(temp.keys()):
                    toggle = input('Select one of the provided options: ')
                print()

                if temp[toggle] in list(structFf.keys()):
                    param.append('SCATTERING FACTOR')
                    param.append('STRUCTURAL')
                    param.append(structFf[temp[toggle]])
                    selected_ele['STRUCTURAL'] = temp[toggle]
                    ffBounds = True
                elif temp[toggle] == 'Return':
                    pass
                elif toggle == 'Exit':
                    cont = False
            elif toggle == '2':
                print('POLYMORPHOUS SCATTERING FACTOR \n')
                print('Select which element scattering factor to shift: ')
                temp = dict()
                idx = 1
                for key in list(polyFf.keys()):
                    print('\t ' + str(idx) + ': ' + key)
                    temp[str(idx)] = key
                    idx = idx + 1
                print('\t ' + str(idx) + ': Return')
                temp[str(idx)] = 'Return'
                idx = idx + 1
                print('\t ' + str(idx) + ': Exit')
                temp[str(idx)] = 'Exit'
                toggle = input('\n -> ')
                print()

                while toggle not in list(temp.keys()):
                    toggle = input('Select one of the provided options: ')
                print()

                if temp[toggle] in list(polyFf.keys()):
                    param.append('SCATTERING FACTOR')
                    param.append('STRUCTURAL')
                    param.append(polyFf[temp[toggle]])
                    selected_ele['POLYMORPHOUS'] = temp[toggle]

                    ffBounds = True
                    ffSelect=False
                elif temp[toggle] == 'Return':
                    pass
                elif toggle == 'Exit':
                    cont = False
            elif toggle == '3':

                print('MAGNETIC SCATTERING FACTOR \n')
                print('Select which element magnetic scattering factor to shift: ')
                temp = dict()
                idx = 1
                for key in list(magFf.keys()):
                    print('\t ' + str(idx) + ': ' + key)
                    temp[str(idx)] = key
                    idx = idx + 1
                print('\t ' + str(idx) + ': Return')
                temp[str(idx)] = 'Return'
                idx = idx + 1
                print('\t ' + str(idx) + ': Exit')
                temp[str(idx)] = 'Exit'
                toggle = input('\n -> ')
                print()

                while toggle not in list(temp.keys()):
                    toggle = input('Select one of the provided options: ')
                print()

                if temp[toggle] in list(magFf.keys()):
                    param.append('SCATTERING FACTOR')
                    param.append('MAGNETIC')
                    param.append(magFf[temp[toggle]])
                    selected_ele['MAGNETIC'] = temp[toggle]
                    ffBounds = True
                elif temp[toggle] == 'Return':
                    pass
                elif toggle == 'Exit':
                    cont = False
            elif toggle == '4':
                ffSelect = False
                paramType = True
                param = []
            elif toggle == '5':
                cont = False

            if ffBounds:
                print('FORM FACTOR BOUNDS \n')
                print('Select an option: ')
                print('\t 1: Select energy bounds')
                print('\t 2: Use default energy bounds')
                print('\t 3: Return')
                print('\t 4: Exit')
                toggle = input('\n -> ')
                print()

                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                    toggle = input('Select one of the provided options:')
                print()

                if toggle == '1':
                    bound = input('Select the energy shift bounds in eV as a tuple(lower, upper): ')
                    bound = ast.literal_eval(bound)
                    bound = (float(bound[0]), float(bound[1]))
                    while bound[0] > bound[1]:
                        bound = input('Make sure the lowerbound is smaller than the upperbound (lower, upper): ')
                        bound = ast.literal_eval(bound)
                        bound = (float(bound[0]), float(bound[1]))

                    upperbound.append(bound[0])
                    lowerbound.append(bound[1])

                elif toggle == '2':
                    lowerbound.append(-0.5)
                    upperbound.append(0.5)
                elif toggle == '3':
                    ffbounds = False
                    selected_ele = dict()
                    param = []
                elif toggle == '4':
                    cont = False

                if toggle == '1' or toggle == '2':
                    print('FINISH SCATTERING FACTOR \n')
                    print('Select an option: ')
                    print('\t 1: Select another scattering factor')
                    print('\t 2: Select another parameter')
                    print('\t 3: Return')
                    print('\t 4: Exit/Finish')
                    toggle = input('\n -> ')
                    print()
                    while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                        toggle = input('Select one of the provided options')
                    print()

                    if toggle == '1':
                        parameters.append(param.copy())
                        ffBounds = False
                        ffSelect = True
                        param = []
                        key = list(selected_ele.keys())[0]
                        if key == 'STRUCTURAL':
                            del structFf[selected_ele[key]]
                        elif key == 'POLYMORPHOUS':
                            del polyFf[selected_ele[key]]
                        elif key == 'MAGNETIC':
                            del magFf[selected_ele[key]]

                        selected_ele = dict()
                    elif toggle == '2':
                        parameters.append(param.copy())
                        ffBounds = False
                        ffSelect = False
                        paramType = True
                        param = []
                        key = list(selected_ele.keys())[0]
                        if key == 'STRUCTURAL':
                            del structFf[selected_ele[key]]
                        elif key == 'POLYMORPHOUS':
                            del polyFf[selected_ele[key]]
                        elif key == 'MAGNETIC':
                            del magFf[selected_ele[key]]

                        selected_ele = dict()
                    elif toggle == '3':
                        ffBounds = False
                        ffSelect = True
                        upperbound.pop()
                        lowerbound.pop()
                        selected_ele = dict()
                        param = []
                    elif toggle == '4':
                        parameters.append(param.copy())
                        cont = False
                        isFinished = True


        elif layerSelect:
            num_layers = len(sample.structure)
            num_layers_list = [i for i in range(num_layers)]
            print('LAYER SELECTION \n')
            toggle = input('Select which layer you want to optimize (0-'+str(num_layers-1)+'): ')
            while int(toggle) not in num_layers_list:
                toggle = print('Enter the layer number you want to select (0-' + str(num_layers - 1) + '): ')

            print()
            layerSelect = False
            sampleParam = True
            param.append(int(toggle))
        # Select sample params ------------------------------------------------------------------
        elif sampleParam:
            print(param)
            layer = sample.structure[param[0]]
            poly_elements = list(sample.poly_elements.keys())
            poly_exists = False

            for ele in list(layer.keys()):
                if ele in poly_elements:
                    poly_exists = True

            temp = dict() # keeps track of which selection is made
            idx = 1
            print('SAMPLE PARAMETERS \n')
            print('Choose an option: ')
            print('\t '+ str(idx) + ': Structural')
            temp[str(idx)] = 'Structural'
            idx = idx + 1
            if poly_exists:
                print('\t '+ str(idx) +': Polymorphous')
                temp[str(idx)] = 'Polymorphous'
                idx = idx + 1
            if sample.layer_magnetized[param[0]]:
                print('\t '+ str(idx) +': Magnetic')
                temp[str(idx)] = 'Magnetic'
                idx = idx + 1
            print('\t '+ str(idx) +': Return to layer selection')
            temp[str(idx)] = 'Return to layer selection'
            idx = idx + 1
            print('\t ' + str(idx) + ': Return to parameter selection')
            temp[str(idx)] = 'Return to parameter selection'
            idx = idx + 1
            print('\t '+ str(idx) +': Exit')
            temp[str(idx)] = 'Exit'
            toggle = input('\n -> ')
            print()

            while toggle not in list(temp.keys()):
                toggle = input('Please input the number of your selection: ')


            if temp[toggle] == 'Structural':    # structural

                modeSelect = True
                sampleParam = False
                param.append('STRUCTURAL')
            elif temp[toggle] == 'Polymorphous':  # polymorphous
                polySelect = True
                sampleParam = False
                param.append('POLYMORPHOUS')
            elif temp[toggle] == 'Magnetic':  # magnetic
                magSelect = True
                sampleParam = False
                param.append('MAGNETIC')
            elif temp[toggle] == 'Return to layer selection':  # return
                layerSelect = True
                sampleParam = False
                param.pop()
            elif temp[toggle] == 'Return to parameter selection':  # return
                paramType = True
                sampleParam = False
                param.pop()
            elif temp[toggle] == 'Exit':  # exit
                cont = False


        elif modeSelect:  # Select Mode --------------------------------------------------------------------------------
            print('STRUCTURAL MODE SELECTION \n')


            print('Select an option:')
            print('\t 1: Compound')
            print('\t 2: Element')
            print('\t 3: Return')
            print('\t 4: Exit')
            toggle = input('\n -> ')
            print()
            while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                toggle = input("Select one of the provided options: ")
            print()
            if toggle == '1':
                compound_mode = True
                modeSelect = False
                param.append('COMPOUND')
            elif toggle == '2':
                element_mode = True
                modeSelect = False
                param.append('ELEMENT')
            elif toggle == '3':
                modeSelect = False
                sampleParam = True
            elif toggle == '4':

                cont = False

        # Polymorphouse case ------------------------------------------------------------------------------------------
        elif polySelect:

            print('POLYMORPHOUS SELECTION \n')

            print('Select which polymorph you would like to optimize: ')
            idx = 1
            temp = dict()
            for ele in polyDict[param[0]]:
                print('\t ' + str(idx) + ': ' + ele)
                temp[str(idx)] = ele
                idx = idx + 1
            print('\t ' + str(idx) + ': Return')
            temp[str(idx)] = 'Return'
            idx = idx + 1
            temp[str(idx)] = 'Exit'
            print('\t ' + str(idx) + ': Exit')

            toggle = input('\n -> ')
            print()

            while toggle not in list(temp.keys()):
                toggle = input('Select one of the provided options: ')
            print()

            if temp[toggle] == 'Return':
                param.pop()
                polySelect = False
                sampleParam = True
            elif temp[toggle] == 'Exit':
                cont = False
            else:
                param.append(temp[toggle])
                polyRatio = True
                polySelect = False

        elif polyRatio:
            polymorphs = sample.structure[param[0]][param[-1]].polymorph
            poly_ratio = sample.structure[param[0]][param[-1]].poly_ratio
            print('SELECT POLYMORPH RATIO \n')
            if len(polymorphs) > 3:
                raise RuntimeError('This code cannot handle more than three polymorphs.')

            print('Select polymorph that you would like to control: ')
            idx = 1
            temp = dict()
            for poly in polymorphs:
                print('\t ' + str(idx)+ ': ' + poly)
                temp[str(idx)] = poly
                idx = idx + 1
            print('\t ' + str(idx) + ': Return')
            temp[str(idx)] = 'Return'
            idx = idx + 1
            print('\t ' + str(idx) + ': Exit')
            temp[str(idx)] = 'Exit'

            toggle = input('\n -> ')
            print()

            while toggle not in list(temp.keys()):
                toggle = input('Select one of the provided options: ')

            if temp[toggle] == 'Return':
                param.pop()
                polyRatio = False
                polySelect = True
            elif temp[toggle] == 'Exit':
                cont = False
            else:
                selected_poly = temp[toggle]
                param.append(temp[toggle])
                if len(polymorphs) == 3:
                    polymorphs.remove(selected_poly)
                    ratio = input('Enter the ratio relation for ' + polymorphs[0] + ' and ' + polymorphs[1] + ': ')
                    ratio = ratio.split()
                    good = True
                    ratio = []
                    while len(ratio) != 2 and not(good):
                        if good and len(ratio) != 2:
                            ratio = input('Separate the numbers with a space: ')
                        elif not(good):
                            ratio = input('Ratio must be an integer or float type: ')
                        if len(ratio) == 2:
                            ratio = [float(ratio[0]), float(ratio[1])]
                            if type(ratio[0]) == float and type(ratio[1]) == float:
                                good = True
                    print()
                    param.append(polymorphs)
                    param.append(ratio)

                polyBounds = True
                polyRatio = False

        elif polyBounds:
            my_ele = param[2]
            my_poly = param[3]
            print('POLYMORPH BOUNDARIES \n')
            print('Select an option: ')
            print('\t 1: Select boundaries')
            print('\t 2: Use default boundaries')
            print('\t 3: Return')
            print('\t 4: Exit')

            toggle = input('\n -> ')
            while toggle not in ['1','2','3','4']:
                toggle = input('Select one of the provided options: ')
            print()

            finishPoly = False
            if toggle == '1':
                bound = input('Input the polmorph ratio bound between 0 and 1 for ' + my_poly + ' (lower, upper): ')
                bound = ast.literal_eval(bound)
                bound = (float(bound[0]), float(bound[1]))
                while bound[0] > bound[1] and bound[0]<0 and bound[1]>1:
                    bound = input('Make sure upper and lower bounds are input in ascending order and are found between 0 and 1: ')
                    bound = ast.literal_eval(bound)
                    bound = (float(bound[0]), float(bound[1]))

                lowerbound.append(bound[0])
                upperbound.append(bound[1])
                finishPoly = True
                polyBounds = False
            elif toggle == '2':
                # find index
                poly = sample.structure[param[0]][my_ele].polymorph
                if type(poly) == np.ndarray:
                    idx = np.where(poly == my_poly)
                    idx = idx[0][0]
                elif type(poly) == list:
                    idx = poly.index(my_poly)

                var = sample.structure[param[0]][my_ele].poly_ratio[idx]

                lw = var - 0.2
                up = var + 0.2

                if lw < 0:
                    lw = 0
                if up > 1:
                    up = 1

                upperbound.append(up)
                lowerbound.append(lw)
                finishPoly = True
                polyBounds = False
            elif toggle == '3':
                polyRatio = True
                if len(param) == 4:
                    param.pop()
                    #param.pop()
                else:
                    #param.pop()
                    param.pop()
                    param.pop()
                    param.pop()

            elif toggle == '4':
               cont = False

            if finishPoly:
                print('FINISH POLYMORPH \n')
                print('Make a selection:')
                print('\t 1: Select another parameter for current layer')
                print('\t 2: Select another layer')
                print('\t 3: Select another parameter type')
                print('\t 4: Return')
                print('\t 5: Exit/Finish')

                toggle = input('\n -> ')
                print()
                while toggle not in ['1','2','3','4','5']:
                    toggle = input('Select one of the provided options: ')
                print()
                if toggle == '1':
                    parameters.append(param.copy())
                    finishPoly = False
                    sampleParam = True
                    polyDict[param[0]].remove(my_ele)
                    param = [param[0]]

                elif toggle =='2':
                    parameters.append(param.copy())

                    finishPoly = False
                    layerSelect = True
                    polyDict[param[0]].remove(my_ele)
                    param = []

                elif toggle =='3':
                    parameters.append(param.copy())
                    finishPoly = False
                    paramType = True
                    polyDict[param[0]].remove(my_ele)
                    param = []
                elif toggle =='4':
                    finishPoly = False
                    polyBounds = True
                    lowerbound.pop()
                    upperbound.pop()
                elif toggle =='5':
                    cont = False
                    isFinished = True
                    parameters.append(param.copy())

        elif magSelect:
            print('MAGNETIC SELECTION \n')
            print('Select which magnetic element you would like to vary: ')
            temp = dict()
            selected_mag = ''
            selected_mag_poly = ''
            idx = 1
            for ele in list(magDict[param[0]].keys()):
                print('\t ' + str(idx) + ': ' + ele)
                temp[str(idx)] = ele
                idx = idx + 1
            print('\t ' + str(idx) + ': Return')
            temp[str(idx)] = 'Return'
            idx = idx + 1
            print('\t ' + str(idx) + ': Exit')
            temp[str(idx)] = 'Exit'
            toggle = input('\n -> ')
            print()

            isPoly = False
            while toggle not in list(temp.keys()):
                toggle = input('Select one of the provided options: ')
            print()

            if temp[toggle] == 'Return':
                magSelect = False
                sampleParam = True
            elif toggle == 'Exit':
                cont = False
            else:
                param.append(temp[toggle])
                magSelect = False
                selected_mag = temp[toggle]
                if len(sample.structure[param[0]][selected_mag].polymorph) > 0:
                    isPoly = True

            if isPoly:
                print('SELECT MAGNETIC POLYMORPH \n')
                idx = 1
                temp = dict()
                for poly in magDict[param[0]][selected_mag]:
                    print('\t ' + str(idx) + ': ' + poly)
                    temp[str(idx)] = poly
                    idx = idx + 1
                print('\t ' + str(idx) + ': Return')
                temp[str(idx)] = 'Return'
                idx = idx + 1
                print('\t ' + str(idx) + ': Exit')
                temp[str(idx)] = 'Exit'

                toggle = input('\n -> ')
                print()
                while toggle not in list(temp.keys()):
                    toggle = input('Select one of the provided options: ')
                print()

                if temp[toggle] == 'Return':
                    param.pop()
                    isPoly = False
                    magSelect = True

                elif temp[toggle] == 'Exit':
                    cont = False
                else:
                    selected_mag_poly = temp[toggle]
                    param.append(selected_mag_poly)
                    magBounds = True
                    selectMag = False
                    isPoly = False

        elif magBounds:
            print('MAGNETIC BOUNDS \n')
            print('Select one of the boundary options: ')
            print('\t 1: Select boundaries')
            print('\t 2: Default boundaries')
            print('\t 3: Return')
            print('\t 4: Exit')
            toggle = input('\n -> ')
            print()
            while toggle not in ['1','2','3','4']:
                toggle = input('Select one of the provided options: ')
            print()

            contMag = False
            if toggle == '1':
                bd = input('Enter the magnetic density boundary in terms of mol/cm^3 (lower, upper): ')
                bd = ast.literal_eval(bd)
                bd = (float(bd[0]), float(bd[1]))
                while bd[0] > bd[1] and bd[0] < 0 and type(bd) != tuple:
                    bd = input('Enter the parameter optimization boundary as a tuple (lower, upper) in ascending order: ')
                    bd = ast.literal_eval(bd)
                    bd = (float(bd[0]), float(bd[1]))

                lowerbound.append(bd[0])
                upperbound.append(bd[1])
                contMag = True
            elif toggle == '2':
                lw = 0
                up = 0
                if len(param) == 2:
                    # non polymorphous case
                    val = sample.structure[param[0]][param[2]].mag_density
                    if type(val) == list or type(val) == np.ndarray:
                        val = val[0]
                    lw = val - 1e-5
                    up = val + 1e-5
                else:
                    # polymorphous case
                    print(param)
                    var = sample.structure[param[0]][param[2]].mag_density
                    polymorph = sample.structure[param[0]][param[2]].polymorph
                    if type(polymorph) == list:
                        idx = polymorph.index(param[3])
                        idx = idx[0]
                    elif type(polymorph) == np.ndarray:
                        idx = np.where(param[3] == polymorph)
                        idx = idx[0][0]
                    val = var[idx]
                    up = val + 1e-5
                    lw = val - 1e-5
                upperbound.append(up)
                lowerbound.append(lw)
                contMag = True
            elif toggle == '3':
                magSelect = True
                magBounds = False
                if len(param) == 2:
                    param.pop()
                else:
                    param.pop()
                    param.pop()
            elif toggle == '4':
                cont = False

            if contMag:
                print('FINISH MAGNETIC DENSITY \n')
                print('Select an option: ')
                print('\t 1: Select another magnetic density for the same layer')
                print('\t 2: Select another parameter for the same layer')
                print('\t 3: Select another layer')
                print('\t 4: Select another parameter type')
                print('\t 5: Return')
                print('\t 6: Exit/Finish')
                toggle = input('\n -> ')
                print()
                while toggle not in ['1','2','3','4','5','6']:
                    toggle = input('Select one of the provided options: ')
                if toggle in ['1','2','3']:
                    parameters.append(param.copy())
                    magBounds = False
                    magLayer = magDict[param[0]][param[2]]
                    if len(magLayer) > 1:
                        if type(magLayer) == np.ndarray:
                            print(magDict[param[0]][param[2]])
                            magDict[param[0]][param[2]] = np.delete(magDict[param[0]][param[2]], np.where(magDict[param[0]][param[2]] == param[3])[0])
                        else:
                            magDict[param[0]][param[2]].remove(param[3])
                    else:
                        del magDict[param[0]][param[2]]
                    print(magDict)
                    if toggle == '1':
                        param = param[0:2]
                        magSelect = True
                    elif toggle == '2':
                        param = [param[0]]
                        sampleParam = True
                    elif toggle == '3':
                        param = list()
                        paramType = True
                    print(magDict)
                elif toggle == '4':
                    print()
                elif toggle == '5':
                    lowerbound.pop()
                    upperbound.pop()
                elif toggle == '6':
                    cont = False
                    isFinished = True

        # Compound Mode -----------------------------------------------------------------------------------------------
        elif compound_mode:
            temp = dict()
            print('COMPOUND MODE \n')
            print('Select an option: ')
            idx = 1
            for char in structDict[param[0]]['compound']:
                print('\t '+str(idx)+': ' + char)
                temp[str(idx)] = char
                idx = idx + 1

            print('\t '+str(idx)+': Return')
            temp[str(idx)] = 'Return'
            idx = idx + 1
            print('\t ' + str(idx) + ': Exit')
            temp[str(idx)] = 'Exit'

            toggle = input('\n -> ')
            print()
            while toggle not in list(temp.keys()):
                toggle = input('Select one of the options provided: ')
            print()

            if temp[toggle] == 'Thickness':
                if 'Thickness' in structDict[param[0]]['compound']:
                    structDict[param[0]]['compound'].remove('Thickness')
                for ele in list(structDict[param[0]]['element'].keys()):
                    structDict[param[0]]['element'][ele].remove('Thickness')
                    val = sample.structure[param[0]][ele].thickness
                param.append('THICKNESS')

                compound_mode = False
                compoundBounds = True

            elif temp[toggle] == 'Density':
                if 'Density' in structDict[param[0]]['compound']:
                    structDict[param[0]]['compound'].remove('Density')
                for ele in list(structDict[param[0]]['element'].keys()):
                    structDict[param[0]]['element'][ele].remove('Density')
                    val = sample.structure[param[0]][ele].density
                param.append('DENSITY')

                compound_mode = False
                compoundBounds = True

            elif temp[toggle] == 'Roughness':
                if 'Roughness' in structDict[param[0]]['compound']:
                    structDict[param[0]]['compound'].remove('Roughness')
                for ele in list(structDict[param[0]]['element'].keys()):
                    structDict[param[0]]['element'][ele].remove('Roughness')
                    val = sample.structure[param[0]][ele].roughness
                param.append('ROUGHNESS')

                compound_mode = False
                compoundBounds = True

            elif temp[toggle] == 'Linked Roughness':
                if 'Linked Roughness' in structDict[param[0]]['compound']:
                    structDict[param[0]]['compound'].remove('Linked Roughness')
                for ele in list(structDict[param[0]]['element'].keys()):
                    structDict[param[0]]['element'][ele].remove('Linked Roughness')
                    val = sample.structure[param[0]][ele].linked_roughness

                param.append('LINKED ROUGHNESS')

                compound_mode = False
                compoundBounds = True

            elif temp[toggle] == 'Return':
                compound_mode = False
                modeSelect = True
                param.pop()  # remove mode selection
            elif temp[toggle] == 'Exit':
                cont = False

        elif element_mode:
            print('ELEMENT MODE \n')
            print('Select element: ')
            idx = 1
            elements = list(sample.structure[param[0]].keys())
            for ele in elements:
                print('\t ' + str(idx) + ': ' + ele)
                idx = idx + 1
            print('\t ' + str(idx) + ': Return')
            idx = idx+1
            print('\t ' + str(idx) + ': Exit')
            toggle = input('\n -> ')
            print()

            while toggle != '1' and toggle !='2' and toggle != '3' and toggle != '4' and toggle !='5':
                toggle = input('Select one of elements provided: ')
            print()

            ele = ''
            if toggle == '1':
                ele = elements[0]
                param.append(ele)
            elif toggle =='2':
                ele = elements[1]
                param.append(ele)
            elif toggle == '3':
                ele = elements[2]
                param.append(ele)
            elif toggle == '4':
                param.pop()
                element_mode = False
                modeSelect = True
            elif toggle =='5':
                cont = False

            if toggle == '1' or toggle =='2' or toggle =='3':
                print('ELEMENT PARAMETER SELECTION \n')
                print('Select parameter you want to vary for ' + ele + ': ')
                idx = 1
                temp = dict()
                for char in structDict[param[0]]['element'][ele]:
                    print('\t ' + str(idx) + ': '+char)
                    temp[str(idx)] = char
                    idx = idx + 1
                print('\t ' + str(idx) +': Return to element selection')
                temp[str(idx)] = 'Return to element selection'
                idx = idx + 1
                print('\t ' + str(idx) + ': Return to mode selection')
                temp[str(idx)] = 'Return to mode selection'
                idx = idx + 1
                print('\t ' + str(idx) + ': Exit')
                temp[str(idx)] = 'Exit'

                toggle = input("\n -> ")
                print()
                while toggle not in list(temp.keys()):
                    toggle = input('Select one of the provided options: ')

                print()
                if temp[toggle] == 'Thickness':
                    structDict[param[0]]['element'][ele].remove('Thickness')
                    if 'Thickness' in structDict[param[0]]['compound']:
                        structDict[param[0]]['compound'].remove('Thickness')
                    param.append('THICKNESS')
                    elementBounds = True
                    element_mode = False
                elif temp[toggle] == 'Density':
                    structDict[param[0]]['element'][ele].remove('Density')
                    if 'Density' in structDict[param[0]]['compound']:
                        structDict[param[0]]['compound'].remove('Density')
                    param.append('DENSITY')
                    elementBounds = True
                    element_mode = False
                elif temp[toggle] == 'Roughness':
                    structDict[param[0]]['element'][ele].remove('Roughness')
                    if 'Roughness' in structDict[param[0]]['compound']:
                        structDict[param[0]]['compound'].remove('Roughness')
                    param.append('ROUGHNESS')
                    elementBounds = True
                    element_mode = False
                elif temp[toggle] == 'Linked Roughness':
                    structDict[param[0]]['element'][ele].remove('Linked Roughness')
                    if 'Linked Roughness' in structDict[param[0]]['compound']:
                        structDict[param[0]]['compound'].remove('Linked Roughness')
                    param.append('LINKED ROUGHNESS')
                    elementBounds = True
                    element_mode = False
                elif temp[toggle] == 'Return to element selection':
                    param.pop()
                elif temp[toggle] == 'Return to mode selection':
                    element_mode = False
                    modeSelect = True
                    param.pop()
                    param.pop()

                elif temp[toggle] == 'Exit':
                    cont = False



        elif compoundBounds:  # compound bounds ---------------------------------------------------------
            print('COMPOUND BOUND SELECTION \n')
            print('Select an option: ')
            print('\t 1: Select parameter boundaries')
            print('\t 2: Use default boundaries')
            print('\t 3: Return')
            print('\t 4: Exit')
            toggle = input('\n -> ')
            print()
            while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                toggle = input('Select one of the provided options:')

            print()
            if toggle == '1':
                bd = input('Enter the parameter optimization boundary as a tuple (lower, upper): ')
                bd = ast.literal_eval(bd)
                bd = (float(bd[0]), float(bd[1]))
                while bd[0]>bd[1] and bd[0] < 0 and type(bd) != tuple:
                    bd = input('Enter the parameter optimization boundary as a tuple (lower, upper) in ascending order: ')
                bd = ast.literal_eval(bd)
                bd = (float(bd[0]), float(bd[1]))
                lowerbound.append(bd[0])
                upperbound.append(bd[1])

            elif toggle == '2':
                characteristic = param[-1]
                if characteristic == 'THICKNESS':
                    lw = val - 5
                    up = val + 5
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)

                elif characteristic == 'DENSITY':
                    lw = val - 0.01
                    up = val + 0.01
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)
                elif characteristic == 'ROUGHNESS':
                    lw = val - 2
                    up = val + 2
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)
                elif characteristic == 'LINKED ROUGHNESS':
                    lw = val - 2
                    up = val + 2
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)

            elif toggle == '3':
                compoundBounds = False
                compound_mode = True
                removed_char = param.pop()
                if removed_char == 'THICKNESS':
                    removed_char = 'Thickness'
                elif removed_char == 'DENSITY':
                    removed_char = 'Density'
                elif removed_char == 'ROUGHNESS':
                    removed_char = 'Roughness'
                elif removed_char == 'LINKED ROUGHNESS':
                    removed_char = 'Linked Roughness'

                structDict[param[0]]['compound'].append(removed_char)
                for ele in list(structDict[param[0]]['element'].keys()):
                    structDict[param[0]]['element'][ele].append(removed_char)
            elif toggle == '4':
                cont = False

            if toggle == '1' or toggle == '2':

                print('Select an option: ')
                print('\t 1: Select new parameter for same layer in compound mode')
                print('\t 2: Select new parameter for same layer in element mode')
                print('\t 3: Select a new layer')
                print('\t 4: Select a new parameter type')
                print('\t 5: Return')
                print('\t 6: Exit/Finish')
                toggle = input('\n ->')
                print()
                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4' and toggle != '5' and toggle != '6':
                    toggle = input('Select one of the provided options: ')
                print()
                if toggle == '1':
                    param1 = param.copy()
                    parameters.append(param1)
                    compoundBounds = False
                    compound_mode = True
                    param.pop()
                elif toggle == '2':
                    param1 = param.copy()
                    parameters.append(param1)
                    compoundBounds = False
                    element_mode = True
                    param.pop()
                    param.pop()
                    param.append('ELEMENT')
                elif toggle == '3':
                    param1 = param.copy()
                    parameters.append(param1)
                    param = list()
                    compoundBounds = False
                    layerSelect = True

                elif toggle == '4':
                    parameters.append(param.copy)
                    param = []
                    compoundBounds = False
                    paramType = True
                elif toggle == '5':
                    upperbound.pop()
                    lowerbound.pop()
                elif toggle == '6':
                    parameters.append(param)
                    cont = False
                    isFinished = True

        elif elementBounds:  # compound bounds ---------------------------------------------------------
            print('ELEMENT BOUND SELECTION \n')
            print('Select an option: ')
            print('\t 1: Select parameter boundaries')
            print('\t 2: Use default boundaries')
            print('\t 3: Return')
            print('\t 4: Exit')
            toggle = input('\n -> ')
            print()
            while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4':
                toggle = input('Select one of the provided options:')

            print()
            if toggle == '1':
                bd = input('Enter the parameter optimization boundary as a tuple (lower, upper): ')
                bd = ast.literal_eval(bd)
                bd = (float(bd[0]), float(bd[1]))
                while bd[0]>bd[1] and bd[0] < 0 and type(bd) != tuple:
                    bd = input('Enter the parameter optimization boundary as a tuple (lower, upper) in ascending order: ')
                    bd = ast.literal_eval(bd)
                    bd = (float(bd[0]), float(bd[1]))

                lowerbound.append(bd[0])
                upperbound.append(bd[1])

            elif toggle == '2':
                characteristic = param[-1]
                if characteristic == 'THICKNESS':
                    lw = val - 5
                    up = val + 5
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)

                elif characteristic == 'DENSITY':
                    lw = val - 0.01
                    up = val + 0.01
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)
                elif characteristic == 'ROUGHNESS':
                    lw = val - 2
                    up = val + 2
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)
                elif characteristic == 'LINKED ROUGHNESS':
                    lw = val - 2
                    up = val + 2
                    if lw < 0:
                        lw = 0

                    lowerbound.append(lw)
                    upperbound.append(up)

            elif toggle == '3':
                elementBounds = False
                element_mode = True
                removed_char = param.pop()
                if removed_char == 'THICKNESS':
                    removed_char = 'Thickness'
                elif removed_char == 'DENSITY':
                    removed_char = 'Density'
                elif removed_char == 'ROUGHNESS':
                    removed_char = 'Roughness'
                elif removed_char == 'LINKED ROUGHNESS':
                    removed_char = 'Linked Roughness'

                structDict[param[0]]['compound'].append(removed_char)
                for ele in list(structDict[param[0]]['element'].keys()):
                    structDict[param[0]]['element'][ele].append(removed_char)
            elif toggle == '4':
                cont = False

            if toggle == '1' or toggle == '2':

                print('Select an option: ')
                print('\t 1: Select new parameter for same layer in compound mode')
                print('\t 2: Select new parameter for same layer in element mode')
                print('\t 3: Select a new layer')
                print('\t 4: Select a new parameter type')
                print('\t 5: Return')
                print('\t 6: Exit/Finish')
                toggle = input('\n ->')
                print()
                while toggle != '1' and toggle != '2' and toggle != '3' and toggle != '4' and toggle != '5' and toggle != '6':
                    toggle = input('Select one of the provided options: ')
                print()
                if toggle == '1':
                    parameters.append(param.copy())
                    elementBounds = False
                    compound_mode = True
                    param.pop()  # removes property
                    param.pop()  # removes element
                    param.pop()  # removes mode
                    param.append('COMPOUND')  # replaces element mode with compound mode
                elif toggle == '2':
                    parameters.append(param.copy())
                    elementBounds = False
                    element_mode = True
                    param.pop()  # removes property
                    param.pop()  # removes element

                elif toggle == '3':
                    parameters.append(param.copy())
                    param = list()
                    elementBounds = False
                    layerSelect = True
                elif toggle == '4':
                    parameters.append(param.copy())
                    param = list()
                    elementBounds = False
                    paramType = True

                elif toggle == '5':
                    upperbound.pop()
                    lowerbound.pop()
                elif toggle == '6':
                    parameters.append(param.copy())
                    cont = False
                    isFinished = True

        elif constSelect:
            print('THICKNESS CONSTRAINTS \n')
            print('\t 1: Compound Mode')
            print('\t 2: Element Mode')
            print('\t 3: Return')
            print('\t 4: Exit')
            toggle = input('\n -> ')
            print()
            while toggle not in ['1','2','3','4']:
                toggle = input('Select one of the provided options: ')
            print()
            if toggle == '1':
                constSelect = False
                constCompound = True
            elif toggle == '2':
                constSelect = False
                constElement = True
            elif toggle == '3':
                pass
            elif toggle == '4':
                cont = False
        elif constCompound:
            myList = list()
            for key in list(constDict.keys()):
                if len(constDict[key]) == 3:
                    myList.append(key)

            print('COMPOUND THICKNESS CONSTRAINT \n')
            print('The layers not already selected are -> ' + str(myList))
            toggle = input('Select the range of layers you want to keep a constant thickness: ')
            toggle = toggle.split()
            good = False

            while len(toggle) != 2 and good:
                if len(toggle) != 2:
                    toggle = input('Enter upper and lower range separated by a space: ')
                    toggle = [int(toggle[0]), int(toggle[1])]
                else:
                    toggle = [int(toggle[0]), int(toggle[1])]
                    if range(toggle[0],toggle[1]+1) not in myList:
                        toggle = input('Select a range that exists in -> ' + str(myList))


            param.append('CONSTRAINT')
            param.append('COMPOUND')

            param.append(toggle)

            constCompound = False
            constThick = True

        elif constElement:
            print('ELEMENT THICKNESS CONSTRAINT \n')
            print('Select which element you want to add a constraint: ')
            idx = 1
            temp = dict()
            contEleConst = False
            selected_ele_const = ''
            for ele in sample.myelements:
                printList = []
                for key in list(constDict.keys()):
                    if ele in constDict[key]:
                        printList.append(key)
                if len(printList) != 0:
                    print('\t ' + str(idx) + ': ' + str(ele))
                    temp[str(idx)] = ele
                    idx = idx + 1
            print('\t ' + str(idx) + ': Return')
            temp[str(idx)] = 'Return'
            idx = idx + 1
            print('\t ' + str(idx) + ': Exit')
            temp[str(idx)] = 'Exit'

            toggle = input('\n -> ')
            print()
            while toggle not in list(temp.keys()):
                toggle = input('Select one of the provided options: ')
            print()
            if temp[toggle] == 'Return':
                constElement = False
                constSelect = True
            elif temp[toggle] == 'Exit':
                cont = False
            else:
                contEleConst = True
                selected_ele_const = temp[toggle]

            if contEleConst:
                myList = []
                for key in list(constDict.keys()):
                    if selected_ele_const in constDict[key]:
                        myList.append(key)

                print('The layers not already selected are -> ' + str(myList))
                toggle = input('Select the range of layers you want to keep a constant thickness: ')
                toggle = toggle.split()
                good = False
                while len(toggle) != 2 and good:
                    if len(toggle) != 2:
                        toggle = input('Enter upper and lower range separated by a space: ')
                        toggle = [int(toggle[0]), int(toggle[1])]
                    else:
                        toggle = [int(toggle[0]), int(toggle[1])]
                        if range(toggle[0], toggle[1] + 1) not in myList:
                            toggle = input('Select a range that exists in -> ' + str(myList))

                param.append('CONSTRAINT')
                param.append('ELEMENT')
                param.append(selected_ele_const)

                param.append(toggle)



                constElement = False
                constThick = True

        elif constThick:
            print('SELECT CONSTRAINT THICKNESS \n')

            toggle = input('Select the thickness constraint in Angstrom for layers ' + str(param[-1][0])+' to ' + str(param[-1][1]) + ': ')
            toggle = float(toggle)
            while type(toggle) != float:
                toggle = input('Input must be a number: ')
            param.append(toggle)

            print('\n FINISH CONSTRAINT')
            print()
            print('\t 1: Select another constraint')
            print('\t 2: Select another parameter type')
            print('\t 3: Return')
            print('\t 4: Exit/Finish')

            toggle = input('\n -> ')
            while toggle not in ['1','2','3','4']:
                toggle = input('Select one of the provided options: ')
            print()
            if toggle == '1' or toggle == '2':
                if param[1] == 'COMPOUND':
                    layers = param[2]
                    for i in range(int(layers[0]), int(layers[1])+1):
                        del constDict[i]
                elif param[1] == 'ELEMENT':
                    ele = param[2]
                    layers = param[3]
                    for i in range(int(layers[0]), int(layers[1]) + 1):
                        constDict[i].remove(ele)
                        if len(constDict[i]) == 0:
                            del constDict
                constraints.append(param.copy())

                constThick = False
                if toggle == '1':
                    constSelect = True
                else:
                    paramType = True

            elif toggle == '3':
                constThick = False
                if param[1] == 'COMPOUND':
                    constCompound = True
                elif param[1] == 'ELEMENT':
                    constElement = True
            elif toggle == '4':
                isFinished = True
                cont = False
            param = []

    if not isFinished:
        sys.exit()

    bounds = list(zip(lowerbound, upperbound))
    queue.put([parameters, constraints, bounds])
    return


def getGlobOptParams(fname):
    # Load in the sample information
    sample = ReadSampleHDF5(fname)

    # load in the data and simulation data
    data, data_dict, sim_dict = ReadDataHDF5(fname)  # file must remain open as we process the dataset

    # Running the two tasks at once
    f1 = functools.partial(createTable, data)
    queue = mp.Queue()

    p1 = mp.Process(target=f1)
    p1.start()

    p2 = mp.Process(target=getScans(data, data_dict, sim_dict, queue))
    p2.start()

    p2.join()
    p1.terminate()

    step1 = queue.get()
    scans = step1[0]
    scans = [int(scan) for scan in scans]

    f_mid = functools.partial(saveScans, data, data_dict, sim_dict, scans, sample)
    p_mid = mp.Process(target=f_mid)
    p_mid.start()
    p_mid.join()

    scanBounds = step1[1]
    print(scanBounds)
    queue1 = mp.Queue()
    # show the current scans
    # scans = [1,2,3,4,5]

    f3 = functools.partial(plotScansWidget, sample)
    p3 = mp.Process(target=f3)
    p3.start()


    p5 = mp.Process(target=parameterSelection(sample, queue1))
    p5.start()

    p5.join()
    p3.terminate()

    val = queue1.get()
    parameters = val[0]
    constraints = val[1]
    bounds = val[2]


    x, fun = differential_evolution(fname,scans, parameters, bounds,scanBounds,tolerance=0.00001)

    return x, fun, parameters, scans

def createBoundsDatatype(fname, scans, sBounds, sWeights=None):

        scanBounds = dict()

        info, data_dict, sim_dict = ReadDataHDF5(fname)  # need for data info

        # make sure the number of bounds, number of weights, and scans are all the same
        ns = len(scans)
        nB = len(sBounds)
        if sWeights != None:
            nW = len(sWeights)
            if ns != nB or ns != nW:
                raise SyntaxError('Make sure that the number of bounds, scans, and weights all have the same length.')
        else:
            if ns != nB:
                raise SyntaxError('Make sure that the number of bounds and scans all have the same length.')

        for s in range(ns):
            scan = scans[s]
            scanType = info[scan-1][1]  # retrieve the scan type this
            bound = sBounds[s]  # retrieve the scan's proper bounds
            nb = len(bound)
            if sWeights != None:
                weight = sWeights[s]  # retrieve the scan's proper weights
                nw = len(weight)

                if nb != nw:
                    raise SyntaxError('Make sure every scan has the same number of bounds and weights.')

            else:
                weight = [1 for i in range(nb)]



            # check to make sure that the bounds are in the proper range
            for b in bound:
                up = b[0]
                low = b[1]
                if scanType == 'Reflectivity':

                    if up < 0 or up > 1:
                        raise ValueError('Scan ' + str(scan) + ': Upper bound momentum transfer bounds should be found between 0 and 1')
                    if low < 0 or low> 1:
                        raise ValueError('Scan ' + str(scan) + ': Lower bound momentum transfer bounds should be found between 0 and 1')

                elif scanType == 'Energy':
                    if up < 1 or low < 1:
                        raise ValueError('Bounds for energy scan appear to be set for a reflectivity scan ')

            scanBounds[scan] = (bound, weight)

        return scanBounds

def saveComparisonScanPlots(fname, x, parameters, scans):

    dir = 'comparisonPlots'

    sample = ReadSampleHDF5(fname)  # get the previous sample version
    info, data_dict, sim_dict = ReadDataHDF5(fname)  # get the sample data and simulation data
    newSample = changeSampleParams(x, parameters, sample)  # change to globally optimized parameters

    scans = [s-1 for s in scans]
    info = info[scans]  # retrieve only the needed scans

    figNum = 1
    for idx in range(len(info)):
        scanNumber = info[idx][0]
        scanType = info[idx][1]
        scanName = info[idx][2]

        pol = data_dict[scanName]['Polarization']
        dat = data_dict[scanName]['Data']
        sim = sim_dict[scanName]['Data']
        if scanType == 'Reflectivity':

            E = data_dict[scanName]['Energy']

            qz = dat[0]
            R = dat[2]
            Rsim = sim[2]

            qz, Rnew = newSample.reflectivity(E,qz)
            Rnew = Rnew[pol]

            plt.figure(figNum)
            plt.suptitle('Reflectivity Scan: ' + scanName)
            plt.plot(qz, R)
            plt.plot(qz, Rsim)
            plt.plot(qz, Rnew)
            plt.legend(['Data', 'Current Model', 'Optimized Model'])
            plt.xlabel('Momentum Transfer, qz (A^{-1})')
            if pol == 'S' or pol == 'P' or pol=='LC' or pol=='RC':
                plt.yscale('log')
                plt.ylabel('log(R)')
            else:
                plt.ylabel('A')

        elif scanType == 'Energy':
            Theta = data_dict[scanName]['Angle']

            E = dat[3]
            R = dat[2]
            Rsim = sim[2]

            E, Rnew = newSample.energy_scan(Theta, E)
            Rnew = Rnew[pol]

            plt.figure(figNum)
            plt.suptitle('Energy Scan: ' + scanName)
            plt.plot(E, R)
            plt.plot(E, Rsim)
            plt.plot(E, Rnew)
            plt.legend(['Data', 'Current Model', 'Optimized Model'])
            plt.xlabel('Energy, E (eV)')
            if pol == 'S' or pol == 'P' or pol=='LC' or pol=='RC':
                plt.yscale('log')
                plt.ylabel('log(R)')
            else:
                plt.ylabel('A')
        figNum = figNum + 1
        plt.savefig(dir + '/' + scanName + '.png')

    return

def comparisonScanPlots():
    dir = 'comparisonPlots'
    root = Tk()
    root.geometry('900x900')
    root.title('Show Comparison Scans')
    tabControl = ttk.Notebook(root)

    idx = 0
    im = list()
    for filename in os.listdir(dir):
        frame = ttk.Frame(tabControl)
        imageFile = dir + '/' + filename
        im.append(ImageTk.PhotoImage(Image.open(imageFile)))
        label = Label(frame, imag=im[idx])
        label.pack()
        frame.pack()

        tabControl.add(frame, text=filename)

        idx = idx + 1

    tabControl.pack()
    root.mainloop()

if __name__ == "__main__":
    sample = slab(8)

    sample.addlayer(0, 'SrTiO3', 50, density=[0.027904, 0.027904, 0.083712], roughness=[7.58207, False, 5.77093])
    sample.addlayer(1, 'SrTiO3', 6, density=[0, 0.027904, 0], roughness=[7.58207, 4.03102, 5.77093])

    sample.addlayer(2, 'LaMnO3', 4, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2],
                    linked_roughness=[False, 0.5, False])
    sample.polymorphous(2, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(2, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(3, 'LaMnO3', 17.8, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(3, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(3, ['Mn2+', 'Mn3+'], [0.025, 0], ['Co', 'Ni'])

    sample.addlayer(4, 'LaMnO3', 9, density=[0.021798, 0.0209, 0.084], roughness=[3.77764, 2, 2])
    sample.polymorphous(4, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(4, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(5, 'LaMnO3', 2.5, density=[0.025, 0.024, 0.05], roughness=[0.25, 0.25, 2])
    sample.polymorphous(5, 'Mn', ['Mn2+', 'Mn3+'], [1, 0], sf=['Mn', 'Fe'])
    sample.magnetization(5, ['Mn2+', 'Mn3+'], [0.016, 0], ['Co', 'Ni'])

    sample.addlayer(6, 'LaMnO3', 4.5, density=[0.025, 0.042, 0.04], roughness=[0.25, 0.25, 2])
    sample.polymorphous(6, 'Mn', ['Mn2+', 'Mn3+'], [0.4, 0.6], sf=['Mn', 'Fe'])
    sample.magnetization(6, ['Mn2+', 'Mn3+'], [0.0053, 0], ['Co', 'Ni'])

    sample.addlayer(7, 'CCO', 11.1, density=[0.05, 0.05, 0.01], roughness=2, linked_roughness=[3, 1.5, False])

    fname = 'Pim10uc.h5'

    #sample.plot_density_profile(1)
    #plt.show()
    #WriteSampleHDF5(fname, sample)
    #print(ReadDataHDF5(fname))
    scans = [1,2,3,4]
    data, data_dict, sim_dict = ReadDataHDF5(fname)
    #sample = ReadSampleHDF5(fname)
    #saveScans(data, data_dict, sim_dict, scans, sample)
    #plotScansWidget(data, data_dict, sim_dict, scans, sample)
    x, fun, parameters, scans = getGlobOptParams(fname)
    f6 = functools.partial(saveComparisonScanPlots, fname, x, parameters, scans)
    p6 = mp.Process(target=f6)
    p6.start()
    p6.join()

    p7 = mp.Process(target=comparisonScanPlots)
    p7.start()
    p7.join()

    #parameterSelection(sample, queue)
    #getScans(data, data_dict, sim_dict, queue)
    #results = queue.get()
    #scans = results[0]
    #scanBounds = results[1]

    # now we set the parameters we want to evaluate
    """
    parameters = [[1, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  [2, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  [3, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'],
                  [4, 'STRUCTURAL', 'COMPOUND', 'THICKNESS']]


    lw = [3.5,3.5,17.3,8.5]
    up = [6.5,6.5,19.8,11.5]
    bounds = list(zip(lw, up))
    scans = [1,2,3,4,5,6]

    # determines the bounds for the scans
    sBounds = [[(0.1,0.8)],
               [(0.1,0.3),(0.3,0.5),(0.6,0.8)],
               [(0.1,0.6),(0.7,0.8)],
               [(0.1,0.5)],
               [(0.2,0.6),(0.6,0.8)],
               [(0.1,0.8)]]

    # Determines the weights you want to use for each bound
    sWeights = [[1],
                [1,0.2,0.5],
                [1,0.1],
                [0.5],
                [1,0.8],
                [0.7]]

    scanBounds = createBoundsDatatype(fname, scans, sBounds, sWeights=sWeights)
    start = time.time()
    x, fun = differential_evolution(fname, scans, parameters, bounds, scanBounds, mIter=10, display=True, tolerance=1e-6)
    end = time.time()
    print(end-start)
    
    comparisonScanPlots(fname, x, parameters, scans)
    """
