
import time as time
import numpy as np
import material_structure as ms
import material_model as mm
import scipy.optimize as optimize
import matplotlib.pyplot as plt
import data_structure as ds
import matplotlib.ticker as ticker

def cost_function(x, *args):
    t = args[0]
    y = args[1]
    mode = args[2]

    yp = np.sin(t*x[1])*x[0]
    f = 0
    if mode == 'L1':
        f = sum(np.abs(y-yp))
    elif mode == 'L2':
        f = sum((y-yp)**2)
    elif mode == 'chi':
        f = sum((y-yp)**2/np.absolute(yp))
    elif mode == 'arctan':
        f = np.arctan((y-yp)**2)

    return f


def total_variation(R, Rsim):
    """
    Purpose: Calculate the difference in the total variation between R and Rsim
    :param R: Smoothed reflectivity data
    :param Rsim: Reflectivity simulation
    :return: Difference in total variation
    """
    if len(R) == 0:
        totVar = 0
    else:
        totVar = sum([abs(R[idx+1]-R[idx]) for idx in range(len(R)-1)])/len(R)  # total variation in fitted data

    if len(Rsim) == 0:
        totVarSim = 0
    else:
        totVarSim = sum([abs(Rsim[idx+1]-Rsim[idx]) for idx in range(len(Rsim)-1)])/len(Rsim)  # total variation in simulation

    variation = abs(totVar-totVarSim)
    return variation

def plotting_scans(data_dict, sim_dict, scans, axes,sub,offset=0.0,step=2.0, type='R', xmin=0.0, xmax=100000.0, reverse=False, size=0, top=1.0):

    if reverse:
        scans = scans[::-1]

    my_idx = 0
    if sub[0] != 0 and sub[1] == 0 and size == 0:
        my_idx = sub[0]
    if sub[1] != 0 and sub[0] == 0 and size == 0:
        my_idx = sub[1]
    if type == 'R':

        for i,scan in enumerate(scans):
            qz = data_dict[scan]['Data'][0]
            R = data_dict[scan]['Data'][2]
            Rsim = sim_dict[scan]['Data'][2]
            idx = [i for i in range(len(qz)) if qz[i]>xmin and qz[i] < xmax]

            qz = qz[idx]
            R = R[idx]
            Rsim = Rsim[idx]

            if size != 0:
                axes[sub[0],sub[1]].plot(qz, np.log10(R)+offset+i*step, 'b')
                axes[sub[0],sub[1]].plot(qz, np.log10(Rsim)+offset+i*step, 'r')
            else:
                axes[my_idx].plot(qz, np.log10(R) + offset + i * step, 'b')
                axes[my_idx].plot(qz, np.log10(Rsim) + offset + i * step, 'r')

        if size != 0:
            axes[sub[0],sub[1]].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[sub[0], sub[1]].tick_params(axis='y', labelleft=False)
            axes[sub[0],sub[1]].legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'], frameon=False)
            axes[sub[0], sub[1]].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0], sub[1]].yaxis.set_minor_locator(ticker.AutoMinorLocator())
        else:
            axes[my_idx].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[my_idx].tick_params(axis='y', labelleft=False)
            axes[my_idx].legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'], frameon=False)
            axes[my_idx].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].yaxis.set_minor_locator(ticker.AutoMinorLocator())



    elif type == 'RA':
        for i, scan in enumerate(scans):
            qz = data_dict[scan]['Data'][0]
            R = data_dict[scan]['Data'][2]
            Rsim = sim_dict[scan]['Data'][2]
            idx = [i for i in range(len(qz)) if qz[i] > xmin and qz[i] < xmax]

            qz = qz[idx]
            R = R[idx]/np.linalg.norm(Rsim)
            Rsim = Rsim[idx]/np.linalg.norm(Rsim)

            if size != 0:
                axes[sub[0],sub[1]].plot(qz, R + step * i + offset, 'b')
                axes[sub[0],sub[1]].plot(qz, Rsim + step * i + offset, 'r')
            else:
                axes[my_idx].plot(qz, R + step * i + offset, 'b')
                axes[my_idx].plot(qz, Rsim + step * i + offset, 'r')

        ymax = max(sim_dict[scan]['Data'][2][idx]) + step * i + offset

        if size != 0:
            axes[sub[0],sub[1]].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[sub[0],sub[1]].tick_params(axis='y', labelleft=False)
            axes[sub[0],sub[1]].legend([r'Exp', r'Calc'], loc='upper right', frameon=False)
            axes[sub[0], sub[1]].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0], sub[1]].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0], sub[1]].set_ylim([axes[sub[0], sub[1]].get_ylim()[0], 2])
        else:
            axes[my_idx].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[my_idx].tick_params(axis='y', labelleft=False)
            axes[my_idx].legend([r'Exp', r'Calc'], loc='upper right', frameon=False)
            axes[my_idx].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].set_ylim([axes[my_idx].get_ylim()[0], 1.9])
            #axes[my_idx].set_ylim([axes[my_idx].get_ylim()[0], 0.8])


    elif type == 'E':
        xlim = 0
        ymin = 0
        for i, scan in enumerate(scans):

            E = data_dict[scan]['Data'][3]
            R = data_dict[scan]['Data'][2]
            Rsim = sim_dict[scan]['Data'][2]
            idx = [i for i in range(len(E)) if E[i] > xmin and E[i] < xmax]
            if i==0:
                ymin = min(sim_dict[scan]['Data'][2][idx])

            max1 = max(data_dict[scan]['Data'][2][idx])
            max2 = max(sim_dict[scan]['Data'][2][idx])
            my_max = max(max1, max2)

            E = E[idx]
            R = R[idx]/my_max
            Rsim = Rsim[idx]/my_max

            if size != 0:
                axes[sub[0],sub[1]].plot(E, R + step * i + offset, 'b')
                axes[sub[0],sub[1]].plot(E, Rsim + step * i + offset, 'r')
            else:
                axes[my_idx].plot(E, R + step * i + offset, 'b')
                axes[my_idx].plot(E, Rsim + step * i + offset, 'r')


        ymax = max(sim_dict[scan]['Data'][2][idx])+step*i + offset

        if size != 0:
            axes[sub[0],sub[1]].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[sub[0],sub[1]].tick_params(axis='y', labelleft=False)
            axes[sub[0],sub[1]].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0], sub[1]].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0],sub[1]].legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'], frameon=False)
            axes[sub[0],sub[1]].set_ylim([axes[sub[0],sub[1]].get_ylim()[0], ymax + step * top])
        else:
            axes[my_idx].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[my_idx].tick_params(axis='y', labelleft=False)
            axes[my_idx].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'], frameon=False)
            axes[my_idx].set_ylim([axes[my_idx].get_ylim()[0], ymax + step*top])



    elif type == 'EA':
        ymin = 0
        for i, scan in enumerate(scans):


            E = data_dict[scan]['Data'][3]
            R = data_dict[scan]['Data'][2]
            Rsim = sim_dict[scan]['Data'][2]


            idx = [i for i in range(len(E)) if E[i] > xmin and E[i] < xmax]

            if i==0:
                ymin = min(sim_dict[scan]['Data'][2][idx])

            min1 = abs(min(data_dict[scan]['Data'][2][idx]))
            min2 = abs(min(sim_dict[scan]['Data'][2][idx]))
            max1 = max(data_dict[scan]['Data'][2][idx])
            max2 = max(sim_dict[scan]['Data'][2][idx])

            my_max = max(min1, min2, max1, max2)

            E = E[idx]
            R = np.array(R[idx])/my_max
            Rsim = np.array(Rsim[idx])/my_max

            if size!=0:
                axes[sub[0],sub[1]].plot(E, R + step * i + offset, 'b')
                axes[sub[0],sub[1]].plot(E, Rsim + step * i + offset, 'r')
            else:
                axes[my_idx].plot(E, R + step * i + offset, 'b')
                axes[my_idx].plot(E, Rsim + step * i + offset, 'r')


        ymax = max(sim_dict[scan]['Data'][2][idx]) + step*i + offset

        if size != 0:
            axes[sub[0],sub[1]].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[sub[0],sub[1]].tick_params(axis='y', labelleft=False)
            axes[sub[0],sub[1]].legend([r'Exp', r'Calc'], frameon=False)
            axes[sub[0], sub[1]].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0], sub[1]].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[sub[0],sub[1]].set_ylim([axes[sub[0],sub[1]].get_ylim()[0], ymax + step*top])
        else:
            axes[my_idx].tick_params(axis='both', which='both', direction='in', top=True, right=True)
            axes[my_idx].tick_params(axis='y', labelleft=False)
            axes[my_idx].legend([r'Exp', r'Calc'], frameon=False)
            axes[my_idx].xaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].yaxis.set_minor_locator(ticker.AutoMinorLocator())
            axes[my_idx].set_ylim([axes[my_idx].get_ylim()[0], ymax + step*top])


def return_name_from_scanNum(data, scanNum):

    my_scans = []
    data = data[:,2]
    for d in data:
        check = d.split('_')[0]
        if check in scanNum:
            my_scans.append(d)
    return my_scans

def bar_graph_slice(start, slice_width):
    # properly create the half thickness array
    my_thickness = []
    for i, slice in enumerate(slice_width):
        if i == 0:
            my_thickness.append(start + slice / 2)
        else:
            last = my_thickness[i - 1]
            my_thickness.append(last + slice_width[i - 1] / 2 + slice / 2)

    return my_thickness

def func(z, a, b,sigma,zi):
    import scipy.special as special
    return a * special.erf((z-zi)*np.sqrt(2)/sigma**2)+b

def createError(thickness, value, bounds):
    my_errors = np.ones(len(thickness))*0.000000000000000000000001
    for v,bound in enumerate(bounds):
        low = bound[0]
        up = bound[1]
        idx = [i for i in range(len(thickness)) if thickness[i] >= low and thickness[i]<up]

        my_errors[idx] = value[v]
    return my_errors


if __name__ == "__main__":
    # fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"
    #fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_complete.h5"

    #struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')

    #data, data_dict, sim_dict = ds.LoadDataHDF5(fname)
    #sample = ds.ReadSampleHDF5(fname)
    #sample.energy_shift()
    # Unit Cell Model - 13uc 150K
    fname1 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v2.h5"

    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname2)
    # nR_num = ['64', '57', '69', '71', '73']
    nR_num = ['64', '57', '71']
    R_num = ['79', '60', '62', '75']
    # R_num = ['77', '79', '60', '62', '75']
    # Ti_num = ['68', '70', '72', '74']
    Ti_num = ['68', '72', '74']
    #Mn_num = ['46', '53', '56', '61', '63']
    Mn_num = [ '53', '56', '61', '63']
    # La_num = ['78', '80']
    La_num = ['78']
    R_asymm_num = ['47-48', '51-52', '54-55']
    E_asymm_num = ['49-50', '58-59', '65-66']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant2 = return_name_from_scanNum(data2, nR_num)
    resonant2 = return_name_from_scanNum(data2, R_num)
    Ti_resonance2 = return_name_from_scanNum(data2, Ti_num)
    Mn_resonance2 = return_name_from_scanNum(data2, Mn_num)
    La_resonance2 = return_name_from_scanNum(data2, La_num)
    Mn_R_asymm2 = return_name_from_scanNum(data2, R_asymm_num)
    Mn_E_asymm2 = return_name_from_scanNum(data2, E_asymm_num)

    fig, axes = plt.subplots(1, 2)
    plotting_scans(data_dict1, sim_dict1, Mn_resonance2, axes, (0, 0), offset=0, step=1, xmin=635, xmax=660,
                   type='E', reverse=True, size=0, top=1.25)
    plotting_scans(data_dict2, sim_dict2, Mn_resonance2, axes, (0, 1), offset=0, step=1, xmin=635, xmax=660,
                   type='E', reverse=True, size=0, top=1.25)
    #plotting_scans(data_dict2, sim_dict2, Mn_E_asymm2, axes, (0, 1), offset=0, step=3, xmin=635,
    #               type='EA', reverse=True, size=0, top=0.75)

    axes[0].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)
    axes[1].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)

    axes[0].set_xlabel(r'Energy, E (eV)', fontsize=16)
    axes[1].set_xlabel(r'Energy, E (eV)', fontsize=16)
    # shared_x.yaxis.set_label_coords(-0.02, 0.5)
    plt.tight_layout()

    fig, axes = plt.subplots(1, 2)
    plotting_scans(data_dict1, sim_dict1, resonant2, axes, (0, 0), offset=0, step=3, xmin=0.03,
                   type='R', reverse=True, size=0, top=0.5)
    plotting_scans(data_dict2, sim_dict2, resonant2, axes, (0, 1), offset=0, step=3,  xmin=0.03,
                   type='R', reverse=True, size=0, top=0.5)
    # plotting_scans(data_dict2, sim_dict2, Mn_E_asymm2, axes, (0, 1), offset=0, step=3, xmin=635,
    #               type='EA', reverse=True, size=0, top=0.75)

    axes[0].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)
    axes[1].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)

    axes[0].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA}$)', fontsize=16)
    axes[1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA}$)', fontsize=16)
    # shared_x.yaxis.set_label_coords(-0.02, 0.5)
    plt.tight_layout()

    fig, axes = plt.subplots(1, 2)
    plotting_scans(data_dict1, sim_dict1, Mn_R_asymm2, axes, (0, 0), offset=0, step=0.75, xmin=0.03, xmax=0.5,
                   type='RA', reverse=True, size=0, top=1)
    plotting_scans(data_dict2, sim_dict2, Mn_R_asymm2, axes, (0, 1), offset=0, step=0.75, xmin=0.03, xmax=0.5,
                   type='RA', reverse=True, size=0, top=1)
    # plotting_scans(data_dict2, sim_dict2, Mn_E_asymm2, axes, (0, 1), offset=0, step=3, xmin=635,
    #               type='EA', reverse=True, size=0, top=0.75)

    axes[0].set_ylabel(r'Circular Polarized Asymmetry (arb. units)', fontsize=16)
    axes[1].set_ylabel(r'Circular Polarized Asymmetry (arb. units)', fontsize=16)

    axes[0].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA}$)', fontsize=16)
    axes[1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA}$)', fontsize=16)
    # shared_x.yaxis.set_label_coords(-0.02, 0.5)
    plt.tight_layout()
    plt.show()
    """
    fig, axes = plt.subplots(2, 1)
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1.h5"

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_mag = sample.density_profile()

    # thickness_shift = -11.75
    # thickness = thickness + thickness_shift

    zero = np.zeros(len(thickness))

    xlim = [-20, 75]

    slice_width = [3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.8915, 3.878, 3.878, 3.878, 3.878, 3.878,
                   3.878,
                   3.878, 3.878, 3.878, 3.878, 3.878, 3.878, 1.939]
    start = -11.715
    offset = -15.75
    bar_width = np.array(slice_width) - 0.3

    # properly create the half thickness array
    my_thickness = bar_graph_slice(start, slice_width)

    Ti = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Mn2 = [0, 0, 0, 0, 0, 0, 0, 0.06314187353650001, 0.11028773825175128, 0.02275647285959691, 0.3701310880359007,
           0, 0, 0, 0, 0, 0, 0, 0, 0.030718436673865812, 0.19552804020068865, 0.29622505384943193]
    Mn3 = [0, 0, 0, 0, 0, 0, 0, 6.054978809009333e-28, 0.01843395509529002, 0.3974431333219987, 0.5284835745428983,
           0.8096132735726913, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7339447327644624, 0.8044719526546098, 0.703774946150568]
    Mn4 = [0, 0, 0, 0, 0, 0, 0, 1.148911588089163e-07, 0.03322193763257628, 1.8224808898128161e-07, 0.09486627140683755,
           0.19038672642730867, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.23533683056167182, 7.144701336509033e-09, 0]
    mag = [0, 0, 0, 0, 0, 0, 0, 0.0004003354023059758, 0.001816557185740252, 0.00658896207089925, 0.012848764938116344,
           0.016830257935419913, 0.016830257935419913, 0.016830257935419913, 0.016830257935419913, 0.016830257935419913,
           0.016830257935419913, 0.016830257935419913, 0.016830257935419913, 0.01070931743712782, 0.005858933727020148,
           0]

    mag = -np.array(mag) / max(mag)

    Arho = [0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,
            0.028, 0.028,
            0.028, 0.028, 0.028, 0.028, 0.03319719939671565]
    Mn3 = np.array(Mn2) + np.array(Mn3)
    Mn4 = np.array(Mn3) + np.array(Mn4)

    Ti = np.array(Ti) * np.array(Arho) / 0.028
    Mn2 = np.array(Mn2) * np.array(Arho) / 0.028
    Mn3 = np.array(Mn3) * np.array(Arho) / 0.028
    Mn4 = np.array(Mn4) * np.array(Arho) / 0.028

    axes[1].bar(np.array(my_thickness) + offset, Ti, color='yellow', edgecolor='black', width=bar_width, label='Ti')
    axes[1].bar(np.array(my_thickness) + offset, Mn4, color='purple', edgecolor='black', width=bar_width,
                label=r'$\mathrm{Mn^{4+}}$')
    axes[1].bar(np.array(my_thickness) + offset, Mn3, color='grey', edgecolor='black', width=bar_width,
                label=r'$\mathrm{Mn^{3+}}$')
    axes[1].bar(np.array(my_thickness) + offset, Mn2, color='green', edgecolor='black', width=bar_width,
                label=r'$\mathrm{Mn^{2+}}$')

    axes[1].bar(np.array(my_thickness) + offset, mag, color='magenta', edgecolor='black', width=bar_width,
                label='Mag')
    axes[1].axhline(0, color='black', linestyle='--')
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)

    # plt.gca().invert_xaxis()
    axes[1].set_xlim(xlim)
    axes[1].set_ylim(axes[1].get_ylim()[0], 1.95)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].legend(frameon=False)

    # Asites
    Sr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    La = [0, 0, 0, 0, 0, 0, 0.08081520235245354, 0.31540659890843226, 0.6020890684884148, 0.7591457496086605,
          0.6916283063660732, 0.6524714534498322,
          0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.3709629905970221, 0.033554501863506725, 0]

    Brho = [0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,
            0.028,
            0.028, 0.028, 0.028, 0.028, 0.02517789271095681, 0]

    Sr = np.array(Sr) * np.array(Brho) / 0.028
    La = np.array(La) * np.array(Brho) / 0.028
    O = density['O'] / 0.028
    C = density['C2'] / 0.028

    axes[0].bar(np.array(my_thickness) + offset, Sr, color='cyan', edgecolor='black', width=bar_width, label='Sr')
    axes[0].bar(np.array(my_thickness) + offset, La, color='blue', edgecolor='black', width=bar_width, label='La')
    axes[0].plot(thickness + offset + 3.878 + 1.939 * 2, O, 'r', label='O')
    axes[0].plot(thickness + offset + 3.878 + 1.939 * 2, C, color='orange', label='C')
    axes[0].axhline(0, color='black', linestyle='--')
    axes[0].set_xlim(xlim)

    # axes[2,0].gca().set_xticklabels([])
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0].set_ylim(-0.1, axes[0].get_ylim()[1])
    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].legend(frameon=False)
    axes[0].set_xticklabels([])
    axes[1].set_xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=20)
    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Fractional Site Occupation', fontsize=20)

    axes[0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_y_axes().join(axes[0], shared_y)
    axes[0].get_shared_y_axes().join(axes[1], shared_y)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    # axes[0,0].grid(True)
    # plt.gca().invert_xaxis()
    plt.show()
    """
    """

    # Unit Cell Model - 13uc 150K
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1.h5"

    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname)
    # nR_num = ['64', '57', '69', '71', '73']
    nR_num = ['64', '57', '71']
    R_num = ['79', '60', '62', '75']
    # R_num = ['77', '79', '60', '62', '75']
    # Ti_num = ['68', '70', '72', '74']
    Ti_num = ['68', '72', '74']
    Mn_num = ['46', '53', '56', '61', '63']
    # La_num = ['78', '80']
    La_num = ['78']
    R_asymm_num = ['47-48', '51-52', '54-55']
    E_asymm_num = ['49-50', '58-59', '65-66']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant2 = return_name_from_scanNum(data2, nR_num)
    resonant2 = return_name_from_scanNum(data2, R_num)
    Ti_resonance2 = return_name_from_scanNum(data2, Ti_num)
    Mn_resonance2 = return_name_from_scanNum(data2, Mn_num)
    La_resonance2 = return_name_from_scanNum(data2, La_num)
    Mn_R_asymm2 = return_name_from_scanNum(data2, R_asymm_num)
    Mn_E_asymm2 = return_name_from_scanNum(data2, E_asymm_num)

    fig, axes = plt.subplots(1,2)
    plotting_scans(data_dict2, sim_dict2, Mn_resonance2, axes, (0,0), offset=0, step=1.5, xmin=635, xmax=660,
                   type='E', reverse=True, size=0, top=1.5)
    plotting_scans(data_dict2, sim_dict2, Mn_E_asymm2, axes, (0, 1), offset=0, step=3, xmin=635,
                   type='EA', reverse=True, size=0, top=0.75)



    axes[0].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)
    axes[1].set_ylabel(r'Circular Polarized Asymmetry (arb. units)', fontsize=16)

    axes[0].set_xlabel(r'Energy, E (eV)', fontsize=16)
    axes[1].set_xlabel(r'Energy, E (eV)', fontsize=16)
    #shared_x.yaxis.set_label_coords(-0.02, 0.5)
    plt.tight_layout()
    plt.show()
    """
    """
    #   -------------------- Mn2+ removed comparison!!  --------------------------------------- #
    fname1 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1_interface.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1_surface.h5"
    fname3 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1_Mn3.h5"

    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname2)
    data3, data_dict3, sim_dict3 = ds.LoadDataHDF5(fname3)

    my_scans = ['53','61']
    keys = return_name_from_scanNum(data2, my_scans)

    fig, axes = plt.subplots(1, 3)

    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel('Energy, E (eV)', fontsize=16)

    plotting_scans(data_dict1, sim_dict1, keys, axes, (0,0), offset=0, step=1, xmin=632, xmax=662,
                   type='E', reverse=True, size=0, top=1.25)
    plotting_scans(data_dict2, sim_dict2, keys, axes, (0, 1), offset=0, step=1, xmin=632, xmax=662,
                   type='E', reverse=True, size=0, top=1.25)
    plotting_scans(data_dict3, sim_dict3, keys, axes, (0, 2), offset=0, step=1, xmin=632, xmax=662,
                   type='E', reverse=True, size=0, top=1.25)


    axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_x_axes().join(axes[0], shared_x)
    axes[1].get_shared_x_axes().join(axes[1], shared_x)
    axes[2].get_shared_x_axes().join(axes[1], shared_x)

    # axes[0, 0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    # axes[1, 0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=12)
    axes[0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    #shared_x.yaxis.set_label_coords(-0.01, 0.5)

    fname1 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete_interface.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete_surface.h5"
    fname3 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete_Mn3.h5"


    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname2)
    data3, data_dict3, sim_dict3 = ds.LoadDataHDF5(fname3)


    my_scans = ['37', '43']
    keys = return_name_from_scanNum(data2, my_scans)

    fig, axes = plt.subplots(1, 3)

    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel('Energy, E (eV)', fontsize=16)

    plotting_scans(data_dict1, sim_dict1, keys, axes, (0, 0), offset=0, step=1, xmin=632, xmax=662,
                   type='E', reverse=True, size=0, top=1.25)
    plotting_scans(data_dict2, sim_dict2, keys, axes, (0, 1), offset=0, step=1, xmin=632, xmax=662,
                   type='E', reverse=True, size=0, top=1.25)
    plotting_scans(data_dict3, sim_dict3, keys, axes, (0, 2), offset=0, step=1, xmin=632, xmax=662,
                   type='E', reverse=True, size=0, top=1.25)


    axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[2].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_x_axes().join(axes[0], shared_x)
    axes[1].get_shared_x_axes().join(axes[1], shared_x)
    axes[2].get_shared_x_axes().join(axes[2], shared_x)


    # axes[0, 0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    # axes[1, 0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=12)
    axes[0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    # shared_x.yaxis.set_label_coords(-0.01, 0.5)

    fname1 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1_mag.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1.h5"

    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname2)

    my_scans = ['47-48', '51-52', '54-55']
    keys = return_name_from_scanNum(data2, my_scans)

    fig, axes = plt.subplots(1, 2)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel(r'Momentum Transfer,$\mathrm{qz}$ ($\mathrm{\AA}$)', fontsize=16)

    plotting_scans(data_dict1, sim_dict1, keys, axes, (0, 0), offset=0, step=0.75, xmin=0.03, xmax=0.5,
                   type='RA', reverse=True, size=0, top=1)
    plotting_scans(data_dict2, sim_dict2, keys, axes, (0, 1), offset=0, step=0.75, xmin=0.03, xmax=0.5,
                   type='RA', reverse=True, size=0, top=1)

    #axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[0].get_shared_x_axes().join(axes[0], shared_x)
    #axes[1].get_shared_x_axes().join(axes[1], shared_x)

    #axes[0, 0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    axes[0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=16)
    axes[0].set_xlabel('Momentum Transfer,$\mathrm{qz}$ ($\mathrm{\AA}$)', fontsize=16)
    #axes[0].set_ylabel('Circular Polarized Asymmetry (arb. units)', fontsize=16)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    # shared_x.yaxis.set_label_coords(-0.01, 0.5)

    fname1 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete_mag.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"

    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname2)

    my_scans = ['111-112', '117-118']
    keys = return_name_from_scanNum(data2, my_scans)

    fig, axes = plt.subplots(1, 2)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel(r'Momentum Transfer,$\mathrm{qz}$ ($\mathrm{\AA}$)', fontsize=16)

    plotting_scans(data_dict1, sim_dict1, keys, axes, (0, 0), offset=0, step=0.5, xmin=0.03, xmax=0.5,
                   type='RA', reverse=True, size=0, top=1)
    plotting_scans(data_dict2, sim_dict2, keys, axes, (0, 1), offset=0, step=0.5, xmin=0.03, xmax=0.5,
                   type='RA', reverse=True, size=0, top=1)

    #axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[0].get_shared_x_axes().join(axes[0], shared_x)
    #axes[1].get_shared_x_axes().join(axes[1], shared_x)

    # axes[0, 0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    # axes[1, 0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=12)
    axes[0].set_ylabel('Circular Polarized Asymmetry (arb. units)', fontsize=16)
    axes[0].set_xlabel(r'Momentum Transfer,$\mathrm{qz}$ ($\mathrm{\AA}$)', fontsize=16)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    # shared_x.yaxis.set_label_coords(-0.01, 0.5)
    plt.show()
    
    """
    """
    #  ------------------ Remove the magnetic dead layer
    
    fig, axes = plt.subplots(2,1)
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_unitCell_v1.h5"

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_mag = sample.density_profile()

    # thickness_shift = -11.75
    # thickness = thickness + thickness_shift

    zero = np.zeros(len(thickness))

    xlim = [-20, 75]

    slice_width = [3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.8915, 3.878, 3.878, 3.878, 3.878, 3.878, 3.878,
                   3.878, 3.878, 3.878, 3.878, 3.878, 3.878, 1.939]
    start = -11.715
    offset = -15.75
    bar_width = np.array(slice_width) - 0.3

    # properly create the half thickness array
    my_thickness = bar_graph_slice(start, slice_width)

    Ti = [1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    Mn2 = [0,0, 0, 0, 0, 0, 0, 0.06314187353650001, 0.11028773825175128, 0.02275647285959691, 0.3701310880359007,
           0,0,0,0,0,0,0,0,0.030718436673865812,0.19552804020068865,0.29622505384943193]
    Mn3 = [0,0, 0, 0, 0, 0, 0, 1.1489115883284551e-07, 0.05165589272786634, 0.3974433155700877, 0.6233498459497357,
           1, 1,1,1,1,1,1,1,0.9692815633261342,0.8044719597993113,0.703774946150568]
    mag = [0,0, 0, 0, 0, 0, 0, 0.0004003354023059758, 0.001816557185740252, 0.00658896207089925, 0.012848764938116344,
           0.016830257935419913, 0.016830257935419913,0.016830257935419913,0.016830257935419913,0.016830257935419913,
           0.016830257935419913,0.016830257935419913,0.016830257935419913,0.01070931743712782,0.005858933727020148,
           0]

    mag = -np.array(mag) / max(mag)

    Arho = [0.028,0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,0.028,0.028,0.028,0.028,0.028,
            0.028,0.028,0.028,0.028,0.03319719939671565]
    Mn3 = np.array(Mn2) + np.array(Mn3)

    Ti = np.array(Ti) * np.array(Arho) / 0.028
    Mn2 = np.array(Mn2) * np.array(Arho) / 0.028
    Mn3 = np.array(Mn3) * np.array(Arho) / 0.028

    axes[1].bar(np.array(my_thickness) + offset, Ti, color='yellow', edgecolor='black', width=bar_width, label='Ti')
    axes[1].bar(np.array(my_thickness) + offset, Mn3, color='grey', edgecolor='black', width=bar_width,
                   label=r'$\mathrm{Mn^{3.3+}}$')
    axes[1].bar(np.array(my_thickness) + offset, Mn2, color='green', edgecolor='black', width=bar_width,
                   label=r'$\mathrm{Mn^{2+}}$')
    axes[1].bar(np.array(my_thickness) + offset, mag, color='magenta', edgecolor='black', width=bar_width,
                   label='Mag')
    axes[1].axhline(0, color='black', linestyle='--')
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)

    # plt.gca().invert_xaxis()
    axes[1].set_xlim(xlim)
    axes[1].set_ylim(axes[1].get_ylim()[0], 1.95)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].legend(frameon=False)

    # Asites
    Sr = [1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
    La = [0,0, 0, 0, 0, 0, 0.08081520235245354,0.31540659890843226,0.6020890684884148,0.7591457496086605, 0.6916283063660732,0.6524714534498322,
          0.7,0.7,0.7,0.7,0.7,0.7,0.7,0.3709629905970221,0.033554501863506725,0]

    Brho = [0.028,0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,0.028,0.028,0.028,0.028,
            0.028,0.028,0.028,0.028,0.02517789271095681,0]

    Sr = np.array(Sr) * np.array(Brho) / 0.028
    La = np.array(La) * np.array(Brho) / 0.028
    O = density['O'] / 0.028
    C = density['C2'] / 0.028

    axes[0].bar(np.array(my_thickness) + offset, Sr, color='cyan', edgecolor='black', width=bar_width, label='Sr')
    axes[0].bar(np.array(my_thickness) + offset, La, color='blue', edgecolor='black', width=bar_width, label='La')
    axes[0].plot(thickness + offset +3.878 +1.939*2, O, 'r', label='O')
    axes[0].plot(thickness + offset+3.878 +1.939*2, C, color='orange', label='C')
    axes[0].axhline(0, color='black', linestyle='--')
    axes[0].set_xlim(xlim)

    # axes[2,0].gca().set_xticklabels([])
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0].set_ylim(-0.1, axes[0].get_ylim()[1])
    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].legend(frameon=False)
    axes[0].set_xticklabels([])
    axes[1].set_xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=20)
    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Fractional Site Occupation', fontsize=20)

    axes[0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_y_axes().join(axes[0], shared_y)
    axes[0].get_shared_y_axes().join(axes[1], shared_y)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    # axes[0,0].grid(True)
    # plt.gca().invert_xaxis()
    plt.show()
    

    """
    """
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim7uc_unitCell_complete.h5"

    data, data_dict, sim_dict = ds.LoadDataHDF5(fname)

    name1 = "38_455.73_S"
    name2 = "42_459.75_S"

    dat1 = data_dict[name1]['Data']
    dat2 = data_dict[name2]['Data']

    fig, axes = plt.subplots(1,2)
    axes[0].plot(dat1[0], np.log10(dat1[2]))
    axes[1].plot(dat2[0], np.log10(dat2[2]))

    axes[0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)

    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].tick_params(axis='y', labelleft=False)
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].tick_params(axis='y', labelleft=False)
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)

    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA}$)', fontsize=16)

    axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_x_axes().join(axes[0], shared_x)
    axes[1].get_shared_x_axes().join(axes[1], shared_x)
    plt.tight_layout()
    plt.show()
    """
    """
    # ---------------------------------- Demonstrate continuous and uc model ----------------------------------------
    fname1 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_v4.h5"

    sample1 = ds.ReadSampleHDF5(fname1)
    sample2 = ds.ReadSampleHDF5(fname2)

    thickness1, density1, temp = sample1.density_profile()
    thickness2, density2, temp2 = sample2.density_profile()

    Sr1 = density1['Sr']
    La1 = density1['La']
    Sr2 = density2['Sr']
    La2 = density2['La']

    fig, axes = plt.subplots(1,2)

    axes[0].plot(thickness1, Sr1,label='Sr')
    axes[0].plot(thickness1, La1,'--' , label='La')
    axes[1].plot(thickness2, Sr2, label='Sr')
    axes[1].plot(thickness2, La2,'--' , label='La')
    axes[0].set_ylabel(r'Density, $\rho$ ($\mathrm{mol/cm^{3}}$)', fontsize=16)

    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].tick_params(axis='y', labelleft=False)
    axes[0].legend(frameon=False)
    axes[1].legend(frameon=False)

    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=16)

    axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_x_axes().join(axes[0], shared_x)
    axes[1].get_shared_x_axes().join(axes[1], shared_x)

    plt.tight_layout()
    plt.show()
    """

    # ---------------------------------- Saving for Emma --------------------------------------------------------------
    """
    fname1 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_complete.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_complete.h5"
    fname3 = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-300K_complete.h5"

    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname2)
    data3, data_dict3, sim_dict3 = ds.LoadDataHDF5(fname3)

    # EM1 150
    R_num = ['35', '22', '24', '33']
    nR_num = ['26', '19', '31']
    Ti_num = ['29', '30', '32']
    Mn_num = ['8', '15', '18', '23', '25']
    La_num = ['34']
    R_asymm_num = ['9-10', '13-14', '16-17']
    E_asymm_num = ['11-12', '20-21', '27-28']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    resonant1 = return_name_from_scanNum(data1, R_num)
    non_resonant1 = return_name_from_scanNum(data1, nR_num)
    Ti_resonance1 = return_name_from_scanNum(data1, Ti_num)
    Mn_resonance1 = return_name_from_scanNum(data1, Mn_num)
    La_resonance1 = return_name_from_scanNum(data1, La_num)
    Mn_R_asymm1 = return_name_from_scanNum(data1, R_asymm_num)
    Mn_E_asymm1 = return_name_from_scanNum(data1, E_asymm_num)

    # EM2 - 150
    # nR_num = ['64', '57', '69', '71', '73']
    nR_num2 = ['64', '57', '71']
    R_num2 = ['79', '60', '62', '75']
    # R_num = ['77', '79', '60', '62', '75']
    # Ti_num = ['68', '70', '72', '74']
    Ti_num2 = ['68', '72', '74']
    Mn_num2 = ['46', '53', '56', '61', '63']
    # La_num = ['78', '80']
    La_num2 = ['78']
    R_asymm_num2 = ['47-48', '51-52', '54-55']
    E_asymm_num2 = ['49-50', '58-59', '65-66']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant2 = return_name_from_scanNum(data2, nR_num2)
    resonant2 = return_name_from_scanNum(data2, R_num2)
    Ti_resonance2 = return_name_from_scanNum(data2, Ti_num2)
    Mn_resonance2 = return_name_from_scanNum(data2, Mn_num2)
    La_resonance2 = return_name_from_scanNum(data2, La_num2)
    Mn_R_asymm2 = return_name_from_scanNum(data2, R_asymm_num2)
    Mn_E_asymm2 = return_name_from_scanNum(data2, E_asymm_num2)

    # EM2 -300
    # nR_num = ['26', '29', '31', '33', '35']
    nR_num3 = ['26', '29', '33']
    # R_num = ['39', '41', '22', '24', '37']
    R_num3 = ['41', '22', '24', '37']
    # Ti_num = ['30', '32', '34', '36']
    Ti_num3 = ['30', '34', '36']
    Mn_num3 = ['8', '15', '18', '23', '25']
    # La_num = ['40', '42']
    La_num3 = ['40']
    R_asymm_num3 = ['9-10', '13-14', '16-17']
    E_asymm_num3 = ['11-12', '20-21', '27-28']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant3 = return_name_from_scanNum(data3, nR_num3)
    resonant3 = return_name_from_scanNum(data3, R_num3)
    Ti_resonance3 = return_name_from_scanNum(data3, Ti_num3)
    Mn_resonance3 = return_name_from_scanNum(data3, Mn_num3)
    La_resonance3 = return_name_from_scanNum(data3, La_num3)
    Mn_R_asymm3 = return_name_from_scanNum(data3, R_asymm_num3)
    Mn_E_asymm3 = return_name_from_scanNum(data3, E_asymm_num3)

    # find the names of the data used and then save it as a csv file
    import csv
    E1_150 = non_resonant1
    my_data_1 = list()
    for name in E1_150:
        E = data_dict1[name]['Energy']
        #Theta = data_dict1[name]['Angle']
        n1 = ['qz - ' + str(E)]
        n2 = ['Rexp - ' + str(E)]
        n3 = ['Rsim - ' +str(E)]

        test1 = n1 + data_dict1[name]['Data'][0].tolist()
        test2 = n2 + data_dict1[name]['Data'][2].tolist()
        test3 = n3 + sim_dict1[name]['Data'][2].tolist()

        my_data_1.append(test1)
        my_data_1.append(test2)
        my_data_1.append(test3)

    my_data_1 = np.array(my_data_1).transpose()

    with open("LSMO-E1-150-nonresonant.csv", "w", newline='') as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerows(my_data_1)

    E2_150 = non_resonant2
    my_data_2 = list()
    for name in E2_150:
        E = data_dict2[name]['Energy']
        #Theta = data_dict2[name]['Angle']

        n1 = ['qz - ' + str(E)]
        n2 = ['Rexp - ' + str(E)]
        n3 = ['Rsim - ' + str(E)]

        test1 = n1 + data_dict2[name]['Data'][0].tolist()
        test2 = n2 + data_dict2[name]['Data'][2].tolist()
        test3 = n3 + sim_dict2[name]['Data'][2].tolist()

        my_data_2.append(test1)
        my_data_2.append(test2)
        my_data_2.append(test3)

    my_data_2 = np.array(my_data_2).transpose()
    with open("LSMO-E2-150-nonresonant.csv", "w", newline='') as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerows(my_data_2)

    E2_300 = non_resonant3
    my_data_3 = list()
    for name in E2_300:
        E = data_dict3[name]['Energy']
        #Theta = data_dict3[name]['Angle']

        n1 = ['qz - ' + str(E)]
        n2 = ['Rexp - ' + str(E)]
        n3 = ['Rsim - ' + str(E)]

        test1 = n1 + data_dict3[name]['Data'][0].tolist()
        test2 = n2 + data_dict3[name]['Data'][2].tolist()
        test3 = n3 + sim_dict3[name]['Data'][2].tolist()

        my_data_3.append(test1)
        my_data_3.append(test2)
        my_data_3.append(test3)

    my_data_3 = np.array(my_data_3).transpose()

    with open("LSMO-E2-300-nonresonant.csv", "w", newline='') as my_csv:
        csvWriter = csv.writer(my_csv, delimiter=',')
        csvWriter.writerows(my_data_3)
    """
    """

    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_complete.h5"
    sample4 = ds.ReadSampleHDF5(fname)
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim7uc_unitCell_complete.h5"
    sample7 = ds.ReadSampleHDF5(fname)
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"
    sample10 = ds.ReadSampleHDF5(fname)

    Total_Mn2_4 = 0
    vals = [0.00252081, 0.00067035, 0.00916918,0.00252081]
    error = 0
    idx = 0
    for i in range(2,6):
        Total_Mn2_4 = Total_Mn2_4 + sample4.structure[i]['B'].poly_ratio[1]*sample4.structure[i]['B'].thickness*1e-8
        error = error + (sample4.structure[i]['B'].poly_ratio[1]*0.1*0.028*1e-8)**2 + (sample4.structure[i]['B'].thickness*vals[idx]*0.028*1e-8)**2
        idx = idx + 1

    Total_Mn2_4 = Total_Mn2_4*(1e-8)**2*(3.905)**2*6.02214076e23
    error = np.sqrt(error)*(1e-8)**2*(3.905)**2*6.02214076e23
    print('4uc: ' + str(Total_Mn2_4*0.028) + ' +- ' + str(error))

    Total_Mn2_7 = 0
    vals = [0.00614008, 0.00344086, 0.04148664, 0.01746246, 0.02394579]
    error = 0
    idx = 0

    for i in range(2, 7):
        Total_Mn2_7 = Total_Mn2_7 + sample7.structure[i]['B'].poly_ratio[1] * sample7.structure[i]['B'].thickness*1e-8
        error = error + (sample7.structure[i]['B'].poly_ratio[1] * 0.1 * 0.028*1e-8) ** 2 + (
                    sample7.structure[i]['B'].thickness * vals[idx] * 0.028*1e-8) ** 2
        idx = idx + 1

    Total_Mn2_7 = Total_Mn2_7 * (1e-8) ** 2 * (3.905) ** 2 * 6.02214076e23
    error = np.sqrt(error) * (1e-8) ** 2 * (3.905) ** 2 * 6.02214076e23
    print('7uc: ' + str(Total_Mn2_7 * 0.028)+' +- ' + str(error))

    Total_Mn2_10 = 0
    vals = [0.00345445, 0.00229006, 0.01244586, 0.02939831, 0.02136382]
    error = 0
    idx = 0

    for i in range(2, 7):
        Total_Mn2_10 = Total_Mn2_10 + sample10.structure[i]['B'].poly_ratio[1] * sample10.structure[i]['B'].thickness*1e-8
        error = error + (sample7.structure[i]['B'].poly_ratio[1] * 0.1 * 0.028 * 1e-8) ** 2 + (
                sample7.structure[i]['B'].thickness * vals[idx] * 0.028 * 1e-8) ** 2
        idx = idx + 1

    Total_Mn2_10 = Total_Mn2_10 * (1e-8) ** 2 * (3.905) ** 2 * 6.02214076e23
    error = np.sqrt(error) * (1e-8) ** 2 * (3.905) ** 2 * 6.02214076e23
    print('10uc: ' + str(Total_Mn2_10 * 0.028)+' +- ' + str(error))
    """
    """
    # [qz_data,R_data,qz_sim,R_sim, random_points, predictions]

    # ---------------------------------- MAGNETIC ERROR
    
    rho = np.array([0.0001502703621304594, 0.0018479773003216926, 0.00712558371573116, 0.0176497168459156, 0.003261994045717825])/10  # scaled by 10 in GUI
    d = [3.905, 3.94, 3.97, 3.97, 3.97]

    delta_rho = np.array([3.13099625e-05, 1.38530565e-04, 1.91284760e-04, 1.65206812e-04, 7.06123230e-05])/10
    delta_d = [1.69973443, 0.53498416, 1.03517472, 0.4971241, 0.49950181]

    D = 0
    E = 0
    for i in range(len(rho)):
        D = D + rho[i]*d[i]*1e-8
        E = E + (delta_rho[i]*d[i]*1e-8)**2 + (rho[i]*delta_d[i]*1e-8)**2

    E = np.sqrt(E)* (1e-8) ** 2 * (3.97) ** 2 * 6.02214076e23
    D = D* (1e-8) ** 2 * (3.97) ** 2 * 6.02214076e23

    print('Total Magnetism (4uc) = ' + str(D) + ' +- ' + str(E))

    rho = np.array([0.0003066750126945097, 0.0002068731239481137, 0.016932058948664633, 0.02655782307499808, 0.0314329491242263, 0.0314329491242263, 0.0314329491242263, 0.0314329491242263,
                    0.00010936761100236586]) / 10  # scaled by 10 in GUI
    d = [3.905, 3.94, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97]

    #delta_rho = np.array([3.38396732e-05, 3.75170965e-05, 5.94321782e-04, 9.17391043e-04,
    #                      1.42100443e-03, 1.70273803e-03, 2.55402673e-03, 3.23024690e-03,
    #                      4.57742480e-04]) / 10
    delta_rho = np.array([3.38396732e-05, 3.75170965e-05, 5.94321782e-04, 9.17391043e-04,
                          1.42100443e-03, 1.70273803e-03, 2.55402673e-03, 3.23024690e-03,
                          4.57742480e-04]) / 10
    delta_d = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1]

    D1 = 0
    E1 = 0
    for i in range(len(rho)):
        D1 = D1 + rho[i] * d[i]*1e-8
        E1 = E1 + (delta_rho[i] * d[i]*1e-8) ** 2 + (rho[i] * delta_d[i]*1e-8) ** 2

    E1 = np.sqrt(E1) * (1e-8) ** 2 * (3.97) ** 2 * 6.02214076e23
    D1 = D1 * (1e-8) ** 2 * (3.97) ** 2 * 6.02214076e23

    print('Total Magnetism (7uc) = ' + str(D1) + ' +- ' + str(E1))

    rho = np.array([0.00019341730101137163, 0.000186586051446197, 0.014676919509906398, 0.027492584861078195,
                    0.031889765641913245, 0.031889765641913245, 0.031889765641913245, 0.031889765641913245,
                    0.031889765641913245, 0.031889765641913245, 0.012130795599578477, 0.009068991034701678]) / 10  # scaled by 10 in GUI
    d = [3.905, 3.94, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97]

    #delta_rho = np.array([1.49852433e-05, 2.14513851e-05, 1.43998028e-04, 1.15682387e-04,
    #                      1.30548200e-04, 1.30548200e-04, 1.30548200e-04, 1.30548200e-04,
    #                      1.30548200e-04, 1.30548200e-04, 1.42232126e-04, 7.19818715e-05]) / 10

    delta_rho = np.array([5.97088115e-05, 1.38820454e-04, 1.09198799e-04, 1.06191995e-04,
                          0.0022852, 0.0032757, 0.0067792, 0.0061167,
                          0.0036589, 0.0018909, 1.22843358e-04, 4.96285155e-05])/10

    delta_d = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,
               0.1, 0.1, 0.1, 0.1, 0.1]

    D2 = 0
    E2 = 0
    for i in range(len(rho)):
        D2 = D2 + rho[i] * d[i]*1e-8
        E2 = E2 + (delta_rho[i] * d[i]*1e-8) ** 2 + (rho[i] * delta_d[i]*1e-8) ** 2

    E2 = np.sqrt(E2) * (1e-8) ** 2 * (3.97) ** 2 * 6.02214076e23
    D2 = D2 * (1e-8) ** 2 * (3.97) ** 2 * 6.02214076e23

    print('Total Magnetism (10uc) = ' + str(D) + ' +- ' + str(E))
    x = [4, 7, 10]
    y = np.array([D/4,D1/7,D2/10])
    dy = [E/4,E1/7,E2/10]
    plt.figure()
    plt.errorbar(x, y, yerr=dy, fmt='o', capsize=3, linestyle='')
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.minorticks_on()
    plt.grid(True)
    #plt.ticklabel_format(axis='y', style='sci', scilimits=(-10, -10))
    plt.xlabel(r'Film Thickness (uc)', fontsize=13)
    plt.ylabel(r'Average magnetic contribution per unit cell  (arb. units)', fontsize=13)
    plt.ylim([0,0.125])
    plt.show()
    """

    """
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    thickness, density, density_mag = sample.density_profile()

    bounds = [[5.1822, 7.716215188017522], [7.716215188017522, 9.250230376035044], [9.250230376035044, 11.250230376035044],
              [11.250230376035044,17.066530376035044], [17.066530376035044, 22.882830376035044],
              [22.882830376035044,28.699130376035044], [28.699130376035044,7.82253]]
    value = np.array([2.47328055e-03, 1.05487732e-03, 1.49386942e-03, 1.25021663e-04,7.39237444e-05, 6.07739396e-05, 1.51467795e-05])/10

    tot_mag = sum(density_mag['Mn3+']/10)*0.1

    # find the error
    my_error = 0
    for i, bound in enumerate(bounds):

        if i == 0:
            idx = [j for j in range(len(thickness)) if thickness[j]<bound[1]]
            my_error = my_error + np.square(len(idx)*0.1*value[i])
        elif i == len(bounds)-1:
            idx = [j for j in range(len(thickness)) if thickness[j]>=bound[0]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])
        else:
            idx = [j for j in range(len(thickness)) if thickness[j]>=bound[0] and thickness[j]<bound[1]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])

    my_error = np.sqrt(my_error)
    print('10uc LSMO at 150K: ' + str(tot_mag/48.78) + ' +- ' + str(my_error/48.78))

    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_complete.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    thickness, density, density_mag = sample.density_profile()

    bounds = [[4.694941599778237, 7.88732020707415], [7.88732020707415, 9.093009125691446],
              [9.093009125691446,11.9407550580210903], [11.9407550580210903, 34.2809656880210903],
              [34.2809656880210903, 39.8660183460210903], [39.8660183460210903,47.5445163369150563]]
    value = np.array([3.28764707e-04, 3.65348818e-04, 9.56763350e-05, 5.22253313e-05,
                      4.86677504e-05, 1.02973070e-05]) / 10

    tot_mag = sum(density_mag['Mn3+'] / 10) * 0.1

    # find the error
    my_error = 0
    for i, bound in enumerate(bounds):

        if i == 0:
            idx = [j for j in range(len(thickness)) if thickness[j] < bound[1]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])
        elif i == len(bounds) - 1:
            idx = [j for j in range(len(thickness)) if thickness[j] >= bound[0]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])
        else:
            idx = [j for j in range(len(thickness)) if thickness[j] >= bound[0] and thickness[j] < bound[1]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])

    my_error = np.sqrt(my_error)
    print('13uc LSMO at 150K: ' + str(tot_mag/48.78) + ' +- ' + str(my_error/48.78))

    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-300K_complete.h5"

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    thickness, density, density_mag = sample.density_profile()

    bounds = [[4.694941599778237,7.88732020707415], [7.88732020707415, 9.093009125691445], [9.093009125691445,11.94075505802109],
              [11.94075505802109,24.28096568802109], [24.28096568802109,34.28096568802109], [34.28096568802109,39.86601834602109],
              [39.86601834602109,47.54451633691506]]
    value = np.array([4.26828299e-04, 4.86196647e-04, 1.35640551e-04, 5.84756968e-05, 6.12784553e-05,
                      5.40956013e-05, 1.08271152e-05]) / 10

    tot_mag = sum(density_mag['Mn3+'] / 10) * 0.1

    # find the error
    my_error = 0
    for i, bound in enumerate(bounds):

        if i == 0:
            idx = [j for j in range(len(thickness)) if thickness[j] < bound[1]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])
        elif i == len(bounds) - 1:
            idx = [j for j in range(len(thickness)) if thickness[j] >= bound[0]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])
        else:
            idx = [j for j in range(len(thickness)) if thickness[j] >= bound[0] and thickness[j] < bound[1]]
            my_error = my_error + np.square(len(idx) * 0.1 * value[i])

    my_error = np.sqrt(my_error)
    print('13uc LSMO at 300K: ' + str(tot_mag/48.78) + ' +- ' + str(my_error/48.78))
    """
    """
    my_d = [0,0,0,0,0,0,0,0]
    for i in range(1,9):
        idx = i - 1
        d = sample.getThickness(i,'all')
        if idx == 0:
            my_d[idx] = d
        else:
            my_d[idx] = my_d[idx-1] + d

    print(my_d)
    """
    """
    # --------------- ERROR FUNCTION FIT
    
    

    from scipy.optimize import curve_fit


    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_complete.h5"
    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')

    sample = ds.ReadSampleHDF5(fname)

    sample.energy_shift()

    xlim = [-10, 40]


    thickness, density, density_mag = sample.density_profile()

    idx = [i for i in range(len(thickness)) if thickness[i]>xlim[0] and thickness[i]<xlim[1]]



    #value = np.square(np.array([0.07291032,0.10234935,0.10512971,0.0725769]))*0.028  # 4uc Sr
    #bounds = [[0,3.905],[3.905,7.81],[7.81,11.75],[11.75,15.72]]  # 4uc Sr
    #value = np.square(np.array([0.06639064,0.07079927,0.06017485, 0.06487379]))*0.028  # 4uc Ti
    #bounds = [[3.905,7.81],[7.81,11.75],[11.75,15.72], [15.72,19.69]]  # 4uc Ti
    value = np.array([0.0108096,  0.0109932,  0.00976422, 0.01245532])  # 4uc mag
    bounds = [[3.905,7.81],[7.81,11.75],[11.75,15.72],[15.72,19.69]]  # 4uc Sr

    #value = np.square(np.array([0.06960216, 0.09067262, 0.09247652, 0.09119746, 0.08421491, 0.06104902])) * 0.028  # 7uc Sr
    #bounds = [[0, 3.905], [3.905, 7.81], [7.81, 11.75], [11.75, 15.72], [19.69,23.66], [23.66,27.63]]  # 7uc Sr
    #value = np.square(np.array([0.06887634, 0.05580875, 0.08318681, 0.0481571,  0.07007438])) * 0.028  # 7uc Ti
    #bounds = [[3.905, 7.81], [7.81, 11.75], [11.75, 15.72], [19.69,23.66], [23.66,27.63]]  # 7uc Ti
    #value = np.square(np.array([0.00397048, 0.00431724, 0.01569875, 0.01145907, 0.01279585]))  # 7uc mag
    #bounds = [[3.905, 7.81], [7.81, 11.75], [11.75, 15.72], [15.72, 19.69],[23.66,27.63]]  # 7uc Sr

    #value = np.array([0.06622369, 0.08949488, 0.09554424, 0.09757894, 0.09292542, 0.06574169]) * 0.028  # 10uc Sr
    #bounds = [[0, 3.905], [3.905, 7.81], [7.81, 11.75], [11.75, 15.72], [19.69, 23.66], [23.66, 27.63]]  # 10uc Sr
    #value = np.array([0.05235802, 0.05299326, 0.07487215, 0.05409448, 0.04625002]) * 0.028  # 10uc Ti
    #bounds = [ [3.905, 7.81], [7.81, 11.75], [11.75, 15.72], [19.69, 23.66], [23.66, 27.63]]  # 10uc Ti
    #value = np.square(np.array([0.00586707, 0.01305662, 0.01693715, 0.01284786, 0.01216948]))  # 10uc mag
    #bounds = [[3.905, 7.81], [7.81, 11.75], [11.75, 15.72], [15.72, 19.69], [23.66, 27.63]]  # 10uc Sr


    dy = createError(thickness, value, bounds)
    my_idx = [i for i in range(len(thickness)) if thickness[i]>19]
    density_mag['Mn3'][my_idx] = 0.0176497168459156
    popt, pcov = curve_fit(func, thickness[idx], density_mag['Mn3'][idx], sigma=dy[idx], absolute_sigma=True)
    perr = np.sqrt(np.diag(pcov))
    popt, pcov = curve_fit(func, thickness[idx], density_mag['Mn3'][idx])

    fitted = func(thickness, popt[0],popt[1], popt[2],popt[3])

    plt.figure()
    plt.plot(thickness, density_mag['Mn3'])
    plt.plot(thickness, fitted)
    plt.show()

    perr = np.sqrt(np.diag(pcov))
    print(popt)
    print(perr)
    """
    """
    data = np.loadtxt('35_460.76_S.txt')

    qz_data = data[0]
    R_data = data[1]

    qz_sim = data[2]
    R_sim = data[3]

    qz_nr = data[4]
    R_nr = data[5]

    qz_f = np.loadtxt('EM1-460-fourier.csv')[:, 4]
    R_f = np.loadtxt('EM1-460-fourier.csv')[:, 5]

    qz_s = np.loadtxt('EM1-460-spline.csv')[:, 4]
    R_s = np.loadtxt('EM1-460-spline.csv')[:, 5]
    up = 0.45
    fig, axes = plt.subplots(2,3)

    axes[0,0].plot(qz_data, R_data, marker='o', markersize=3, markerfacecolor='none', markeredgecolor='k', color='none')
    axes[0,0].plot(qz_sim, R_sim, 'b')
    axes[0,0].plot(qz_nr, R_nr, 'r')
    axes[0,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,0].tick_params(axis='y', labelleft=False)
    axes[0,0].legend(['Exp', 'Sim', 'NN'], loc='upper right', frameon=False)
    axes[0, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 0].set_xlim([axes[0, 0].get_xlim()[0], up])

    axes[0,1].plot(qz_data, R_data, marker='o', markersize=3, markerfacecolor='none', markeredgecolor='k', color='none')
    axes[0,1].plot(qz_s, R_s, 'r')
    axes[0,1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,1].tick_params(axis='y', labelleft=False)
    axes[0,1].legend(['Exp','Spline'], loc='upper right', frameon=False)
    axes[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].set_xlim([axes[0, 1].get_xlim()[0], up])

    axes[0,2].plot(qz_data, R_data, marker='o', markersize=3, markerfacecolor='none', markeredgecolor='k', color='none')
    axes[0,2].plot(qz_f, R_f, 'r')
    axes[0,2].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,2].tick_params(axis='y', labelleft=False)
    axes[0,2].legend(['Exp', 'Fourier'], loc='upper right', frameon=False)
    axes[0, 2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 2].set_xlim([axes[0, 2].get_xlim()[0], up])


    axes[1,0].plot(qz_data[:-1], np.diff(R_data)/np.diff(qz_data), 'k')
    axes[1,0].plot(qz_nr[:-1], np.diff(R_nr)/np.diff(qz_nr), 'r')
    axes[1,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1,0].tick_params(axis='y', labelleft=False)
    axes[1,0].legend(['Exp', 'NN'], loc='upper right', frameon=False)
    axes[1, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 0].set_xlim([axes[1, 0].get_xlim()[0], up])

    axes[1,1].plot(qz_data[:-1], np.diff(R_data) / np.diff(qz_data),'k')
    axes[1,1].plot(qz_s[:-1], np.diff(R_s) / np.diff(qz_s), 'r')
    axes[1,1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1,1].tick_params(axis='y', labelleft=False)
    axes[1,1].legend(['Exp', 'Spline'], loc='upper right', frameon=False)
    axes[1, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].set_xlim([axes[0, 0].get_xlim()[0], up])

    axes[1,2].plot(qz_data[:-1], np.diff(R_data) / np.diff(qz_data),'k')
    axes[1,2].plot(qz_f[:-1], np.diff(R_f) / np.diff(qz_f), 'r')
    axes[1,2].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1,2].tick_params(axis='y', labelleft=False)
    axes[1,2].legend(['Exp', 'Fourier'], loc='upper right', frameon=False)
    axes[1, 2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 2].set_xlim([axes[0, 0].get_xlim()[0], up])

    axes[1,1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA^{-1}}$)', fontsize=16)

    axes[0,0].set_ylabel(r'log $\left(\frac{R}{R_{0}} \right)$', fontsize=20)
    axes[1,0].set_ylabel(r'$\frac{d}{d q_{z}} \left[ \frac{R(q_{z})}{R_{0}(q_{z})}  \right]$', fontsize=20)

    plt.tight_layout()
    plt.show()
    """


    """
    fig, axes = plt.subplots(2,3)
    offset = -4.694941599778237
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_complete.h5"
    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()


    thickness, density, density_mag = sample.density_profile()
    thickness = thickness + offset

    Asites = ['Sr', 'La', 'O', 'C2']


    axes[0,1].plot(thickness, density['Sr'], color='cyan', label='Sr')
    axes[0,1].plot(thickness, density['La'], color='blue', label='La')
    axes[0,1].plot(thickness, density['O'], color='red', label='O')
    axes[0,1].plot(thickness, density['C2'], color='orange', label='C')

    # plt.fill_between(thickness, density['O'], color='red', alpha=0.001)
    # plt.fill_between(thickness, density['C2'], color='orange', alpha=0.25)
    axes[0,1].fill_between(thickness, density['Sr'], color='cyan', alpha=0.75)
    axes[0,1].fill_between(thickness, density['La'], color='blue', alpha=0.25)
    axes[0,1].set_xlim([axes[0,1].get_xlim()[0], 85])
    axes[0,1].axhline(0, color='black', linestyle='--')
    axes[0,1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].set_yticklabels([])
    axes[0, 1].set_xticklabels([])
    axes[0, 1].legend(frameon=False)


    Bsites = ['Ti', 'Mn2+', 'Mn3+']


    axes[1,1].plot(thickness, density['Ti'], color='yellow', label='Ti')
    axes[1,1].plot(thickness, density['Mn3+'], color='green', label=r'$\mathrm{Mn^{3.3+}}$')
    axes[1,1].plot(thickness, density['Mn2+'], color='grey', label=r'$\mathrm{Mn^{2+}}$')
    axes[1,1].plot(thickness, -density_mag['Mn3+'], color='magenta', label='Mag')

    axes[1,1].fill_between(thickness, density['Ti'], color='yellow', alpha=0.5)
    axes[1,1].fill_between(thickness, density['Mn3+'], color='green', alpha=0.25)
    axes[1,1].fill_between(thickness, density['Mn2+'], color='grey', alpha=1)
    axes[1,1].fill_between(thickness, -density_mag['Mn3+'], color='magenta', alpha=0.5)
    axes[1,1].set_xlim([axes[1,1].get_xlim()[0], 85])
    axes[1,1].set_ylim([-0.025, 0.032])
    axes[1,1].axhline(0, color='black', linestyle='--')
    axes[1,1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].set_yticklabels([])
    axes[1, 1].legend(frameon=False)


    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-300K_complete.h5"
    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    thickness, density, density_mag = sample.density_profile()
    thickness = thickness + offset

    Asites = ['Sr', 'La', 'O', 'C2']

    axes[0,2].plot(thickness, density['Sr'], color='cyan', label='Sr')
    axes[0,2].plot(thickness, density['La'], color='blue', label='La')
    axes[0,2].plot(thickness, density['O'], color='red', label='O')
    axes[0,2].plot(thickness, density['C2'], color='orange', label='C')

    # plt.fill_between(thickness, density['O'], color='red', alpha=0.001)
    # plt.fill_between(thickness, density['C2'], color='orange', alpha=0.25)
    axes[0,2].fill_between(thickness, density['Sr'], color='cyan', alpha=0.75)
    axes[0,2].fill_between(thickness, density['La'], color='blue', alpha=0.25)
    axes[0,2].set_xlim([axes[0,2].get_xlim()[0], 85])
    axes[0,2].axhline(0, color='black', linestyle='--')
    axes[0,2].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0,2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 2].set_yticklabels([])
    axes[0, 2].set_xticklabels([])
    axes[0, 2].legend(frameon=False)


    Bsites = ['Ti', 'Mn2+', 'Mn3+']

    axes[1,2].plot(thickness, density['Ti'], color='yellow', label='Ti')
    axes[1,2].plot(thickness, density['Mn3+'], color='green', label=r'$\mathrm{Mn^{3.3+}}$')
    axes[1,2].plot(thickness, density['Mn2+'], color='grey', label=r'$\mathrm{Mn^{2+}}$')
    axes[1,2].plot(thickness, -density_mag['Mn3+'], color='magenta', label='Mag')

    axes[1,2].fill_between(thickness, density['Ti'], color='yellow', alpha=0.5)
    axes[1,2].fill_between(thickness, density['Mn3+'], color='green', alpha=0.25)
    axes[1,2].fill_between(thickness, density['Mn2+'], color='grey', alpha=1)
    axes[1,2].fill_between(thickness, -density_mag['Mn3+'], color='magenta', alpha=0.5)
    axes[1,2].set_xlim([axes[1,2].get_xlim()[0], 85])
    axes[1,2].set_ylim([-0.025, 0.032])
    axes[1,2].axhline(0, color='black', linestyle='--')
    axes[1,2].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1, 2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 2].set_yticklabels([])
    axes[1, 2].legend(frameon=False)


    #shared_y = fig.add_subplot(111, frame_on=False)
    #shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_y.set_ylabel(r'Density, $\mathrm{\rho}$ ($\mathrm{mol/cm^{3}}$)', fontsize=12)

    #axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    #axes[0, 1].get_shared_y_axes().join(axes[1, 0], shared_y)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel(r'z Position $\mathrm{\AA}$)', fontsize=12)

    #axes[1, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[1, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[1, 0].get_shared_x_axes().join(axes[0, 0], shared_x)
    #axes[1, 1].get_shared_x_axes().join(axes[1, 0], shared_x)

    #plt.tight_layout()
    #plt.subplots_adjust(hspace=0)
    """
    """

    offset = -4.694941599778237
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v10.h5"
    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')

    sample = ds.ReadSampleHDF5(fname)
    sample.energy_shift()

    thickness, density, density_mag = sample.density_profile()
    thickness = thickness + offset

    Asites = ['Sr', 'La', 'O', 'C2']


    axes[0,0].plot(thickness, density['Sr'], color='cyan', label='Sr')
    axes[0,0].plot(thickness, density['La'], color='blue', label='La')
    axes[0,0].plot(thickness, density['O'], color='red', label='O')
    axes[0,0].plot(thickness, density['C2'], color='orange', label='C')

    #plt.fill_between(thickness, density['O'], color='red', alpha=0.001)
    #plt.fill_between(thickness, density['C2'], color='orange', alpha=0.25)
    axes[0,0].fill_between(thickness, density['Sr'], color='cyan', alpha=0.75)
    axes[0,0].fill_between(thickness, density['La'], color='blue', alpha=0.25)
    axes[0,0].set_xlim([axes[0,0].get_xlim()[0], 70])
    axes[0,0].axhline(0, color='black', linestyle='--')
    axes[0,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0,0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0,0].legend(frameon=False)
    axes[0,0].set_xticklabels([])



    Bsites = ['Ti','Mn2+', 'Mn3+']


    axes[1,0].plot(thickness, density['Ti'], color='yellow', label='Ti')
    axes[1,0].plot(thickness, density['Mn3+'], color='green', label=r'$\mathrm{Mn^{3.3+}}$')
    axes[1,0].plot(thickness, density['Mn2+'], color='grey', label=r'$\mathrm{Mn^{2+}}$')
    axes[1,0].plot(thickness, -density_mag['Mn3+'], color='magenta', label='Mag')

    axes[1,0].fill_between(thickness, density['Ti'], color='yellow', alpha=0.5)
    axes[1,0].fill_between(thickness, density['Mn3+'], color='green', alpha=0.25)
    axes[1,0].fill_between(thickness, density['Mn2+'], color='grey', alpha=1)
    axes[1,0].fill_between(thickness, -density_mag['Mn3+'], color='magenta', alpha=0.5)
    axes[1,0].set_xlim([axes[1,0].get_xlim()[0], 70])
    axes[1,0].set_ylim([-0.025,0.032])
    axes[1,0].axhline(0, color='black', linestyle='--')
    axes[1,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1,0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1,0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1,0].legend(frameon=False)

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel(r'Density, $\mathrm{\rho}$ ($\mathrm{mol/cm^{3}}$)', fontsize=16)

    axes[0,0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1,0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0,0].get_shared_y_axes().join(axes[0,0], shared_y)
    axes[1,0].get_shared_y_axes().join(axes[1,0], shared_y)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=12)

    #axes[0].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[1].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0].get_shared_x_axes().join(axes[0], shared_x)
    #axes[1].get_shared_x_axes().join(axes[1], shared_x)

    axes[1,1].set_xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)

    plt.show()
    """
    """
    dq_10 = 2.12

    dEz2 = 0.6*dq_10
    dEx2y2 = 0.6*dq_10
    dExzyz = -0.4*dq_10
    dExy = -0.4*dq_10
    import pickle
    orbitals = {'Ti4': [dExy,dExzyz,dEx2y2,dEz2]}

    #with open('Ti_orbitals.pkl', 'wb') as handle:
    #    pickle.dump(orbitals, handle)
    with open('Ti_orbitals.pkl', 'rb') as handle:
        b = pickle.load(handle)
    print(b)
    """
    """

    fig, axes = plt.subplots(2,3)
    # -------------------------------------------------- Atomic Slice Density Profile - 10uc
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_mag = sample.density_profile()

    # thickness_shift = -11.75
    # thickness = thickness + thickness_shift

    zero = np.zeros(len(thickness))



    xlim = [-30, 75]

    slice_width = [3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.94, 3.97, 3.97, 3.97, 3.97,3.97, 3.97, 3.97, 3.97,
                   3.97, 3.97]
    start = -15.62
    offset = -15.75
    bar_width = np.array(slice_width) - 0.3

    # properly create the half thickness array
    my_thickness = bar_graph_slice(start, slice_width)

    Ti = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
    Mn2 = [0, 0, 0, 0, 0, 0, 0.039197114403670075, 0.1019020202535671,0.19143669416656375,0.0923611842360831,0.09788910206548185,
           0, 0, 0, 0, 0, 0.14728965317765888, 0.3708915848554904]
    Mn3 = [0, 0, 0, 0, 0, 0, 1.2920867046687424e-05, 0.05604048341999848, 0.3462764023579752, 0.7950837175671438,
           0.8386105840135152, 1, 1, 1, 1, 1, 0.8527103468223411, 0.6291084151445097]
    mag = [0, 0, 0, 0, 0, 0, 0.00019341730101137163, 0.000186586051446197, 0.014676919509906398, 0.027492584861078195,
           0.031889765641913245, 0.031889765641913245, 0.031889765641913245, 0.031889765641913245, 0.031889765641913245,
           0.031889765641913245, 0.012130795599578477, 0.009068991034701678]

    mag = -np.array(mag) / max(mag)

    Brho = [0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,
            0.028, 0.028, 0.041592366909668055]
    Mn3 = np.array(Mn2) + np.array(Mn3)

    Ti = np.array(Ti) * np.array(Brho) / 0.028
    Mn2 = np.array(Mn2) * np.array(Brho) / 0.028
    Mn3 = np.array(Mn3) * np.array(Brho) / 0.028

    #plt.figure(1)

    axes[1,2].bar(np.array(my_thickness) + offset, Ti, color='yellow', edgecolor='black', width=bar_width, label='Ti')
    axes[1,2].bar(np.array(my_thickness) + offset, Mn3, color='grey', edgecolor='black', width=bar_width,
            label=r'$\mathrm{Mn^{3+}}$')
    axes[1,2].bar(np.array(my_thickness) + offset, Mn2, color='green', edgecolor='black', width=bar_width,
            label=r'$\mathrm{Mn^{2+}}$')
    axes[1,2].bar(np.array(my_thickness) + offset, mag, color='magenta', edgecolor='black', width=bar_width, label='Mag')
    axes[1,2].axhline(0, color='black', linestyle='--')
    axes[1,2].tick_params(axis='both', which='both', direction='in', top=True, right=True)

    #plt.gca().invert_xaxis()
    axes[1,2].set_xlim(xlim)
    axes[1,2].set_ylim(axes[1,2].set_ylim()[0], 1.95)
    axes[1, 2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 2].set_yticklabels([])
    axes[1,2].legend(frameon=False)
    # Asites
    Sr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0]
    La = [0, 0, 0, 0, 0, 0.05396744660948172, 0.22661826360730397, 0.36421004570555726, 0.6075721465851404, 0.8318896059187244,
          0.922032621683752, 1, 1, 1, 1, 1, 1, 1]

    Arho = [0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,0.028,
            0.028, 0.028, 0.0027315383881216246]

    Sr = np.array(Sr) * np.array(Arho) / 0.028
    La = np.array(La) * np.array(Arho) / 0.028
    O = density['O'] / 0.028
    C = density['C2'] / 0.028


    axes[0,2].bar(np.array(my_thickness) + offset, Sr, color='cyan', edgecolor='black', width=bar_width, label='Sr')
    axes[0,2].bar(np.array(my_thickness) + offset, La, color='blue', edgecolor='black', width=bar_width, label='La')
    axes[0,2].plot(thickness + offset+1.985+3.97, O, 'r', label='O')
    axes[0,2].plot(thickness + offset+1.985+3.97, C, color='orange', label='C')
    axes[0,2].axhline(0, color='black', linestyle='--')
    axes[0,2].set_xlim(xlim)

    #axes[0, 0].gca().set_xticklabels([])
    axes[0,2].set_ylim(-0.1, axes[0,2].get_ylim()[1])
    axes[0,2].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0, 2].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 2].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 2].set_xticklabels([])
    axes[0, 2].set_yticklabels([])
    axes[0,2].legend(frameon=False)
    # plt.gca().invert_xaxis()


    # -------------------------------------------------- Atomic Slice Density Profile - 7uc
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim7uc_unitCell_complete.h5"

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_mag = sample.density_profile()

    # thickness_shift = -11.75
    # thickness = thickness + thickness_shift


    xlim = [-30, 60]

    slice_width = [3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.905, 3.94, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97, 3.97]
    start = -15.62
    offset = -15.75
    bar_width = np.array(slice_width) - 0.3

    # properly create the half thickness array
    my_thickness = bar_graph_slice(start, slice_width)

    Ti = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    Mn2 = [0, 0, 0, 0, 0, 0, 0.12249797231667228, 0.13128036566477286, 0.18529187193983435, 0.05789249282762787, 0.0471180530123988,
           0, 0, 1.6187596931838133e-05, 0.41180846162665774]
    Mn3 = [0, 0, 0, 0, 0, 0, 0.0010238731064084816, 0.1465311417775903, 0.496282522173469, 0.9257043501763164,
           0.950174818457884, 1,1,0.9999838124030682, 0.5881915383733423]
    mag = [0, 0, 0, 0, 0, 0, 0.0003066750126945097,0.0002068731239481137,0.016932058948664633,0.02655782307499808,
           0.0314329491242263, 0.0314329491242263,0.0314329491242263,0.0314329491242263,0.00010936761100236586]

    mag = -np.array(mag) / max(mag)

    Brho = [0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,0.028, 0.028, 0.028,0.039184448609039656]
    Mn3 = np.array(Mn2) + np.array(Mn3)

    Ti = np.array(Ti) * np.array(Brho) / 0.028
    Mn2 = np.array(Mn2) * np.array(Brho) / 0.028
    Mn3 = np.array(Mn3) * np.array(Brho) / 0.028

    axes[1,1].bar(np.array(my_thickness) + offset, Ti, color='yellow', edgecolor='black', width=bar_width, label='Ti')
    axes[1,1].bar(np.array(my_thickness) + offset, Mn3, color='grey', edgecolor='black', width=bar_width,
            label=r'$\mathrm{Mn^{3+}}$')
    axes[1,1].bar(np.array(my_thickness) + offset, Mn2, color='green', edgecolor='black', width=bar_width,
            label=r'$\mathrm{Mn^{2+}}$')
    axes[1,1].bar(np.array(my_thickness) + offset, mag, color='magenta', edgecolor='black', width=bar_width, label='Mag')
    axes[1,1].axhline(0, color='black', linestyle='--')
    axes[1,1].tick_params(axis='both', which='both', direction='in', top=True, right=True)

    #plt.gca().invert_xaxis()
    axes[1,1].set_xlim(xlim)
    axes[1, 1].set_yticklabels([])
    axes[1,1].set_ylim(axes[1,1].get_ylim()[0], 1.95)
    axes[1, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1,1].legend(frameon=False)


    # Asites
    Sr = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0]
    La = [0, 0, 0, 0, 0, 0.09066289857636511, 0.22453131141735938, 0.4223460926698023, 0.7322631754977746, 0.9403263612357609, 0.9299965577372289, 1, 1, 1, 1]

    Arho = [0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028,0.028,0.028,0.028,0.0020259688666223116]

    Sr = np.array(Sr) * np.array(Arho) / 0.028
    La = np.array(La) * np.array(Arho) / 0.028
    O = density['O'] / 0.028
    C = density['C2'] / 0.028

    axes[0,1].bar(np.array(my_thickness) + offset, Sr, color='cyan', edgecolor='black', width=bar_width, label='Sr')
    axes[0,1].bar(np.array(my_thickness) + offset, La, color='blue', edgecolor='black', width=bar_width, label='La')
    axes[0,1].plot(thickness + offset +1.985+3.97, O, 'r', label='O')
    axes[0,1].plot(thickness + offset +1.985+3.97, C, color='orange', label='C')
    axes[0,1].axhline(0, color='black', linestyle='--')
    axes[0,1].set_xlim(xlim)

    #axes[1,0].gca().set_xticklabels([])
    axes[0,1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,1].set_ylim(-0.1, axes[0,1].get_ylim()[1])
    axes[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].set_xticklabels([])
    axes[0, 1].set_yticklabels([])
    axes[0,1].legend(frameon=False)
    # plt.gca().invert_xaxis()


    # -------------------------------------------------- Atomic Slice Density Profile - 4uc
    # fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v9.h5"
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_complete.h5"

    struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/')

    sample = ds.ReadSampleHDF5(fname)

    thickness, density, density_mag = sample.density_profile()

    #thickness_shift = -11.75
    #thickness = thickness + thickness_shift

    zero = np.zeros(len(thickness))

    xlim = [-30,45]

    slice_width = [3.905,3.905,3.905,3.905,3.905,3.905,3.905, 3.94, 3.97, 3.97,3.97,1.985]
    start = -15.62
    offset = -15.75
    bar_width = np.array(slice_width) - 0.3

    # properly create the half thickness array
    my_thickness = bar_graph_slice(start, slice_width)

    Ti  = [1,1,1,1,1, 1, 1, 1, 1, 1, 0, 0]
    Mn2 = [0,0,0,0,0, 0, 0, 0.011852930426155843, 0.2512293921053707, 0.0018988607009098994, 0.49522149100860147, 0.2120438416848024]
    Mn3 = [0,0,0,0,0, 0, 0.047918728033946256, 0.07279242172343539, 0.45243148969544544, 0.849697598400682, 0.5047785089913985, 0.7879561583151976]
    mag = [0,0,0,0,0, 0, 0.0001502703621304594, 0.0018479773003216926, 0.00712558371573116, 0.0176497168459156, 0.003261994045717825, 0]

    mag = -np.array(mag)/max(mag)

    Arho = [0.028,0.028,0.028,0.028,0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.028170951815703227, 0.049876937505378]
    Mn3 = np.array(Mn2) + np.array(Mn3)

    Ti = np.array(Ti)*np.array(Arho)/0.028
    Mn2 = np.array(Mn2) * np.array(Arho) / 0.028
    Mn3 = np.array(Mn3) * np.array(Arho) / 0.028

    axes[1,0].bar(np.array(my_thickness) + offset, Ti, color='yellow', edgecolor='black', width=bar_width, label='Ti')
    axes[1,0].bar(np.array(my_thickness) + offset, Mn3, color='grey',edgecolor='black', width=bar_width, label=r'$\mathrm{Mn^{3+}}$')
    axes[1,0].bar(np.array(my_thickness) + offset, Mn2, color='green', edgecolor='black', width=bar_width, label=r'$\mathrm{Mn^{2+}}$')
    axes[1,0].bar(np.array(my_thickness) + offset, mag, color='magenta', edgecolor='black', width=bar_width, label='Mag' )
    axes[1,0].axhline(0, color='black', linestyle='--')
    axes[1,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)

    #plt.gca().invert_xaxis()
    axes[1,0].set_xlim(xlim)
    axes[1,0].set_ylim(axes[1,0].get_ylim()[0], 1.95)
    axes[1, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1,0].legend(frameon=False)



    # Asites
    Sr = [1,1,1,1,1, 1, 1, 1, 1, 0, 0, 0]
    La = [0,0,0,0,0, 0.037082772361346605, 0.18186144221870582, 0.3851959154495215, 0.8012772709165596, 1, 1, 1]


    Brho = [0.028,0.028,0.028,0.028,0.028, 0.028, 0.028, 0.028, 0.028, 0.028, 0.025591699518419088, 0.00698397624316434]

    Sr = np.array(Sr) * np.array(Brho) / 0.028
    La = np.array(La) * np.array(Brho) / 0.028
    O = density['O']/0.028
    C = density['C2']/0.028

    axes[0,0].bar(np.array(my_thickness) +offset, Sr, color='cyan', edgecolor='black', width=bar_width, label='Sr')
    axes[0,0].bar(np.array(my_thickness) +offset, La, color='blue', edgecolor='black', width=bar_width, label='La')
    axes[0,0].plot(thickness+offset+1.985+3.97, O, 'r', label='O')
    axes[0,0].plot(thickness+offset+1.985+3.97, C, color='orange', label='C')
    axes[0,0].axhline(0, color='black', linestyle='--')
    axes[0,0].set_xlim(xlim)

    #axes[2,0].gca().set_xticklabels([])
    axes[0,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0,0].set_ylim(-0.1, axes[0,0].get_ylim()[1])
    axes[0, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())

    axes[0,0].legend(frameon=False)
    axes[0,0].set_xticklabels([])
    axes[1,1].set_xlabel(r'z Position ($\mathrm{\AA}$)', fontsize=20)
    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Fractional Site Occupation', fontsize=20)

    axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    axes[0, 1].get_shared_y_axes().join(axes[1, 0], shared_y)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    #axes[0,0].grid(True)
    #plt.gca().invert_xaxis()
    plt.show()

    """

    """

    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-300K_complete.h5"

    data1, data_dict1, sim_dict1 = ds.LoadDataHDF5(fname)
    #nR_num = ['26', '29', '31', '33', '35']
    nR_num = ['26', '29', '33']
    #R_num = ['39', '41', '22', '24', '37']
    R_num = ['41', '22', '24', '37']
    #Ti_num = ['30', '32', '34', '36']
    Ti_num = ['30', '34', '36']
    Mn_num = ['8', '15', '18', '23', '25']
    #La_num = ['40', '42']
    La_num = ['40']
    R_asymm_num = ['9-10', '13-14', '16-17']
    E_asymm_num = ['11-12', '20-21', '27-28']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant1 = return_name_from_scanNum(data1, nR_num)
    resonant1 = return_name_from_scanNum(data1, R_num)
    Ti_resonance1 = return_name_from_scanNum(data1, Ti_num)
    Mn_resonance1 = return_name_from_scanNum(data1, Mn_num)
    La_resonance1 = return_name_from_scanNum(data1, La_num)
    Mn_R_asymm1 = return_name_from_scanNum(data1, R_asymm_num)
    Mn_E_asymm1 = return_name_from_scanNum(data1, E_asymm_num)

    #plotting_scans(data_dict, sim_dict, non_resonant, offset=0, step=3, xmin=0.04, type='R', reverse=True, fig=1)
    #plotting_scans(data_dict, sim_dict, resonant, offset=0, step=3, xmin=0.04, type='R', reverse=True, fig=2)
    #plotting_scans(data_dict, sim_dict, Ti_resonance, offset=0, step=1, xmin=455, xmax=470, type='E', reverse=True,
    #               fig=3)
    #plotting_scans(data_dict, sim_dict, Mn_resonance, offset=0, step=1, xmin=635, xmax=660, type='E', reverse=True,
    #               fig=4)
    #plotting_scans(data_dict, sim_dict, La_resonance, offset=0, step=1, xmin=820, xmax=860, type='E', reverse=True,
    #               fig=5)
    #plotting_scans(data_dict, sim_dict, Mn_R_asymm, offset=0, step=1, xmin=0.035, type='RA', reverse=True, fig=6)
    #plotting_scans(data_dict, sim_dict, Mn_E_asymm, offset=0, step=2, xmin=635, xmax=660, type='EA', reverse=True,
    #               fig=7)


    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM2-150K_complete.h5"

    data2, data_dict2, sim_dict2 = ds.LoadDataHDF5(fname)
    #nR_num = ['64', '57', '69', '71', '73']
    nR_num = ['64', '57', '71']
    R_num = ['79', '60', '62', '75']
    # R_num = ['77', '79', '60', '62', '75']
    #Ti_num = ['68', '70', '72', '74']
    Ti_num = ['68', '72', '74']
    Mn_num = ['46', '53', '56', '61', '63']
    #La_num = ['78', '80']
    La_num = ['78']
    R_asymm_num = ['47-48', '51-52','54-55']
    E_asymm_num = ['49-50', '58-59', '65-66']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant2 = return_name_from_scanNum(data2, nR_num)
    resonant2 = return_name_from_scanNum(data2, R_num)
    Ti_resonance2 = return_name_from_scanNum(data2, Ti_num)
    Mn_resonance2 = return_name_from_scanNum(data2, Mn_num)
    La_resonance2 = return_name_from_scanNum(data2, La_num)
    Mn_R_asymm2 = return_name_from_scanNum(data2, R_asymm_num)
    Mn_E_asymm2 = return_name_from_scanNum(data2, E_asymm_num)

    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v10.h5"
    # struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    data, data_dict, sim_dict = ds.LoadDataHDF5(fname)
    R_num = ['35' , '22', '24', '33']
    nR_num = ['26', '19', '31']
    Ti_num = ['29', '30', '32']
    Mn_num = ['8', '15', '18', '23', '25']
    La_num = ['34']
    R_asymm_num = ['9-10', '13-14', '16-17']
    E_asymm_num = ['11-12', '20-21', '27-28']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    resonant = return_name_from_scanNum(data, R_num)
    non_resonant = return_name_from_scanNum(data, nR_num)
    Ti_resonance = return_name_from_scanNum(data, Ti_num)
    Mn_resonance = return_name_from_scanNum(data, Mn_num)
    La_resonance = return_name_from_scanNum(data, La_num)
    Mn_R_asymm = return_name_from_scanNum(data, R_asymm_num)
    Mn_E_asymm = return_name_from_scanNum(data, E_asymm_num)

    #plotting_scans(data_dict, sim_dict, non_resonant, offset=0, step=3, xmin=0.04, type='R', reverse=True, fig=8)
    #plotting_scans(data_dict, sim_dict, resonant, offset=0, step=3, xmin=0.04, type='R', reverse=True, fig=9)
    #plotting_scans(data_dict, sim_dict, Ti_resonance, offset=0, step=1, xmin=455, xmax=470, type='E', reverse=True,
    #               fig=10)
    #plotting_scans(data_dict, sim_dict, Mn_resonance, offset=0, step=1, xmin=635, xmax=660, type='E', reverse=True,
    #               fig=11)
    #plotting_scans(data_dict, sim_dict, La_resonance, offset=0, step=1, xmin=820, xmax=860, type='E', reverse=True,
    #               fig=12)
    #plotting_scans(data_dict, sim_dict, Mn_R_asymm, offset=0, step=1, xmin=0.035, type='RA', reverse=True, fig=13)
    #plotting_scans(data_dict, sim_dict, Mn_E_asymm, offset=0, step=2, xmin=635, xmax=660, type='EA', reverse=True,
    #               fig=14)

    fig, axes = plt.subplots(2,3)
    plotting_scans(data_dict, sim_dict, Mn_resonance, axes, (0, 0), offset=0, step=1.5, xmin=635, xmax=660,
                   type='E', reverse=True, size=1, top=1.5)
    plotting_scans(data_dict, sim_dict, Mn_E_asymm, axes, (1, 0), offset=0, step=3, xmin=635,
                   type='EA', reverse=True, size=1, top=0.75)
    plotting_scans(data_dict2, sim_dict2, Mn_resonance2, axes, (0, 1), offset=0, step=1.5, xmin=635, xmax=660,
                   type='E', reverse=True, size=1, top=1.5)
    plotting_scans(data_dict2, sim_dict2, Mn_E_asymm2, axes, (1, 1), offset=0, step=3, xmin=635,
                   type='EA', reverse=True, size=1, top=0.75)
    plotting_scans(data_dict1, sim_dict1, Mn_resonance1, axes, (0, 2), offset=0, step=1.5, xmin=635, xmax=660,
                   type='E', reverse=True, size=1, top=1.5)
    plotting_scans(data_dict1, sim_dict1, Mn_E_asymm1, axes, (1, 2), offset=0, step=3, xmin=635,
                   type='EA', reverse=True, size=1, top=0.75)


    #shared_y = fig.add_subplot(111, frame_on=False)
    #shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_y.set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)

    #axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    #axes[0, 1].get_shared_y_axes().join(axes[1, 0], shared_y)
    # axes[0, 1].get_shared_x_axes().join(axes[0, 2], shared_x)


    #axes[1, 0].set_xlabel('Energy, E (eV)', fontsize=12)
    #axes[1,1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA^{-1}}$)', fontsize=12)
    #axes[0,0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    axes[1,1].set_xlabel(r'Energy, E (eV)', fontsize=16)
    axes[0,0].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=13)
    axes[1,0].set_ylabel(r'Circular Polarization Asymmetry (arb. units)', fontsize=13)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    #plt.subplots_adjust(hspace=0.08)
    #shared_y.yaxis.set_label_coords(-0.01, 0.5)
    plt.show()

    plt.show()

    fig, axes = plt.subplots(1, 3)

    plotting_scans(data_dict, sim_dict, Mn_R_asymm, axes, (0, 0), offset=0, step=0.75, xmin=0.04, xmax=0.5,
                   type='RA', reverse=True, top=0.5)
    plotting_scans(data_dict2, sim_dict2, Mn_R_asymm2, axes, (0, 1), offset=0, step=0.75, xmin=0.04,xmax=0.5,
                   type='RA', reverse=True, top=0.5)
    plotting_scans(data_dict1, sim_dict1, Mn_R_asymm1, axes, (0, 2), offset=0, step=0.75, xmin=0.04, xmax=0.5,
                   type='RA', reverse=True, top=0.5)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_ylabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)', fontsize=12)

    #axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    #axes[0, 1].get_shared_y_axes().join(axes[0, 1], shared_y)
    #axes[1, 0].set_xlabel('Energy, E (eV)', fontsize=12)
    #axes[1, 1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA^{-1}}$)', fontsize=12)
    # axes[0,0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)

    # shared_ax.xaxis.set_label_coords(0.35, -0.075)

    axes[1].set_xlabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)', fontsize=16)
    axes[0].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.08)
    plt.show()

    # supplemental energy scans!
    fig, axes = plt.subplots(2, 3)
    plotting_scans(data_dict, sim_dict, Ti_resonance, axes, (0, 0), offset=0, step=1, xmin=455, xmax=470,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict2, sim_dict2, Ti_resonance2, axes, (0, 1), offset=0, step=1, xmin=455, xmax=470,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict1, sim_dict1, Ti_resonance1, axes, (0, 2), offset=0, step=1, xmin=455, xmax=470,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict, sim_dict, La_resonance, axes, (1, 0), offset=0, step=1, xmin=830, xmax=860,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict2, sim_dict2, La_resonance2, axes, (1, 1), offset=0, step=1, xmin=830, xmax=860,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict1, sim_dict1, La_resonance1, axes, (1, 2), offset=0, step=1, xmin=830, xmax=860,
                   type='E', reverse=True, size=1, top=1.25)

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)

    axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    axes[0, 1].get_shared_y_axes().join(axes[0, 1], shared_y)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel('Energy, E (eV)', fontsize=12)

    #axes[0, 1].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[1, 1].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 1].get_shared_x_axes().join(axes[0, 1], shared_x)
    #axes[1, 1].get_shared_x_axes().join(axes[1, 0], shared_x)

    axes[1, 1].set_xlabel(r'Energy, E (eV)', fontsize=16)
    # axes[0,0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)

    #x_label = axes[1,2].xaxis.get_label()
    #x_label_pos = x_label.get_position()

    #shared_x.xaxis.set_label_coords(0.35, -0.04)
    shared_y.yaxis.set_label_coords(-0.02, 0.5)

    plt.tight_layout()
    #plt.subplots_adjust(hspace=0.08)
    plt.show()

    fig, axes = plt.subplots(2,3)

    plotting_scans(data_dict, sim_dict, resonant, axes, (0, 0), offset=0, step=3, xmin=0.01, xmax=0.9,
                   type='R', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict2, sim_dict2, resonant2, axes, (0, 1), offset=0, step=3, xmin=0.037, xmax=0.9,
                   type='R', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict1, sim_dict1,resonant1, axes, (0, 2), offset=0, step=3, xmin=0.037, xmax=0.9,
                   type='R', reverse=True, size=1, top=1.25)

    plotting_scans(data_dict, sim_dict, non_resonant, axes, (1, 0), offset=0, step=3, xmin=0.01, xmax=0.9,
                   type='R', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict2, sim_dict2, non_resonant2, axes, (1, 1), offset=0, step=3, xmin=0.037, xmax=0.9,
                   type='R', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict1, sim_dict1, non_resonant1, axes, (1, 2), offset=0, step=3, xmin=0.037, xmax=0.9,
                   type='R', reverse=True, size=1, top=1.25)

    axes[0,0].set_xticks([0.1,0.2,0.3,0.4,0.5,0.6,0.7])
    axes[0, 1].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
    axes[0, 2].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])

    axes[1, 0].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5])
    axes[1, 1].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5])
    axes[1, 2].set_xticks([0.1, 0.2, 0.3, 0.4, 0.5])

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)

    axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    axes[1, 0].get_shared_y_axes().join(axes[1, 0], shared_y)

    axes[1,1].set_xlabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)', fontsize=16)
    shared_y.yaxis.set_label_coords(-0.02, 0.5)
    plt.tight_layout()
    plt.show()
    
    
    fname = "//cabinet/work$/lsk601/My Documents/LSMO_For_Lucas/RXR_Twente-EM1-150K_v10.h5"
    #struct_names, mag_names = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    data, data_dict, sim_dict = ds.LoadDataHDF5(fname)
    R_num = ['26', '35', '19', '31', '22', '24','33']
    Ti_num = ['29', '30', '32']
    Mn_num = ['8', '15', '18', '23', '25']
    La_num = ['34']
    R_asymm_num = ['9-10', '13-14','16-17']
    E_asymm_num = ['11-12', '20-21', '27-28']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    resonant = return_name_from_scanNum(data, R_num)
    Ti_resonance = return_name_from_scanNum(data, Ti_num)
    Mn_resonance = return_name_from_scanNum(data, Mn_num)
    La_resonance = return_name_from_scanNum(data, La_num)
    Mn_R_asymm = return_name_from_scanNum(data, R_asymm_num)
    Mn_E_asymm = return_name_from_scanNum(data, E_asymm_num)

    #plotting_scans(data_dict, sim_dict, non_resonant, offset=0, step=3, xmin=0.01, type='R', reverse=True, fig=1)
    #plotting_scans(data_dict, sim_dict, resonant, offset=0, step=3, xmin=0.005, type='R', reverse=True, fig=1)
    #plotting_scans(data_dict, sim_dict, Ti_resonance, offset=0, step=1, xmin=455, xmax=470, type='E', reverse=True,
    #               fig=2)
    #plotting_scans(data_dict, sim_dict, Mn_resonance, offset=0, step=1, xmin=635, xmax=660, type='E', reverse=True,
    #               fig=3)
    #plotting_scans(data_dict, sim_dict, La_resonance, offset=0, step=0.75, xmin=820, xmax=860, type='E', reverse=True,
    #               fig=4)
    #plotting_scans(data_dict, sim_dict, Mn_R_asymm, offset=0, step=1, xmin=0.01, type='RA', reverse=True, fig=6)
    #plotting_scans(data_dict, sim_dict, Mn_E_asymm, offset=0, step=2, xmin=635, xmax=660, type='EA', reverse=True,
    #               fig=5)
    #plt.show()
    from matplotlib import gridspec

   
    

    # shared_x = fig.add_subplot(111, frame_on=False)
    # shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    # shared_x.set_xlabel('Energy, E (eV)')

    # axes[0, 0].set_xlabel('')  # Remove existing x-axis label for this subplot
    # axes[0, 1].set_xlabel('')  # Remove existing x-axis label for this subplot
    # axes[0, 0].get_shared_x_axes().join(axes[0, 0], shared_x)
    # axes[0, 1].get_shared_x_axes().join(axes[0, 1], shared_x)
    # axes[0, 1].get_shared_x_axes().join(axes[0, 2], shared_x)

    """
    """
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim10uc_unitCell_complete.h5"
    struct_names_10, mag_names_10 = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    data_10, data_dict_10, sim_dict_10 = ds.LoadDataHDF5(fname)
    nR_num = ['9', '14', '19', '24', '29', '34']
    R_num = ['39', '45', '51', '57', '63', '64']
    #Ti_num = ['59', '47', '35', '25', '15', '10', '20', '30', '41', '53', '65']
    Ti_num = ['47', '15']
    #Mn_num = ['61', '49', '37', '27', '17', '12', '22', '32', '43', '55', '67']
    Mn_num = ['61', '37', '17', '22', '43', '67']
    La_num = ['69', '70']
    R_asymm_num = ['111-112', '117-118']
    E_asymm_num = ['127-128', '125-126', '121-122', '115-116', '109-110']
    #E_asymm_num = ['127-128', '125-126', '121-122', '115-116', '109-110', '107-108', '113-114', '119-120', '123-124']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant_10 = return_name_from_scanNum(data_10, nR_num)
    resonant_10 = return_name_from_scanNum(data_10, R_num)
    Ti_resonance_10 = return_name_from_scanNum(data_10, Ti_num)
    Mn_resonance_10 = return_name_from_scanNum(data_10, Mn_num)
    La_resonance_10 = return_name_from_scanNum(data_10, La_num)
    Mn_R_asymm_10 = return_name_from_scanNum(data_10, R_asymm_num)
    Mn_E_asymm_10 = return_name_from_scanNum(data_10, E_asymm_num)

    #plotting_scans(data_dict, sim_dict, non_resonant, offset=0, step=3, xmin=0.01, type='R', reverse=True, fig=1)
    #plotting_scans(data_dict, sim_dict, resonant, offset=0, step=3, xmin=0.005, type='R', reverse=True, fig=2)
    #plotting_scans(data_dict, sim_dict, Ti_resonance, offset=0, step=1, xmin=455, xmax=470, type='E', reverse=True,
    #               fig=3)
    #plotting_scans(data_dict, sim_dict, Mn_resonance, offset=0, step=1, xmin=635, xmax=660, type='E', reverse=True,
    #               fig=4)
    #plotting_scans(data_dict, sim_dict, La_resonance, offset=0, step=1, xmin=830, xmax=860, type='E', reverse=True,
    #               fig=5)
    #plotting_scans(data_dict, sim_dict, Mn_R_asymm, offset=0, step=0.5, xmin=0.01, type='RA', reverse=True, fig=6)
    #plotting_scans(data_dict, sim_dict, Mn_E_asymm, offset=0, step=2, xmin=635, xmax=660, type='EA', reverse=True,
    #               fig=7)
    #plt.show()

    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim7uc_unitCell_complete.h5"
    struct_names_7, mag_names_7 = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    data_7, data_dict_7, sim_dict_7 = ds.LoadDataHDF5(fname)
    nR_num = ['8','13','18','23','28','33']
    R_num = ['38','42','48','52','58','59']
    #Ti_num = ['54','44','34','24','14','9','19','29','40','50','60']
    Ti_num = ['44', '14']
    #Mn_num = ['56','46','36','26','16','11','21','31','41','51','61']
    Mn_num = ['56', '36', '16', '21', '41', '61']
    La_num = ['62','63']
    R_asymm_num = ['104-105','110-111']
    E_asymm_num = ['120-121','118-119','114-115','108-109','102-103']
    #E_asymm_num = ['120-121','118-119','114-115','108-109','102-103','100-101','106-107','112-113','116-117']

    # non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant_7 = return_name_from_scanNum(data_7, nR_num)
    resonant_7 = return_name_from_scanNum(data_7, R_num)
    Ti_resonance_7 = return_name_from_scanNum(data_7, Ti_num)
    Mn_resonance_7 = return_name_from_scanNum(data_7, Mn_num)
    La_resonance_7 = return_name_from_scanNum(data_7, La_num)
    Mn_R_asymm_7 = return_name_from_scanNum(data_7, R_asymm_num)
    Mn_E_asymm_7 = return_name_from_scanNum(data_7, E_asymm_num)

    #plotting_scans(data_dict, sim_dict, non_resonant, offset=0, step=3, xmin=0.01, type='R', reverse=True, fig=8)
    #plotting_scans(data_dict, sim_dict, resonant, offset=0, step=3, xmin=0.005, type='R', reverse=True, fig=9)
    #plotting_scans(data_dict, sim_dict, Ti_resonance, offset=0, step=1, xmin=455, xmax=470, type='E', reverse=True,
    #               fig=10)
    #plotting_scans(data_dict, sim_dict, Mn_resonance, offset=0, step=1, xmin=635, xmax=660, type='E', reverse=True,
    #               fig=11)
    #plotting_scans(data_dict, sim_dict, La_resonance, offset=0, step=1, xmin=830, xmax=860, type='E', reverse=True,
    #               fig=12)
    #plotting_scans(data_dict, sim_dict, Mn_R_asymm, offset=0, step=0.5, xmin=0.01, type='RA', reverse=True, fig=13)
    #plotting_scans(data_dict, sim_dict, Mn_E_asymm, offset=0, step=2, xmin=635, xmax=660, type='EA', reverse=True,
    #               fig=14)
    #plt.show()


    # plots for the 4uc sample
    fname = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_complete.h5"
    struct_names_4, mag_names_4 = mm._use_given_ff('//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/')
    data_4, data_dict_4, sim_dict_4 = ds.LoadDataHDF5(fname)

    #non_resonant = ['9_399.39_S','14_499.93_S','19_600.18_S','24_700.14_S','29_799.82_S','34_899.22_S']
    non_resonant_4 = ['9_399.39_S', '14_499.93_S', '19_600.18_S', '24_700.14_S', '29_799.82_S']
    resonant_4 = ['39_455.73_S', '45_459.75_S','51_640.2_S','57_642.2_S','63_833.85_S','64_836.04_S']
    Ti_resonance_4 = ['47_E429.58_Th10.0_S','15_E429.58_Th25.0_S']
    # Ti_resonance = ['59_E429.58_Th5.0_S', '47_E429.58_Th10.0_S', '35_E429.58_Th15.0_S', '25_E429.58_Th20.0_S',
    #                     '15_E429.58_Th25.0_S', '10_E429.58_Th30.0_S', '20_E429.58_Th35.0_S', '30_E429.58_Th40.0_S',
    #                     '41_E429.58_Th45.0_S', '53_E429.58_Th50.0_S', '65_E429.58_Th55.0_S']
    #Mn_resonance = ['61_E600.18_Th5.0_S', '49_E600.18_Th10.0_S', '37_E600.18_Th15.0_S', '27_E600.18_Th20.0_S',
    #                '17_E600.18_Th25.0_S', '12_E600.18_Th30.0_S', '22_E600.18_Th35.0_S', '32_E600.18_Th40.0_S',
    #                '43_E600.18_Th45.0_S', '55_E600.18_Th50.0_S', '67_E600.18_Th55.0_S']
    Mn_resonance_4 = ['61_E600.18_Th5.0_S', '37_E600.18_Th15.0_S', '17_E600.18_Th25.0_S', '22_E600.18_Th35.0_S',
                    '43_E600.18_Th45.0_S', '67_E600.18_Th55.0_S']

    La_resonance_4 = ['69_E814.75_Th10.0_S', '70_E814.75_Th25.0_S']

    Mn_R_asymm_4 = ['112-113_640.2_AC_Asymm', '118-119_642.2_AC_Asymm']
    Mn_E_asymm_4 = ['128-129_E600.18_Th5.0_AC_Asymm', '126-127_E600.18_Th10.0_AC_Asymm', '122-123_E600.18_Th15.0_AC_Asymm',
                  '116-117_E600.18_Th20.0_AC_Asymm', '110-111_E600.18_Th25.0_AC_Asymm']
    # Mn_E_asymm = ['128-129_E600.18_Th5.0_AC_Asymm', '126-127_E600.18_Th10.0_AC_Asymm', '122-123_E600.18_Th15.0_AC_Asymm',
    #                   '116-117_E600.18_Th20.0_AC_Asymm', '110-111_E600.18_Th25.0_AC_Asymm', '108-109_E600.18_Th30.0_AC_Asymm',
    #                   '114-115_E600.18_Th35.0_AC_Asymm', '120-121_E600.18_Th40.0_AC_Asymm', '124-125_E600.18_Th45.0_AC_Asymm']

    #plotting_scans(data_dict, sim_dict, non_resonant, offset=0, step=3, xmin=0.01, type='R', reverse=True, fig=15)
    #plotting_scans(data_dict, sim_dict, resonant, offset=0, step=3, xmin=0.01, type='R', reverse=True, fig=16)
    #plotting_scans(data_dict, sim_dict, Ti_resonance, offset=0, step=1, xmin=455, xmax=470, type='E', reverse=True, fig=17)
    #plotting_scans(data_dict, sim_dict, Mn_resonance, offset=0, step=1,xmin=635,xmax=660, type='E', reverse=True, fig=18)
    #plotting_scans(data_dict, sim_dict, La_resonance, offset=0, step=1, xmin=830, xmax=860, type='E', reverse=True, fig=19)
    #plotting_scans(data_dict, sim_dict, Mn_R_asymm, offset=0, step=0.5, xmin=0.01, type='RA', reverse=True,fig=20)
    #plotting_scans(data_dict, sim_dict, Mn_E_asymm, offset=0, step=2, xmin=635, xmax=660, type='EA', reverse=True,fig=21)
    #plt.show()

    fig, axes = plt.subplots(2,3)

    plotting_scans(data_dict_4, sim_dict_4, Mn_resonance_4,axes,(0,0), offset=0, step=1,xmin=635,xmax=660,
                   type='E', reverse=True, size=1, top=2.3)
    plotting_scans(data_dict_7, sim_dict_7, Mn_resonance_7, axes, (0, 1), offset=0, step=1, xmin=635, xmax=660,
                   type='E', reverse=True, size=1, top=2.3)
    plotting_scans(data_dict_10, sim_dict_10, Mn_resonance_10, axes, (0, 2), offset=0, step=1, xmin=635, xmax=660,
                   type='E', reverse=True, size=1, top=2.3)
    plotting_scans(data_dict_4, sim_dict_4, Mn_E_asymm_4, axes, (1, 0), offset=0, step=2, xmin=635, xmax=660,
                   type='EA', reverse=True, size=1)
    plotting_scans(data_dict_7, sim_dict_7, Mn_E_asymm_7, axes, (1, 1), offset=0, step=2, xmin=635, xmax=660,
                   type='EA', reverse=True, size=1)
    plotting_scans(data_dict_10, sim_dict_10, Mn_E_asymm_10, axes, (1, 2), offset=0, step=2, xmin=635, xmax=660,
                   type='EA', reverse=True, size=1)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel('Energy, E (eV)')


    #axes[0, 0].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 1].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[0, 0].get_shared_x_axes().join(axes[0, 0], shared_x)
    #axes[0, 1].get_shared_x_axes().join(axes[0, 1], shared_x)
    #axes[0, 1].get_shared_x_axes().join(axes[0, 2], shared_x)

    axes[0,0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=13)
    axes[1, 0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=13)
    axes[1,1].set_xlabel('Energy, E (eV)', fontsize=16)
    #shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    plt.show()

    fig, axes = plt.subplots(1,3)

    plotting_scans(data_dict_4, sim_dict_4, Mn_R_asymm_4, axes, (0, 0), offset=0, step=0.5, xmin=0.01,
                   type='RA', reverse=True, top=0)
    plotting_scans(data_dict_7, sim_dict_7, Mn_R_asymm_7, axes, (0, 1), offset=0, step=0.5, xmin=0.01,
                   type='RA', reverse=True, top=-0.25)
    plotting_scans(data_dict_10, sim_dict_10, Mn_R_asymm_10, axes, (0, 2), offset=0, step=0.5, xmin=0.01,
                   type='RA', reverse=True, top=-0.25)

    axes[0].set_ylabel('Circular Polarized Asymmetry (arb. units)', fontsize=14)
    axes[1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA^{-1}}$)', fontsize=14)
    plt.tight_layout()
    plt.show()

    fig, axes = plt.subplots(2, 3)

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)

    my_data = data_dict_7['38_455.73_S']['Data']
    my_idx = [int(i) for i in range(len(my_data[0])) if my_data[0][i] < 0.37]

    data_dict_7['38_455.73_S']['Data'][0] = data_dict_7['38_455.73_S']['Data'][0][my_idx]
    data_dict_7['38_455.73_S']['Data'][2] = data_dict_7['38_455.73_S']['Data'][2][my_idx]
    sim_dict_7['38_455.73_S']['Data'][0] = sim_dict_7['38_455.73_S']['Data'][0][my_idx]
    sim_dict_7['38_455.73_S']['Data'][2] = sim_dict_7['38_455.73_S']['Data'][2][my_idx]

    plotting_scans(data_dict_4, sim_dict_4, resonant_4, axes, (0, 0), offset=0, step=3, xmin=0.01,
                   type='R', reverse=True, size=1, top=0)
    plotting_scans(data_dict_7, sim_dict_7, resonant_7, axes, (0, 1), offset=0, step=3, xmin=0.01,
                   type='R', reverse=True, size=1, top=0)
    plotting_scans(data_dict_10, sim_dict_10, resonant_10, axes, (0, 2), offset=0, step=3, xmin=0.01,
                   type='R', reverse=True, size=1, top=0)
    plotting_scans(data_dict_4, sim_dict_4, non_resonant_4, axes, (1, 0), offset=0, step=3, xmin=0.01,
                   type='R', reverse=True, size=1, top=0)
    plotting_scans(data_dict_7, sim_dict_7, non_resonant_7, axes, (1, 1), offset=0, step=3, xmin=0.01,
                   type='R', reverse=True, size=1, top=0)
    plotting_scans(data_dict_10, sim_dict_10, non_resonant_10, axes, (1, 2), offset=0, step=3, xmin=0.01,
                   type='R', reverse=True, size=1, top=0)

    axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    axes[0, 1].get_shared_y_axes().join(axes[1, 0], shared_y)

    #axes[0, 0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    #axes[1, 0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=12)
    axes[1, 1].set_xlabel(r'Momentum Transfer, $\mathrm{q_{z}}$ ($\mathrm{\AA^{-1}}$)', fontsize=16)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    shared_y.yaxis.set_label_coords(-0.01, 0.5)
    plt.show()

    fig, axes = plt.subplots(2, 3)

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=16)

    plotting_scans(data_dict_4, sim_dict_4, Ti_resonance_4, axes, (0, 0), offset=0, step=1, xmin=455, xmax=470,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict_7, sim_dict_7, Ti_resonance_7, axes, (0, 1), offset=0, step=1, xmin=455, xmax=470,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict_10, sim_dict_10, Ti_resonance_10, axes, (0, 2), offset=0, step=1, xmin=455, xmax=470,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict_4, sim_dict_4, La_resonance_4, axes, (1, 0), offset=0, step=1, xmin=830, xmax=860,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict_7, sim_dict_7, La_resonance_7, axes, (1, 1), offset=0, step=1, xmin=830, xmax=860,
                   type='E', reverse=True, size=1, top=1.25)
    plotting_scans(data_dict_10, sim_dict_10, La_resonance_10, axes, (1, 2), offset=0, step=1, xmin=830, xmax=860,
                   type='E', reverse=True, size=1, top=1.25)

    axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_y)
    axes[0, 1].get_shared_y_axes().join(axes[1, 0], shared_y)

    # axes[0, 0].set_ylabel('Normalized Reflected Intensity (arb. units)', fontsize=12)
    # axes[1, 0].set_ylabel('Circular Polarization Asymmetry (arb. units)', fontsize=12)
    axes[1, 1].set_xlabel('Energy, E (eV)', fontsize=16)
    # shared_ax.xaxis.set_label_coords(0.35, -0.075)
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.075)
    shared_y.yaxis.set_label_coords(-0.01, 0.5)
    plt.show()
    """
    """
    name1 = '47_E429.58_Th10.0_S'

    E = data_dict[name1]['Data'][3]
    R = data_dict[name1]['Data'][2]

    plt.figure()
    plt.plot(E,R)
    plt.xlabel('Energy, E (eV)')
    plt.ylabel('Reflectivity (arb. units)')

    name2 = '27_E600.18_Th20.0_S'
    Theta = data_dict[name2]['Angle']
    E = data_dict[name2]['Data'][3]
    R = data_dict[name2]['Data'][2]

    #sample.eShift['Mn2'] = 0
    #sample.eShift['Mn3'] = 0

    E, Rnew = sample.energy_scan(Theta, E)

    Rnew = Rnew['S']


    plt.figure(2)
    plt.plot(E, R)
    plt.plot(E, Rnew)
    plt.xlabel('Energy, E (eV)')
    plt.ylabel('Reflectivity (arb. units)')
    plt.legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'])
    plt.show()
    #print(sample.eShift)
    """
    """
    data = np.loadtxt('31_635.99_S.txt')

    plt.figure(1)
    plt.plot(data[0],data[1])
    plt.plot(data[2],data[3])
    plt.plot(data[4],data[5])

    plt.figure(2)
    plt.plot(np.diff(data[1]))
    plt.plot(np.diff(data[5]))
    plt.show()
    """
    """

    # Total Variation

    fname1 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_test1-0.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_test1-100.h5"
    fname3 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_test1-1000000.h5"

    data1, data_dict1, sim_dict1 = ds.ReadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.ReadDataHDF5(fname2)
    data3, data_dict3, sim_dict3 = ds.ReadDataHDF5(fname3)

    keys = ["64_836.04_S", "63_833.85_S"]
    #keys = ["63_833.85_S"]
    fig, axes = plt.subplots(1,2)

    for i, key in enumerate(keys):
        E = data_dict1[key]['Energy']
        qz = data_dict1[key]['Data'][0]
        R = np.log10(data_dict1[key]['Data'][2])

        Rs = np.log10(sim_dict1[key]['Data'][2])

        l1_norm = sum(np.abs(R - Rs)) / len(R)
        tv = total_variation(R, Rs)
        idx = [i for i in range(len(qz)) if qz[i] >= 0.01]

        axes[0].plot(qz[idx], R[idx] + i * 2, 'b')
        axes[0].plot(qz[idx], Rs[idx] + i * 2, 'r')

        print('E=' + str(E) + ': ' + str(l1_norm))

    for i, key in enumerate(keys):
        E = data_dict2[key]['Energy']
        qz = data_dict2[key]['Data'][0]
        R = np.log10(data_dict2[key]['Data'][2])

        Rs = np.log10(sim_dict2[key]['Data'][2])

        l1_norm = sum(np.abs(R - Rs)) / len(R)
        tv = total_variation(R, Rs)
        idx = [i for i in range(len(qz)) if qz[i] >= 0.01]

        axes[1].plot(qz[idx], R[idx] + i * 2, 'b')
        axes[1].plot(qz[idx], Rs[idx] + i * 2, 'r')

        print('E=' + str(E) + ': ' + str(l1_norm))

    #axes[0].legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'], frameon=False, loc='lower left')

    #axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    #axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    #axes[1].legend([r'Exp ($\sigma$)', r'Calc ($\sigma$)'], frameon=False, loc='lower left')
    #axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[1].tick_params(axis='y', labelleft=False)
    #axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    #axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    #axes[0].set_ylabel(r'Normalized Reflected Intensity (arb. units)', fontsize=16)

    #shared_x = fig.add_subplot(111, frame_on=False)
    #shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    #shared_x.set_xlabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)', fontsize=16)

    #axes[0].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[1].set_xlabel('')  # Remove existing x-axis label for this subplot
    #axes[0].get_shared_x_axes().join(axes[0], shared_x)
    #axes[1].get_shared_x_axes().join(axes[1], shared_x)
    # lt.xlabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)')
    #shared_x.yaxis.set_label_coords(-0.01, 0.5)
    axes[0].set_xlim([0.01,0.1])
    axes[1].set_xlim([0.01, 0.1])
    axes[0].set_ylim([0, 1.60])
    axes[1].set_ylim([0, 1.60])

    axes[0].tick_params(axis='both', which='both', direction='in', bottom=False, left=True)
    axes[0].tick_params(axis='y', labelleft=False, left=False)
    axes[0].tick_params(axis='x', labelbottom=False, bottom=False)
    axes[1].tick_params(axis='y', labelleft=False, left=False)
    axes[1].tick_params(axis='x', labelbottom=False, bottom=False)
    plt.tight_layout()
    plt.show()
    """
    """
    example1 = np.loadtxt("//cabinet/work$/lsk601/My Documents\Master-data/test2-data.csv", skiprows=1)
    example2 = np.loadtxt("//cabinet/work$/lsk601/My Documents\Master-data/test3-data.csv", skiprows=1)

    plt.figure(2)
    plt.plot(example1[:, 0], example1[:, 1], 'r')
    plt.plot(example1[:, 4], example1[:, 5] + 0.1, 'b')
    plt.plot(example2[:, 0], example2[:, 1] - 0.8, 'r')
    plt.plot(example2[:, 4], example2[:, 5] - 0.5, 'b')
    plt.xticks([i for i in range(2, 22, 2)])
    plt.yticks([])
    # plt.ylim([0.1,0.9])
    # plt.gca().xaxis.grid(True)
    plt.legend(['Total', '600.18 eV'], loc='upper right')
    plt.xlabel('Differential Evolution Iteration')
    plt.ylabel('Cost Function (arb. units)')

    plt.show()
    """
    """
    # Boundary and Weights ----------------------------------------------------------------------------------------
    fname1 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_test2.h5"
    fname2 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_test3.h5"
    fname3 = "//cabinet/work$/lsk601/My Documents/SrTiO3-LaMnO3/Pim4uc_unitCell_test4.h5"

    data1, data_dict1, sim_dict1 = ds.ReadDataHDF5(fname1)
    data2, data_dict2, sim_dict2 = ds.ReadDataHDF5(fname2)
    data3, data_dict3, sim_dict3 = ds.ReadDataHDF5(fname3)


    keys = ["29_799.82_S", "24_700.14_S", "19_600.18_S", "14_499.93_S"]

    fig, axes = plt.subplots(2,1)

    for i, key in enumerate(keys):

        E = data_dict1[key]['Energy']
        qz = data_dict1[key]['Data'][0]
        R = np.log10(data_dict1[key]['Data'][2])

        Rs = np.log10(sim_dict1[key]['Data'][2])
        
        l1_norm = sum(np.abs(R - Rs))/len(R)

        idx = [i for i in range(len(qz)) if qz[i] >= 0.01 ]

        axes[0].plot(qz[idx], R[idx]+i*2.5+8, 'b')
        axes[0].plot(qz[idx], Rs[idx]+i*2.5+8, 'r')


        print('E=' + str(E) + ': ' + str(l1_norm))

    for i, key in enumerate(keys):
        E = data_dict2[key]['Energy']
        qz = data_dict2[key]['Data'][0]
        R = np.log10(data_dict2[key]['Data'][2])

        Rs = np.log10(sim_dict2[key]['Data'][2])

        l1_norm = sum(np.abs(R - Rs)) / len(R)

        idx = [i for i in range(len(qz)) if qz[i] >= 0.01]

        axes[1].plot(qz[idx], R[idx] + i * 2.5 + 8, 'b')
        axes[1].plot(qz[idx], Rs[idx] + i * 2.5 + 8, 'r')

        print('E=' + str(E) + ': ' + str(l1_norm))

    axes[0].legend([r'Exp ($\sigma$)',r'Calc ($\sigma$)'],frameon=False)
    axes[0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[0].tick_params(axis='y', labelleft=False)
    axes[0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].legend([r'Exp ($\sigma$)',r'Calc ($\sigma$)'],frameon=False)
    axes[1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    axes[1].tick_params(axis='y', labelleft=False)
    axes[1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1].set_xlabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)', fontsize=16)

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel(r'Normalized Reflected Intensity (arb.units)', fontsize=16)

    axes[0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0].get_shared_y_axes().join(axes[0], shared_y)
    axes[1].get_shared_y_axes().join(axes[1], shared_y)
    #lt.xlabel(r'Momentum Transfer, $q_{z}$ ($\AA^{-1}$)')
    shared_y.yaxis.set_label_coords(-0.01, 0.5)
    plt.tight_layout()

    #plt.yticks([])


    example1 = np.loadtxt("//cabinet/work$/lsk601/My Documents\Master-data/test2-data-v2.csv", skiprows=1)
    example2 = np.loadtxt("//cabinet/work$/lsk601/My Documents\Master-data/test3-data-v2.csv", skiprows=1)

    plt.figure(2)
    plt.plot(example1[:, 0], example1[:, 1], 'r')
    plt.plot(example2[:, 0], example2[:, 1], 'b')  # -0.8
    plt.plot(example1[:, 4], example1[:, 5], 'r')  # +0.1
    plt.plot(example2[:, 4], example2[:, 5], 'b')  # -0.5
    plt.xticks([i for i in range(2,22, 2)])
    #plt.yticks([0,0.2,0.4,0.6,0.8])
    #plt.ylim([0.1,0.9])
    #plt.gca().xaxis.grid(True)
    plt.tick_params(axis='both', which='both', direction='in', top=True, right=True)
    plt.tick_params(axis='y', labelleft=False)
    plt.minorticks_on()
    plt.legend(['Case (a)', 'Case (b)'])
    plt.xlabel('Differential Evolution Iteration', fontsize=16)
    plt.ylabel('Cost Function (arb. units)', fontsize=16)




    plt.show()
    """
    """
    sample = ds.ReadSampleHDF5(fname)
    data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
    sample.energy_shift()

    name1 = '38_455.73_S'
    name2 = '42_459.75_S'

    data1 = data_dict[name1]['Data']
    data2 = data_dict[name2]['Data']
    print(data_dict[name1]['Energy'],data_dict[name2]['Energy'])
    plt.figure(1)
    plt.plot(data1[0], data1[2])
    plt.xlabel(r'Momentum Transfer, $q_{z}$ ($\mathrm{\AA^{-1}}$)')
    plt.ylabel('Reflectivity')
    plt.yscale('log')

    plt.figure(2)
    plt.plot(data2[0], data2[2])
    plt.xlabel(r'Momentum Transfer, $q_{z}$ ($\mathrm{\AA^{-1}}$)')
    plt.ylabel('Reflectivity')
    plt.yscale('log')

    plt.show()

    """
    """

    # ----------------------------------- Cost Function Testing ------------------------------------------------------#

    t = np.linspace(0.5, 3, num=31)
    y = np.sin(t)
    #y = y + np.random.uniform(-0.25, 0.25, size=(1, 31))[0]
    y[15] = 5*np.sin(t[10])
    y[5] = 3*np.sin(t[5])
    #y[22] = -1*np.sin(t[22])
    #y[27] = 6*np.sin(t[27])

    #for i in range(3,28):
    #    y[i] = y[i]+20

    # L1-norm
    params = [t, y, 'L1']

    guess = [0.75, 1.25, 0]

    res1 = optimize.least_squares(cost_function,guess,args=params)
    y1 = np.sin(t*res1.x[1])*res1.x[0] + res1.x[2]

    # L2- norm

    params = [t, y, 'L2']


    res2 = optimize.least_squares(cost_function, guess, args=params)
    y2 = np.sin(t * res2.x[1]) * res2.x[0] + res2.x[2]

    # chi

    params = [t, y, 'chi']


    res3 = optimize.least_squares(cost_function, guess, args=params)
    y3 = np.sin(t * res3.x[1]) * res3.x[0] + res3.x[2]

    # arctan

    params = [t, y, 'arctan']


    res4 = optimize.least_squares(cost_function, guess, args=params)
    y4 = np.sin(t * res4.x[1]) * res4.x[0] + res4.x[2]

    fig, axes = plt.subplots(2,2)
    axes[0,0].plot(t, y, 'x')
    axes[0,0].plot(t, y1)
    axes[0,0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[0,0].tick_params(axis='y', labelleft=False)
    #axes[0, 0].tick_params(axis='x', labelleft=False)
    axes[0,0].legend(['Data', 'L1-norm'], frameon=False, loc='upper right')
    axes[0,0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0,0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0,0].set_ylim([axes[0,0].get_ylim()[0], 5])
    #axes[0,0].set_xlabel('t (arb. units)')
    #axes[0,0].set_ylabel(r'$Asin(\omega t)$')


    axes[0,1].plot(t, y, 'x')
    axes[0,1].plot(t, y2)
    axes[0, 1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[0, 1].tick_params(axis='y', labelleft=False)
    #axes[0, 1].tick_params(axis='x', labelleft=False)
    axes[0, 1].legend(['Data', 'L2-norm'], frameon=False, loc='upper right')
    axes[0, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[0, 1].set_ylim([axes[0, 0].get_ylim()[0], 5])


    axes[1,0].plot(t, y, 'x')
    axes[1,0].plot(t, y3)
    axes[1,0].legend(['Data', 'Chi-Square'])
    axes[1, 0].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[1, 0].tick_params(axis='y', labelleft=False)
    axes[1, 0].legend(['Data', 'Chi-Square'], frameon=False, loc='upper right')
    axes[1, 0].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 0].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 0].set_ylim([axes[0, 0].get_ylim()[0], 5])

    axes[1,1].plot(t, y, 'x')
    axes[1,1].plot(t, y4)
    axes[1,1].legend(['Data','Arctan'])
    axes[1, 1].tick_params(axis='both', which='both', direction='in', top=True, right=True)
    #axes[1, 1].tick_params(axis='y', labelleft=False)
    axes[1, 1].legend(['Data', 'Arctan'], frameon=False, loc='upper right')
    axes[1, 1].xaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].yaxis.set_minor_locator(ticker.AutoMinorLocator())
    axes[1, 1].set_ylim([axes[0, 0].get_ylim()[0], 5])
    #axes[1,1].set_xlabel('t (arb. units)')
    #axes[1,1].set_ylabel(r'$Asin(\omega t)$')

    shared_x = fig.add_subplot(111, frame_on=False)
    shared_x.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_x.set_xlabel('t (arb. units)', fontsize=12)

    axes[1, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1, 1].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1, 0].get_shared_x_axes().join(axes[1, 0], shared_x)
    axes[1, 1].get_shared_x_axes().join(axes[1, 1], shared_x)

    shared_y = fig.add_subplot(111, frame_on=False)
    shared_y.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    shared_y.set_ylabel(r'$Asin(\omega t)$', fontsize=12)

    axes[0, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[1, 0].set_ylabel('')  # Remove existing x-axis label for this subplot
    axes[0, 0].get_shared_y_axes().join(axes[0, 0], shared_x)
    axes[1, 0].get_shared_y_axes().join(axes[1, 0], shared_x)

    plt.tight_layout()

    plt.show()

    plt.figure(5)
    plt.plot(t,y,'kx')
    plt.plot(t,y1)
    plt.plot(t, y2)
    plt.plot(t, y3,'--')
    plt.plot(t,y4, '--')
    plt.legend(['Data', 'L1-norm','L2-norm','Chi-Square', 'Arctan'],loc='upper right')
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', top=True, right=True)
    plt.xlabel('t (arb. units)', fontsize=12)
    plt.ylabel(r'$Asin(\omega t)+b$', fontsize=12)
    plt.xlim(0.25,4)
    plt.show()

    # Global Minimum Example
    """

