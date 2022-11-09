from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt
import numpy as np
import time
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.figure as Figure
import sys
import material_structure as ms
import os
import pyqtgraph as pg
import data_structure as ds
import copy
import global_optimization as go
from PyQt5.QtCore import QObject, QThread, pyqtSignal
from scipy import signal
import h5py
import multiprocessing as mp

global stop
stop = False


def stringcheck(string):

    # checks to make sure that the roughness and whatnot is in the correct format
    num = 0
    correctFormat = True
    if len(string) == 0:
        correctFormat = False
    else:
        first = True
        for char in string:
            if first:
                if char == '.':
                    correctFormat = False
                elif not char.isdigit():
                    correctFormat = False
                first = False

            elif not char.isdigit():
                if char == '.':
                    num = num + 1
                    if num > 1:
                        correctFormat = False
                else:
                    correctFormat = False


    return correctFormat


class compoundInput(QDialog):
    def __init__(self):
        super().__init__()
        self.val = []
        pagelayout = QVBoxLayout()
        infolayout = QGridLayout()
        formula = QLabel('Formula: ')
        self.formula = QLineEdit()
        self.formula.editingFinished.connect(self.formulaDone)

        thickness = QLabel('Thickness (A): ')
        self.thickness = QLineEdit()
        self.thickness.setText('10')

        density = QLabel('Density (g/cm^3): ')
        self.density = QLineEdit()

        roughness = QLabel('Roughness (A): ')
        self.roughness = QLineEdit()
        self.roughness.setText('2')

        linkedroughnesslayout = QHBoxLayout()
        #consider including a tab to select linked roughness
        linkedroughness = QLabel('Linked Roughness (A): ')
        self.linkedroughness = QLineEdit()
        self.linkedroughness.setHidden(True)

        self.checkbox = QCheckBox()
        self.checkbox.stateChanged.connect(self.linkedroughnessbox)
        self.checkboxstate = 0

        linkedroughnesslayout.addWidget(self.checkbox)
        linkedroughnesslayout.addWidget(linkedroughness)



        infolayout.addWidget(formula, 0,0)
        infolayout.addWidget(thickness,1,0)
        infolayout.addWidget(density,2,0)
        infolayout.addWidget(roughness,3,0)
        infolayout.addLayout(linkedroughnesslayout,4,0)

        infolayout.addWidget(self.formula,0,1)
        infolayout.addWidget(self.thickness,1,1)
        infolayout.addWidget(self.density,2,1)
        infolayout.addWidget(self.roughness,3,1)
        infolayout.addWidget(self.linkedroughness,4,1)

        enterButton = QPushButton('Enter')
        enterButton.clicked.connect(self.inputComplete)

        self.errorMessage = QLabel('')

        pagelayout.addLayout(infolayout)
        pagelayout.addWidget(enterButton)
        pagelayout.addWidget(self.errorMessage)
        self.setLayout(pagelayout)

    def formulaDone(self):
        cwd = os.getcwd()
        filename =  'Perovskite_Density.txt'

        found = False
        with open(filename) as file:
            for line in file:
                myformula = line.split()[0]
                mydensity = line.split()[1]
                if not found:
                    if self.formula.text() == myformula:
                        self.density.setText(mydensity)
                        found = True
                    else:
                        self.density.clear()

    def linkedroughnessbox(self):
        self.checkboxstate = self.checkbox.checkState()
        if self.checkboxstate > 0:
            self.linkedroughness.setHidden(False)
        else:
            self.linkedroughness.setHidden(True)

    def inputComplete(self):

        finished = True
        # gets the elements and their stoichiometry
        myElements = ms.find_stoichiometry(self.formula.text())  # gets the elements and their stoichiometry

        # gets the density
        myThickness = self.thickness.text()
        thicknessCorrect = stringcheck(myThickness)

        myDensity = self.density.text()
        densityCorrect = stringcheck(myDensity)

        myRoughness = self.roughness.text()
        roughnessCorrect = stringcheck(myRoughness)

        myLinkedroughness = self.linkedroughness.text()

        linkedroughnessCorrect = True
        if myLinkedroughness != '':
            linkedroughnessCorrect = stringcheck(myLinkedroughness)

        if not(thicknessCorrect) or not(densityCorrect) or not(roughnessCorrect) or not(linkedroughnessCorrect):
            if not(thicknessCorrect):
                self.errorMessage.setText('Please check thickness!')
            elif not(densityCorrect):
                self.errorMessage.setText('Please check density!')
            elif not (roughnessCorrect):
                self.errorMessage.setText('Please check roughness!')
            elif not (linkedroughnessCorrect):
                self.errorMessage.setText('Please check linked roughness!')
        else:
            molar_mass = 0
            elements = list(myElements[0].keys())
            # gets the molar mass of the compound
            for ele in elements:
                stoich = myElements[0][ele].stoichiometry

                if ele[-1].isdigit():
                    ele = ele.rstrip(ele[-1])

                molar_mass = molar_mass + ms.atomic_mass(ele)*stoich

            tempArray = []
            if myLinkedroughness == '':
                for ele in elements:
                    stoich = myElements[0][ele].stoichiometry
                    density = float(myDensity)
                    if ele[-1].isdigit():
                        ff_ele = ele.rstrip(ele[-1])
                    else:
                        ff_ele = ele

                    tempArray.append([ele, myThickness, str(density*float(stoich)/molar_mass), myRoughness, False, ff_ele, stoich])
            else:
                for ele in elements:
                    stoich = myElements[0][ele].stoichiometry
                    density = float(myDensity)
                    if ele[-1].isdigit():
                        ff_ele = ele.rstrip(ele[-1])
                    else:
                        ff_ele = ele

                    tempArray.append([ele, myThickness, str(density*float(stoich)/molar_mass), myRoughness, myLinkedroughness, ff_ele, stoich])

            self.val = np.array(tempArray)
            self.accept()






class variationWidget(QDialog):
    def __init__(self, mainWidget, sample):
        super().__init__()

        pagelayout = QHBoxLayout()

        self.elelayout = QVBoxLayout()
        self.mainWidget = mainWidget
        self.mainWidget.layerBox.currentIndexChanged.connect(self.changeElements)


        addButton = QPushButton('Add')
        addButton.clicked.connect(self.addVarEle)
        # Create the connect function
        deleteButton = QPushButton('Delete')
        deleteButton.clicked.connect(self.deleteVarEle)
        # create the connect function

        # buttons for adding and removing

        self.sample = sample

        self.radiobutton = QRadioButton()
        idx = self.mainWidget.layerBox.currentIndex()


        # Setting up the element variation check boxes
        for j in range(len(list(self.sample.structure[idx].keys()))):
            ele = list(self.sample.structure[idx].keys())[j]
            self.mainWidget.elementBox.addItem(ele)


        self.mainWidget.elementBox.currentIndexChanged.connect(self.mainWidget.setTableVar)
        self.elelayout.addWidget(self.mainWidget.elementBox)

        self.mainWidget.setTableVar()

        self.elelayout.addWidget(addButton)
        self.elelayout.addWidget(deleteButton)

        self.mainWidget.varTable.setRowCount(2)
        self.mainWidget.varTable.setColumnCount(3)
        self.mainWidget.varTable.setHorizontalHeaderLabels(
            ['Name', 'Ratio', 'Scattering Factor'])

        pagelayout.addLayout(self.elelayout)
        pagelayout.addWidget(self.mainWidget.varTable)

        self.setLayout(pagelayout)


    def changeElements(self):

        # this function simply changes the elements based on current layer
        # needs to be changed when adding and removing layers

        self.mainWidget.change_elements = True
        idx = self.mainWidget.layerBox.currentIndex()
        self.mainWidget.elementBox.clear()

        for j in range(len(self.mainWidget.structTableInfo[idx])):
            ele = self.mainWidget.structTableInfo[idx][j][0]
            self.mainWidget.elementBox.addItem(ele)

        self.mainWidget.change_elements = False
        self.mainWidget.elementBox.setCurrentIndex(self.mainWidget.element_index)

    def addVarEle(self):
        current_layer = self.mainWidget.layerBox.currentIndex()
        current_element = self.mainWidget.elementBox.currentIndex()

        element = self.mainWidget.structTableInfo[current_layer][current_element][0]

        row = len(self.mainWidget.varData[element][current_layer][0])
        for lay in range(len(self.mainWidget.varData[element])):
            if type(self.mainWidget.varData[element][lay][0]) == np.ndarray:
                self.mainWidget.varData[element][lay][0] = np.append(self.mainWidget.varData[element][lay][0], '')
                self.mainWidget.varData[element][lay][1] = np.append(self.mainWidget.varData[element][lay][1], '')
                self.mainWidget.varData[element][lay][2] = np.append(self.mainWidget.varData[element][lay][2], '')

                self.mainWidget.magData[element][lay][0] = np.append(self.mainWidget.magData[element][lay][0], '')
                self.mainWidget.magData[element][lay][1] = np.append(self.mainWidget.magData[element][lay][1], '')
                self.mainWidget.magData[element][lay][2] = np.append(self.mainWidget.magData[element][lay][2], '')
            else:
                self.mainWidget.varData[element][lay][0].append('')  # add another element to name list
                self.mainWidget.varData[element][lay][1].append('')  # add another element to name list
                self.mainWidget.varData[element][lay][2].append('')  # add another element to name list

                self.mainWidget.magData[element][lay][0].append('')  # make appropriate changes to magnetic data
                self.mainWidget.magData[element][lay][1].append('')
                self.mainWidget.magData[element][lay][2].append('')


        #row = self.mainWidget.varTable.rowCount()
        self.mainWidget.varTable.setRowCount(row + 1)

    def deleteVarEle(self):

        current_layer = self.mainWidget.layerBox.currentIndex()
        current_element = self.mainWidget.elementBox.currentIndex()

        element = self.mainWidget.structTableInfo[current_layer][current_element][0]

        row = len(self.mainWidget.varData[element][current_layer][0])

        if row != 2:
            for lay in range(len(self.mainWidget.varData[element])):

                if type(self.mainWidget.varData[element][lay][0]) == np.ndarray:
                    self.mainWidget.varData[element][lay][0] = self.mainWidget.varData[element][lay][0][:-1]
                    self.mainWidget.varData[element][lay][1] = self.mainWidget.varData[element][lay][1][:-1]
                    self.mainWidget.varData[element][lay][2] = self.mainWidget.varData[element][lay][2][:-1]

                    self.mainWidget.magData[element][lay][0] = self.mainWidget.magData[element][lay][0][:-1]
                    self.mainWidget.magData[element][lay][1] = self.mainWidget.magData[element][lay][1][:-1]
                    self.mainWidget.magData[element][lay][2] = self.mainWidget.magData[element][lay][2][:-1]
                else:
                    self.mainWidget.varData[element][lay][0].pop()  # add another element to name list
                    self.mainWidget.varData[element][lay][1].pop()  # add another element to name list
                    self.mainWidget.varData[element][lay][2].pop()  # add another element to name list

                    self.mainWidget.magData[element][lay][0].pop() # make changes to magnetic data
                    self.mainWidget.magData[element][lay][1].pop()
                    self.mainWidget.magData[element][lay][2].pop()

            self.mainWidget.varTable.setRowCount(row-1)

class ReadOnlyDelegate(QStyledItemDelegate):
    def createEditor(self, parent, option, index):
        return

class magneticWidget(QDialog):
    def __init__(self, mainWidget, sample):
        super().__init__()

        pagelayout = QHBoxLayout()

        self.mainWidget = mainWidget


        # buttons for adding and removing

        self.sample = sample

        idx = self.mainWidget.layerBox.currentIndex()
        self.mainWidget.layerBox.currentIndexChanged.connect(self.mainWidget.setTableMag)
        # Magnetization direction Widget format
        magLabel = QLabel('Magnetization Direction')
        magLayout = QVBoxLayout()

        self.mainWidget.magDirBox.addItem('x-direction')
        self.mainWidget.magDirBox.addItem('y-direction')
        self.mainWidget.magDirBox.addItem('z-direction')

        magLayout.addWidget(magLabel)
        magLayout.addWidget(self.mainWidget.magDirBox)
        magLayout.addStretch(1)

        self.mainWidget.magDirBox.currentIndexChanged.connect(self.magDirectionChange)

        self.mainWidget.magTable.setRowCount(3)
        self.mainWidget.magTable.setColumnCount(2)
        self.mainWidget.magTable.setHorizontalHeaderLabels(
            ['Magnetic Density (mol/cm^3)', 'Scattering Factor'])

        self.mainWidget.setTableMag()
        #self.mainWidget.magButton.clicked.connect(self.setTable)
        pagelayout.addWidget(self.mainWidget.magTable)
        pagelayout.addLayout(magLayout)

        self.setLayout(pagelayout)


    def magDirectionChange(self):
        lay = self.mainWidget.layerBox.currentIndex()
        mag = self.mainWidget.magDirBox.currentIndex()

        if mag == 0:
            self.mainWidget.magDirection[lay] = 'x'
        elif mag == 1:
            self.mainWidget.magDirection[lay] = 'y'
        elif mag == 2:
            self.mainWidget.magDirection[lay] = 'z'


class sampleWidget(QWidget):
    def __init__(self, sample):
        super(sampleWidget,self).__init__()
        self.sample = sample  # variable used to define sample info
        self.structTableInfo = []  # used to keep track of the table info instead of constantly switching
        self.parameterFit = []
        self.varData = {ele: [[['',''],['',''],['','']] for i in range(len(sample.structure))] for ele in sample.myelements}
        self.eShift = dict()  # keep track of the energy shift
        self.currentVal = []
        self.change_eShift = True
        self.varTable = QTableWidget()
        self.elementBox = QComboBox()
        self.variationElements = sample.poly_elements

        self.magData = {ele: [[[''], [''], ['']] for i in range(len(sample.structure))] for ele in
                        sample.myelements}
        self.magGo = True
        self.magDirection = ['z' for i in range(len(sample.structure))]
        self.magDirBox = QComboBox()
        self.getData() # gets the element variation and magnetic information

        self.magTable = QTableWidget()


        self.change_elements = False
        self.element_index = 0

        self.previousLayer = 0
        self.changeLayer = False
        self.firstStruct = True

        self._step_size = '0.1'
        self.step_size = QLineEdit()
        self.step_size.textChanged.connect(self.changeStepSize)
        self.step_size.setText(self._step_size)
        self.step_size.setMaximumWidth(100)
        step_size_label = QLabel('Step Size (Å):')
        step_size_label.setMaximumWidth(65)
        step_size_layout = QHBoxLayout()
        step_size_layout.addWidget(step_size_label)
        step_size_layout.addWidget(self.step_size)

        # Initializes the sample definition widget
        pagelayout = QHBoxLayout()  # page layout

        cblayout = QVBoxLayout()  # combobox and button layout


        # bottons for adding, copying, and deleteing layers
        addlayerButton = QPushButton('Add Layer')  # add layer
        addlayerButton.clicked.connect(self._addLayer)

        copylayerButton = QPushButton('Copy Current Layer')  # copy layer
        copylayerButton.clicked.connect(self._copyLayer)

        deletelayerButton = QPushButton('Remove Current Layer')  # delete layer
        deletelayerButton.clicked.connect(self._removeLayer)

        # Layer Box
        self.structInfo = self.sample.structure
        self.layerBox = QComboBox(self)


        layerList = []
        for i in range(len(self.sample.structure)):
            if i == 0:
                layerList.append('Substrate')
            else:
                layerList.append('Layer '+str(i))

        # change this for an arbitrary sample model
        self.layerBox.addItems(layerList)
        self.layerBox.currentIndexChanged.connect(self.setTable)
        # changes the table on the screen when new layer selected


        # buttons for adding and removing layers
        cblayout.addWidget(addlayerButton)
        cblayout.addWidget(copylayerButton)
        cblayout.addWidget(deletelayerButton)
        cblayout.addLayout(step_size_layout)

        # layer combo box
        cblayout.addWidget(self.layerBox)

        self.sampleInfoLayout = QStackedLayout() # stacked layout for the different parameter types

        self.structTable = QTableWidget()

        self.structTable.setRowCount(3)
        self.structTable.setColumnCount(6)
        self.structTable.setHorizontalHeaderLabels(
            ['Element', 'Thickness', 'Density', 'Roughness', 'Linked Roughness', 'Scattering Factor'])

        self._setStructFromSample(sample)


        # setTable
        self.setTable()

        self.energyShiftTable = QTableWidget()
        self.energyShiftTable.setColumnCount(1)
        self.energyShiftTable.setHorizontalHeaderLabels(['Energy Shift (eV)'])


        # Element Variation Stuff
        self.elementVariation = variationWidget(self,sample)
        self.elementMagnetic = magneticWidget(self,sample)

        self.sampleInfoLayout.addWidget(self.structTable)
        self.sampleInfoLayout.addWidget(self.elementVariation)
        self.sampleInfoLayout.addWidget(self.elementMagnetic)
        self.sampleInfoLayout.addWidget(self.energyShiftTable)

        self.structTable.viewport().installEventFilter(self)
        self.varTable.viewport().installEventFilter(self)
        self.magTable.viewport().installEventFilter(self)
        self.energyShiftTable.viewport().installEventFilter(self)
        # need to add fitting option to energy shift

        selectlayout = QVBoxLayout()

        # buttons for choosing which parameters to choose
        self.structButton = QPushButton('Structure')
        self.structButton.setStyleSheet('background: blue; color: white')
        self.structButton.clicked.connect(self._structural)
        self.structButton.clicked.connect(self.setTableVar)
        self.structButton.clicked.connect(self.setTableMag)
        selectlayout.addWidget(self.structButton)

        self.polyButton = QPushButton('Element Variation')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.polyButton.clicked.connect(self._elementVariation)
        self.polyButton.clicked.connect(self.setTableVar)
        selectlayout.addWidget(self.polyButton)

        self.magButton = QPushButton('Magnetic')
        self.magButton.setStyleSheet('background: lightGrey')
        self.magButton.clicked.connect(self._magnetic)
        self.magButton.clicked.connect(self.setTableMag)
        selectlayout.addWidget(self.magButton)

        self.shiftButton = QPushButton('Energy Shift')  # energy shift button
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.shiftButton.clicked.connect(self._energy_shift)
        selectlayout.addWidget(self.shiftButton)

        self.structBool = True
        self.polyBool = False
        self.magBool = False

        dpButton = QPushButton('Density Profile')
        dpButton.clicked.connect(self._densityprofile)
        dpButton.setStyleSheet("background-color : cyan")
        selectlayout.addWidget(dpButton)

        pagelayout.addLayout(cblayout)
        pagelayout.addLayout(self.sampleInfoLayout)
        pagelayout.addLayout(selectlayout)

        mylayout = QVBoxLayout()
        mylayout.addLayout(pagelayout)

        # Adding the plotting Widget
        self.densityWidget = pg.PlotWidget()
        self.densityWidget.setBackground('w')

        self.densityWidget.addLegend()

        mylayout.addWidget(self.densityWidget)

        # change parameters when clicked
        self.structTable.itemChanged.connect(self.changeStructValues)
        self.varTable.itemChanged.connect(self.changeVarValues)
        self.magTable.itemChanged.connect(self.changeMagValues)
        self.energyShiftTable.itemChanged.connect(self.changeEShiftValues)
        delegate = ReadOnlyDelegate()
        self.structTable.setItemDelegateForColumn(0,delegate)
        self.setLayout(mylayout)

    def _energy_shift(self):
        self.sampleInfoLayout.setCurrentIndex(3)
        self.structButton.setStyleSheet('background: lightGrey')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.magButton.setStyleSheet('background: lightGrey')
        self.shiftButton.setStyleSheet('background: blue; color: white')
        self.setTableEShift()

    def setTableEShift(self):
        self.energyShiftTable.blockSignals(True)
        keys = list(self.eShift.keys())
        self.energyShiftTable.setColumnCount(len(keys))
        self.energyShiftTable.setRowCount(1)
        self.energyShiftTable.setHorizontalHeaderLabels(keys)
        for column,key in enumerate(keys):
            item = QTableWidgetItem(str(self.eShift[key]))
            self.energyShiftTable.setItem(0, column, item)


        copy_of_list = copy.deepcopy(self.parameterFit)
        my_fits = []
        for fit in copy_of_list:
            for column, key in enumerate(keys):
                if len(fit) == 3:
                    if key[:3] == 'ff-':
                        if fit[0] == 'SCATTERING FACTOR' and fit[1] == 'STRUCTURAL' and fit[2] == key[3:]:
                            my_fits.append(column)
                        else:
                            self.energyShiftTable.item(0, column).setBackground(QtGui.QColor(255, 255, 255))
                    elif key[:3] == 'ffm':
                        if fit[0] == 'SCATTERING FACTOR' and fit[1] == 'MAGNETIC' and fit[2] == key[4:]:
                            my_fits.append(column)
                        else:
                            self.energyShiftTable.item(0, column).setBackground(QtGui.QColor(255, 255, 255))

        for col in my_fits:
            self.energyShiftTable.item(0, col).setBackground(QtGui.QColor(150, 255, 150))

        self.energyShiftTable.blockSignals(False)

    def changeEShiftValues(self):
        column = self.energyShiftTable.currentColumn()
        key = self.energyShiftTable.horizontalHeaderItem(column).text()
        value = self.energyShiftTable.item(0,column).text()
        self.eShift[key] = float(value)

        # make changes to the sample
        if key[:3] == 'ff-':
            sf = key[3:]
            self.sample.eShift[sf] = float(value)
        elif key[:3] == 'ffm':
            sf = key[4:]
            self.sample.mag_eShift[sf] = float(value)

        for idx, fit in enumerate(self.parameterFit):
            if key[:3] == 'ff-':
                if fit == ['SCATTERING FACTOR', 'STRUCTURAL', key[3:]]:
                    upper = str(float(value) + 0.5)
                    lower = str(float(value) - 0.5)
                    self.currentVal[idx] = [value, [lower, upper]]
            elif key[:4] == 'ffm':
                if fit == ['SCATTERING FACTOR', 'MAGNETIC', key[4:]]:
                    upper = str(float(value) + 0.5)
                    lower = str(float(value) - 0.5)
                    self.currentVal[idx] = [value, [lower, upper]]

    def changeMagValues(self):

        layer = self.layerBox.currentIndex()

        column = self.magTable.currentColumn()
        row = self.magTable.currentRow()

        if self.magTable.item(row, column) is not None and self.magTable.verticalHeaderItem(row) is not None:
            name = self.magTable.verticalHeaderItem(row).text()
            element = ''
            idx = 0
            for i in range(len(self.structTableInfo[layer])):

                ele = self.structTableInfo[layer][i][0]
                if name in self.magData[ele][layer][0]:
                    idx = list(self.magData[ele][layer][0]).index(name)
                    element = copy.copy(ele)

            value = self.magTable.item(row, column).text()


            prev_value = self.magData[element][layer][column + 1][idx]
            if column == 1:
                prev_dict_name = 'ffm-'+prev_value
                if prev_value != '':

                    del self.eShift[prev_dict_name]
                dict_name = 'ffm-' + value
                self.eShift[dict_name] = 0

                if prev_value != value:
                    for i in range(len(self.magData[element])):
                        inLayer = False
                        for j in range(len(self.structTableInfo[i])):
                            if element == self.structTableInfo[i][j][0]:
                                inLayer = True
                        if inLayer:
                            self.magData[element][i][column+1][idx] = value
                            if self.magData[element][i][1][idx] == '':
                                self.magData[element][i][1][idx] = '0'


            self.magData[element][layer][column+1][idx] = value

            copy_of_list = copy.deepcopy(self.parameterFit)
            for fit in copy_of_list:
                if column == 0: # density
                    if layer == fit[0] and fit[1] == 'MAGNETIC' and fit[-1] == name:
                        idx = self.parameterFit.index(fit)
                        lower = float(value) - 0.01
                        if lower < 0:
                            lower = 0

                        upper = str(float(value) + 0.01)
                        self.currentVal[idx] = [value, [str(lower), upper]]
                elif column == 1 and fit[0] == 'SCATTERING FACTOR' and fit[1] == 'MAGNETIC': # magnetic scattering factor
                    if value != prev_value and prev_value == fit[2]:
                        self.parameterFit.remove(fit)
                        self.magTable.item(row, column).setBackground(QtGui.QColor(255, 255, 255))


    def changeVarValues(self):

        # change the polymorphous parameters
        layer = self.layerBox.currentIndex()
        ele_idx = self.elementBox.currentIndex()
        element = self.structTableInfo[layer][ele_idx][0]
        column = self.varTable.currentColumn()
        row = self.varTable.currentRow()



        if self.varTable.item(row, column) is not None and not(self.change_elements):

            copy_of_list = copy.deepcopy(self.parameterFit)

            # check to see if element array has already been initialized
            empty = True
            for i in range(self.varTable.rowCount()):
                for j in range(self.varTable.colorCount()):
                    if self.varTable.item(i,j) is None:
                        item = ''
                    else:
                        item = self.varTable.item(i, j).text()
                    if item != '':
                        empty = False

            if empty:
                self.magData[element][layer] = [[element], [''], ['']]
            else:
                if len(self.magData[element][layer][0]) == 1:
                    self.magData[element][layer] = [['', ''], ['', ''], ['', '']]

            value = self.varTable.item(row, column).text()  # setting varData correctly depending on user input
            prev_value = self.varData[element][layer][column][row]

            self.varData[element][layer][column][row] = value
            # everytime we need to check is all values in the table are nothing
            # else we need to initialize a new set of polymorphs


            if column == 0:
                self.magData[element][layer][column][row] = value
                if prev_value != value and prev_value!='':
                    for i in range(len(self.varData[element])):
                        inLayer = False

                        for j in range(len(self.structTableInfo[i])):
                            if element == self.structTableInfo[i][j][0]:
                                inLayer = True

                        if inLayer:
                            self.varData[element][i][0][row] = value
                            self.magData[element][i][0][row] = value


            # changing the scattering factor
            if column == 2:
                if prev_value != value:

                    # takes into account the scattering factor change
                    prev_dict_name = 'ff-' + prev_value
                    if prev_value != '':
                        del self.eShift[prev_dict_name]
                    dict_name = 'ff-' + value
                    self.eShift[dict_name] = 0


            for fit in copy_of_list:
                if fit[0] == layer and fit[1] == 'POLYMORPHOUS' and column != 2:
                    idx = self.parameterFit.index(fit)
                    if column == 0:  # case where we changed the element variation name
                        if prev_value != value and fit[3] == prev_value:
                            self.parameterFit[idx][3] = value
                    elif column == 1:  # case for change in ratio
                        if fit[2] == element and fit[3] == self.varData[element][layer][0][row]:
                            lower = float(value)-0.2
                            if lower < 0:
                                lower = 0
                            upper = float(value) + 0.2
                            if upper > 1:
                                upper = 1

                            self.currentVal[idx] = [value, [str(lower), str(upper)]]

                elif fit[0] == 'SCATTERING FACTOR' and fit[1] == 'STRUCTURAL' and column == 2:
                    if fit[2] == prev_value and value != prev_value:
                        self.parameterFit.remove(fit)
                        self.varTable.item(row, column).setBackground(QtGui.QColor(255, 255, 255))


    def changeStructValues(self):
        layer = self.layerBox.currentIndex()
        row = self.structTable.currentRow()
        column = self.structTable.currentColumn()
        if self.structTable.item(row, column) is not None:
            copy_of_list = self.parameterFit
            value = self.structTable.item(row, column).text()
            prev_value = copy.copy(self.structTableInfo[layer][row][column])

            # change name or scattering factor if it is changed (eShift name change)
            if column == 5:
                name = self.structTable.item(row, 0).text()
                if value != prev_value:
                    # changing scattering factor info
                    prev_dict_name = 'ff-'+prev_value
                    del self.eShift[prev_dict_name]
                    dict_name = 'ff-'+ value
                    self.eShift[dict_name] = 0

                    self.eShift[dict_name] = value
                    for i in range(len(self.structTableInfo)):
                        for j in range(len(self.structTableInfo[i])):
                            if self.structTableInfo[i][j][0] == name:
                                self.structTableInfo[i][j][5] = value  # changes the scattering factor for all

            self.structTableInfo[layer][row][column] = self.structTable.item(row, column).text()

            for fit in copy_of_list:
                element = self.structTableInfo[layer][row][0]
                if fit[0] == layer and fit[1] == 'STRUCTURAL':
                    if column == 1:  # thickness
                        if fit[2] == 'COMPOUND' and fit[3] == 'THICKNESS':
                            idx = self.parameterFit.index([layer,'STRUCTURAL','COMPOUND', 'THICKNESS'])
                            lower = float(value) - 5
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 5)
                            self.currentVal[idx] = [value,[str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == element and fit[4] == 'THICKNESS':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'THICKNESS'])
                            lower = float(value) - 5
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 5)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                    elif column == 2: # density
                        if fit[2] == 'COMPOUND' and fit[3] == 'DENSITY':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'DENSITY'])
                            lower = float(value) - 0.01
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 0.01)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == element and fit[4] == 'DENSITY':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'DENSITY'])
                            lower = float(value) - 0.01
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 0.01)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                    elif column == 3: # roughness
                        if fit[2] == 'COMPOUND' and fit[3] == 'ROUGHNESS':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'ROUGHNESS'])
                            lower = float(value) - 1
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 1)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == element and fit[3] == 'ROUGHNESS':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'ROUGHNESS'])
                            lower = float(value) - 1
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 1)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                    elif column == 4: # linked roughness
                        if fit[2] == 'COMPOUND' and fit[3] == 'LINKED ROUGHNESS':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'LINKED ROUGHNESS'])
                            lower = float(value) - 1
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 1)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == 'LINKED ROUGHNESS':
                            idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'LINKED ROUGHNESS'])
                            lower = float(value) - 1
                            if lower < 0:
                                lower = 0
                            upper = str(float(value) + 1)
                            self.currentVal[idx] = [value, [str(lower), upper]]
                elif fit[0] == 'SCATTERING FACTOR' and fit[2] == element and fit[1] == 'STRUCTURAL':
                    if value != prev_value and prev_value == fit[2]:
                        self.parameterFit.remove(fit)
                        self.structTable.item(row, column).setBackground(QtGui.QColor(255, 255, 255))


    def eventFilter(self, source, event):
        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.button() == Qt.RightButton:
                idx = self.sampleInfoLayout.currentIndex()
                if idx == 0:
                    self.structure_handler()
                elif idx == 1:
                    self.var_handler()
                elif idx == 2:
                    self.mag_handler()
                elif idx == 3:
                    self.eShift_handler()
        return False

    def eShift_handler(self):
        column = self.energyShiftTable.currentColumn()
        copy_fit_list = copy.deepcopy(self.parameterFit)
        name = self.energyShiftTable.horizontalHeaderItem(column).text()
        val = self.energyShiftTable.currentItem().text()

        top_menu = QMenu()

        menu = top_menu.addMenu("Menu")

        _fit = menu.addAction('Fit')
        _remove_fit = menu.addAction('Remove Fit')

        action = menu.exec_(QtGui.QCursor.pos())

        if action == _fit:
            fit = []
            if name[:3] == 'ff-':
                fit = ['SCATTERING FACTOR', 'STRUCTURAL', name[3:]]
            elif name[:3] == 'ffm':
                fit = ['SCATTERING FACTOR', 'MAGNETIC', name[4:]]

            if fit != []:
                self.parameterFit.append(fit)
                lower = str(float(val) - 0.5)
                upper = str(float(val) + 0.5)
                self.currentVal.append([val, [lower, upper]])

        elif action == _remove_fit:
            fit = []
            if name[:3] == 'ff-':
                fit = ['SCATTERING FACTOR', 'STRUCTURAL', name[3:]]
            elif name[:3] == 'ffm':
                fit = ['SCATTERING FACTOR', 'MAGNETIC', name[4:]]

            if fit in self.parameterFit:
                idx = self.parameterFit.index(fit)
                self.parameterFit.remove(fit)
                self.currentVal.pop(idx)

        self.setTableEShift()
        self.setTable()
        self.setTableMag()
        self.setTableVar()

    def mag_handler(self):
        idx = self.sampleInfoLayout.currentIndex()  # keeps track of which parameters are to be fit
        my_layer = self.layerBox.currentIndex()  # layer index
        copy_fit_list = copy.deepcopy(self.parameterFit)
        row = self.magTable.currentRow()
        column = self.magTable.currentColumn()

        # find the element that the variation belongs too
        name = self.magTable.verticalHeaderItem(row).text()
        element = ''
        if name in list(self.magData.keys()):
            element = name
        else:
            for key in list(self.magData.keys()): # loop through all the keys
                if name in self.magData[key][my_layer][0]:
                    element = key

        top_menu = QMenu()

        menu = top_menu.addMenu("Menu")

        _fit = menu.addAction('Fit')
        _remove_fit = menu.addAction('Remove Fit')

        action = menu.exec_(QtGui.QCursor.pos())
        if action == _fit:

            # determines if the element has already been selected
            alreadySelected = False
            for fit in copy_fit_list:
                if column == 0:
                    if len(fit) == 3:
                        if fit[0] == my_layer and fit[1] == 'MAGNETIC' and fit[2] == element:
                            alreadySelected = True
                    elif len(fit) == 4:
                        if fit[0] == my_layer and fit[1] == 'MAGNETIC' and fit[2] == element and fit[2] == name:
                            alreadySelected = True
                elif column == 1:
                    scattering_factor = self.magTable.item(row,1)
                    if fit[0] == 'SCATTERING FACTOR' and fit[1] == 'MAGNETIC' and fit[2] ==  scattering_factor:
                        alreadySelected = True

            # Check to make sure that parameter not already selected
            if not alreadySelected:
                if column == 1:
                    scattering_factor = self.magTable.item(row, column).text()
                    self.currentVal.append([0,[-0.5,0.5]])
                    self.parameterFit.append(['SCATTERING FACTOR', 'MAGNETIC', scattering_factor])

                elif column == 0:
                    value = self.magTable.item(row,column).text()
                    lower = float(value) - 0.01
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 0.01
                    self.currentVal.append([float(value),[lower, upper]])
                    if name == element:
                        self.parameterFit.append([my_layer, 'MAGNETIC', element])
                    else:
                        self.parameterFit.append([my_layer, 'MAGNETIC', element, name])

        if action == _remove_fit:
            for fit in copy_fit_list:
                if column == 0:  # magnetic density
                    if len(fit) == 3 and fit[1] == 'MAGNETIC':  # not polymorphous case
                        if fit[0] == my_layer and fit[2] == name:
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                    elif len(fit) == 4 and fit[1] == 'MAGNETIC': # polymorphous case
                        if fit[0] == my_layer and fit[3] == name:
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)

                elif column == 1:  # Magnetic Scattering Factor
                    if len(fit) == 3 and fit[0] == 'SCATTERING FACTOR':
                        scattering_factor = self.magTable.item(row, 1).text()
                        if scattering_factor == fit[2] and fit[1] == 'MAGNETIC':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)

        self.setTableMag()

    def var_handler(self):
        idx = self.sampleInfoLayout.currentIndex()  # keeps track of which parameters are to be fit
        my_layer = self.layerBox.currentIndex()  # layer index
        my_ele = self.elementBox.currentIndex()  # element index
        copy_fit_list = copy.deepcopy(self.parameterFit)  # fit list

        element = self.structTableInfo[my_layer][my_ele][0]  # retrieves element in layer selected

        top_menu = QMenu()

        menu = top_menu.addMenu("Menu")

        _fit = menu.addAction('Fit')
        _remove_fit = menu.addAction('Remove Fit')

        row = self.varTable.currentRow()
        column = self.varTable.currentColumn()

        if column == 0:
            _fit.setDisabled(True)
            _remove_fit.setDisabled(True)

        action = menu.exec_(QtGui.QCursor.pos())
        if action == _fit:
            alreadySelected = False
            for fit in copy_fit_list:
                if column == 1:
                    if fit[1] == 'POLYMORPHOUS' and fit[0] == my_layer:
                        name = self.varTable.item(row, 0).text()
                        if element == fit[2] and name == fit[3]:
                            alreadySelected = True
                elif column == 2:
                    scattering_factor = self.varTable.item(row,column).text()
                    if len(fit) == 3 and scattering_factor == fit[2]:
                        alreadySelected = True

            if not alreadySelected:
                if column == 1:
                    name = self.varTable.item(row, 0).text()
                    ratio = self.varTable.item(row, 1).text()

                    if ratio != '':
                        self.parameterFit.append([my_layer, 'POLYMORPHOUS', element, name])
                        lower = float(ratio) - 0.2
                        upper = float(ratio) + 0.2
                        if lower < 0:
                            lower = 0
                        if upper > 1:
                            upper = 1

                        self.currentVal.append([float(ratio), [lower, upper]])

                elif column == 2:
                    scattering_factor = self.varTable.item(row,2).text()
                    if scattering_factor != '':
                        self.parameterFit.append(['SCATTERING FACTOR', 'STRUCTURAL',scattering_factor])
                        self.currentVal.append([0, [-0.5,0.5]])

        elif action == _remove_fit:

            for fit in copy_fit_list:
                if column == 1:
                    if len(fit) == 4 and fit[1] == 'POLYMORPHOUS':
                        name = self.varTable.item(row,0).text()
                        if my_layer == fit[0] and element == fit[2] and name == fit[3]:
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                elif column == 2:
                    if len(fit) == 3:
                        scattering_factor = self.varTable.item(row,2).text()
                        if scattering_factor == fit[2] and fit[1] == 'STRUCTURAL':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
        self.setTableVar()

    def structure_handler(self):

        idx = self.sampleInfoLayout.currentIndex()
        my_layer = self.layerBox.currentIndex()
        copy_fit_list = copy.deepcopy(self.parameterFit)
        top_menu = QMenu()

        menu = top_menu.addMenu("Menu")

        _element_fit = menu.addAction('Element Fit')
        _compound_fit = menu.addAction('Compound Fit')
        _remove_fit = menu.addAction('Remove Fit')

        row = self.structTable.currentRow()
        column = self.structTable.currentColumn()

        # only allows user to fit certain parameters
        if column == 5:
            _compound_fit.setDisabled(True)
        elif column == 0 or (column == 1 and my_layer == 0):
            _element_fit.setDisabled(True)
            _compound_fit.setDisabled(True)
            _remove_fit.setDisabled(True)
        elif column == 4: # linked roughness case
            LRough = self.structTableInfo[my_layer][row][4]
            if str(LRough).upper() == 'FALSE' or LRough == '':
                _element_fit.setDisabled(True)
                _compound_fit.setDisabled(True)
                _remove_fit.setDisabled(True)
            compound = True
            for cool in self.structTableInfo[my_layer]:
                if str(cool[4]).upper() == 'FALSE' or cool[4] == '':
                    compound = False

            if not compound:
                _compound_fit.setDisabled(True)



        action = menu.exec_(QtGui.QCursor.pos())


        # Element Mode
        if action == _element_fit:

            value = self.structTable.currentItem().text()
            element = self.structTable.item(row, 0).text()
            alreadySelected = False

            # Check to make sure parameter is not already selected
            for fit in copy_fit_list:
                # check if layer and parameter in compound mode
                n = len(fit)
                if n == 3:  # scattering factor
                    if column == 5:  # trying to fit scattering factor
                        item = self.structTable.currentItem().text()
                        if item == fit[2]:
                            alreadySelected = True

                elif n == 5:  # compound mode
                    layer = fit[0]
                    param = fit[3]
                    param_num = 0
                    if param == 'THICKNESS':
                        param_num = 1
                    elif param == 'DENSITY':
                        param_num = 2
                    elif param == 'ROUGHNESS':
                        param_num = 3
                    elif param == 'LINKED ROUGHNESS':
                        param_num = 4

                    if layer == my_layer and column == param_num:
                        idx = self.parameterFit.index(fit)
                        self.parameterFit.remove(fit)
                        self.currentVal.pop(idx)

                elif n == 5:  # element mode
                    layer = fit[0]
                    ele = fit[3]
                    param = fit[4]
                    my_ele = self.structTableInfo[idx][row][0]

                    param_num = 1
                    if param == 'THICKNESS':
                        param_num = 1
                    elif param == 'DENSITY':
                        param_num = 2
                    elif param == 'ROUGHNESS':
                        param_num = 3
                    elif param == 'LINKED ROUGHNESS':
                        param_num = 4

                    if layer == my_layer and column == param_num and ele == my_ele:
                        alreadySelected = True

            if not alreadySelected:
                if column == 1:  # thickness
                    self.parameterFit.append([my_layer, 'STRUCTURAL', 'ELEMENT', element, 'THICKNESS'])
                    lower = float(value) - 5
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 5
                    self.currentVal.append([float(value), [lower, upper]])
                elif column == 2:  # density
                    self.parameterFit.append([my_layer, 'STRUCTURAL' , 'ELEMENT', element, 'DENSITY'])
                    lower = float(value) - 0.01
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 0.01
                    self.currentVal.append([float(value), [lower, upper]])
                elif column == 3:  # roughness
                    self.parameterFit.append([my_layer, 'STRUCTURAL' , 'ELEMENT', element, 'ROUGHNESS'])
                    lower = float(value) - 1
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 1
                    self.currentVal.append([float(value), [lower, upper]])
                elif column == 4:  # linked roughness
                    self.parameterFit.append([my_layer, 'STRUCTURAL', 'ELEMENT', element, 'LINKED ROUGHNESS'])
                    lower = float(value) - 1
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 1
                    self.currentVal.append([float(value), [lower, upper]])
                elif column == 5:  # scattering factor
                    scattering_factor = self.structTable.item(row, 5).text()
                    if scattering_factor[0] != '[':
                        self.parameterFit.append(['SCATTERING FACTOR', 'STRUCTURAL', scattering_factor])
                        self.currentVal.append([0,[-0.5,0.5]])

        elif action == _compound_fit:

            # retrieve minimum value in the row
            my_vals = list()
            for i in range(self.structTable.rowCount()):
                my_vals.append(float(self.structTable.item(i,column).text()))

            my_row = my_vals.index(min(my_vals))
            value = self.structTable.item(my_row, column).text()  # minimum value

            alreadySelected = False
            for fit in copy_fit_list:
                mode = fit[2]
                if mode == 'COMPOUND':  # compound check
                    layer = fit[0]
                    param = fit[3]

                    param_n = 0
                    if param == 'THICKNESS':
                        param_n = 1
                    elif param == 'DENSITY':
                        param_n = 2
                    elif param == 'ROUGHNESS':
                        param_n = 3
                    elif param == 'LINKED ROUGHNESS':
                        param_n = 4

                    if layer == my_layer and param_n == column:
                        alreadySelected = True

                elif mode == 'ELEMENT':  # element check
                    layer = fit[0]
                    param = fit[4]
                    param_n = 1
                    if param == 'THICKNESS':
                        param_n = 1
                    elif param == 'DENSITY':
                        param_n = 2
                    elif param == 'ROUGHNESS':
                        param_n = 3
                    elif param == 'LINKED ROUGHNESS':
                        param_n = 4


                    if param_n == column and layer == my_layer:
                        idx = self.parameterFit.index(fit)
                        self.parameterFit.remove(fit)
                        self.currentVal.pop(idx)

            if not alreadySelected:
                if column == 1:  # thickness
                    if my_layer != 0:
                        my_fit = [my_layer, 'STRUCTURAL','COMPOUND', 'THICKNESS', my_row]
                        self.parameterFit.append(my_fit)
                        lower = float(value) - 5
                        if lower < 0:
                            lower = 0
                        upper = float(value) + 5
                        self.currentVal.append([float(value), [lower, upper]])
                elif column == 2:  # density
                    my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'DENSITY', my_row]
                    self.parameterFit.append(my_fit)
                    lower = float(value) - 0.01
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 0.01
                    self.currentVal.append([float(value), [lower, upper]])
                elif column == 3:  # roughness
                    my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'ROUGHNESS', my_row]
                    self.parameterFit.append(my_fit)
                    lower = float(value) - 1
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 1
                    self.currentVal.append([float(value), [lower, upper]])
                elif column == 4:  # linked roughness
                    my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'LINKED ROUGHNESS', my_row]
                    self.parameterFit.append(my_fit)
                    lower = float(value) - 1
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 1
                    self.currentVal.append([float(value), [lower, upper]])

        elif action == _remove_fit:
            element = self.structTableInfo[my_layer][row][0]
            scattering_factor = self.structTableInfo[my_layer][row][5]
            for fit in copy_fit_list:
                n = len(fit)
                if column == 1:
                    mode = fit[2]
                    if mode == "ELEMENT":
                        ele = fit[3]
                        if ele == element and fit[4] == 'THICKNESS':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)


                    elif mode == 'COMPOUND':
                        if fit[3] == 'THICKNESS':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                elif column == 2:
                    mode = fit[2]
                    if mode == "ELEMENT":
                        ele = fit[3]
                        if ele == element and fit[4] == 'DENSITY':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                    elif mode == 'COMPOUND':
                        if fit[3] == 'DENSITY':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                elif column == 3:
                    mode = fit[2]
                    if mode == "ELEMENT":
                        ele = fit[3]
                        if ele == element and fit[4] == 'ROUGHNESS':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                    elif mode == 'COMPOUND':
                        if fit[3] == 'ROUGHNESS':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                elif column == 4:
                    mode = fit[2]
                    if mode == "ELEMENT":
                        ele = fit[3]
                        if ele == element and fit[4] == 'LINKED ROUGHNESS':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                    elif mode == 'COMPOUND':
                        if fit[3] == 'LINKED ROUGHNESS':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)
                elif column == 5 and n == 3:
                    if scattering_factor == fit[2] and fit[1] == 'STRUCTURAL':
                        idx = self.parameterFit.index(fit)
                        self.parameterFit.remove(fit)
                        self.currentVal.pop(idx)

        self.setTable()

    def changeStepSize(self):
        self._step_size = self.step_size.text()

    def _setVarMagFromSample(self, sample):

        # setting up the varData and magData dictionaries
        num_layers = len(sample.structure)
        self.magDirection = ['z' for i in range(num_layers)]
        for ele in sample.myelements:
            self.varData[ele] = [[['',''],['',''],['','']] for i in range(num_layers)]
            self.magData[ele] = [[[ele], [''], ['']] for i in range(num_layers)]

        for idx, layer in enumerate(sample.structure):

            for ele in list(layer.keys()):
                element = layer[ele]

                if len(element.polymorph) != 0:
                    self.varData[ele][idx][0] = element.polymorph
                    self.magData[ele][idx][0] = element.polymorph

                    self.varData[ele][idx][1] = element.poly_ratio
                    self.magData[ele][idx][1] = ['' for i in range(len(element.polymorph))]

                    self.varData[ele][idx][2] = element.scattering_factor
                    self.magData[ele][idx][2] = ['' for i in range(len(element.polymorph))]

                # element writing
                if len(element.mag_density) != 0:
                    self.magData[ele][idx][1] = element.mag_density
                    self.magData[ele][idx][2] = element.mag_scattering_factor

                # retrieve the correct magnetization direction
                gamma = layer[ele].gamma
                phi = layer[ele].phi
                if gamma == 90 and phi == 90:
                    self.magDirection[idx] = 'y'
                elif gamma == 0 and phi == 90:
                    self.magDirection[idx] = 'x'
                elif gamma == 0 and phi == 0:
                    self.magDirection[idx] = 'z'


    def _setStructFromSample(self, sample):
        self.sample = sample
        self.structTableInfo = []
        for idx in range(len(sample.structure)):
            # this function will take the upon initialization and load in all parameters
            structInfo = sample.structure
            rows = len(structInfo[idx])
            tempArray = np.empty((rows, 7), dtype=object)
            for row in range(rows):
                for col in range(7):
                    if col == 0:
                        element = list(structInfo[idx].keys())[row]
                        tempArray[row, col] = str(element)
                    elif col == 1:
                        element = list(structInfo[idx].keys())[row]
                        thickness = structInfo[idx][element].thickness
                        tempArray[row, col] = str(thickness)
                    elif col == 2:
                        element = list(structInfo[idx].keys())[row]
                        density = structInfo[idx][element].density
                        tempArray[row, col] = str(density)
                    elif col == 3:
                        element = list(structInfo[idx].keys())[row]
                        roughness = structInfo[idx][element].roughness
                        tempArray[row, col] = str(roughness)
                    elif col == 4:
                        element = list(structInfo[idx].keys())[row]
                        linked_roughness = structInfo[idx][element].linked_roughness
                        tempArray[row, col] = str(linked_roughness)
                    elif col == 5:
                        element = list(structInfo[idx].keys())[row]
                        scattering_factor = structInfo[idx][element].scattering_factor

                        # keeps track of the scattering factors - implementation will change when loading data
                        if type(scattering_factor) is not list and type(scattering_factor) is not np.ndarray:
                            name = 'ff-' + scattering_factor
                            if name not in list(self.eShift.keys()):
                                self.eShift[name] = 0

                        tempArray[row, col] = str(scattering_factor)
                    elif col == 6:
                        element = list(structInfo[idx].keys())[row]
                        stoichiometry = structInfo[idx][element].stoichiometry
                        tempArray[row, col] = str(stoichiometry)

            self.structTableInfo.append(tempArray)


    def setTable(self):
        self.structTable.blockSignals(True)
        idx = self.layerBox.currentIndex()
        tableInfo = self.structTableInfo[idx]
        num_rows = len(tableInfo)

        for col in range(6):
            for row in range(num_rows):
                item = QTableWidgetItem(str(tableInfo[row][col]))
                self.structTable.setItem(row, col, item)
                if col == 0:
                    self.structTable.item(row,col).setBackground(QtGui.QColor('lightGray'))
                if col == 1 and idx == 0:
                    self.structTable.item(row, col).setBackground(QtGui.QColor('lightGray'))

        for fit in self.parameterFit:
            layer = fit[0]
            n = len(fit)
            if layer == idx:  # not scattering factor parameters
                mode = fit[2]
                if mode == 'COMPOUND': # compound mode
                    param = fit[3]
                    param_n = 0

                    if param == 'THICKNESS':
                        param_n = 1
                    elif param == 'DENSITY':
                        param_n = 2
                    elif param == 'ROUGHNESS':
                        param_n = 3
                    elif param == 'LINKED ROUGHNESS':
                        param_n = 4

                    for row in range(num_rows):
                        if param_n != 0:
                            self.structTable.item(row, param_n).setBackground(QtGui.QColor(150, 150, 255))
                elif mode == 'ELEMENT':  # element mode
                    ele = fit[3]
                    param = fit[4]

                    param_n = 0
                    if param == 'THICKNESS':
                        param_n = 1
                    elif param == 'DENSITY':
                        param_n = 2
                    elif param == 'ROUGHNESS':
                        param_n = 3
                    elif param == 'LINKED ROUGHNESS':
                        param_n = 4

                    for row in range(num_rows):
                        my_ele = self.structTableInfo[idx][row][0]
                        if my_ele == ele:
                            if param_n != 0:
                                self.structTable.item(row, param_n).setBackground(QtGui.QColor(150, 255, 150))
            if layer == 'SCATTERING FACTOR':
                for row in range(num_rows):
                    if fit[2] == self.structTableInfo[idx][row][5] and fit[1] == 'STRUCTURAL':
                        self.structTable.item(row,5).setBackground(QtGui.QColor(150, 255, 150))

        self.structTable.blockSignals(False)

    def setTableVar(self):
        self.varTable.blockSignals(True)
        idx = self.sampleInfoLayout.currentIndex()
        layer_idx = self.layerBox.currentIndex()
        ele_idx = self.elementBox.currentIndex()
        # makes sure that when we switch layers we show the same positional element
        if not self.change_elements:
            ele_idx = self.elementBox.currentIndex()
            self.element_index = copy.deepcopy(ele_idx)

        if ele_idx != -1:
            ele = self.structTableInfo[layer_idx][ele_idx][0]

            info = self.varData[ele][layer_idx]

            # Might need to change this implementation
            if len(info[0]) != 0:
                self.varTable.setRowCount(len(info[0]))

                for row in range(len(info[0])):

                    # Element Name
                    item = QTableWidgetItem(info[0][row])
                    self.varTable.setItem(row, 0, item)

                    # Ratio
                    item = QTableWidgetItem(str(info[1][row]))
                    self.varTable.setItem(row, 1, item)
                    # Scattering Factor
                    item = QTableWidgetItem(info[2][row])
                    self.varTable.setItem(row, 2, item)

                    if idx == 1:
                        for fit in self.parameterFit:
                            if len(fit) == 3 and fit[2] == info[2][row] and fit[1] == 'STRUCTURAL':
                                self.varTable.item(row,2).setBackground(QtGui.QColor(150, 255, 150))
                            elif len(fit) == 4:
                                if fit[0] == layer_idx and fit[1] == "POLYMORPHOUS" and fit[2] == ele and fit[3] == info[0][row]:
                                    self.varTable.item(row,1).setBackground(QtGui.QColor(150, 255, 150))
            else:
                for row in range(self.varTable.rowCount()):
                    item = QTableWidgetItem('')
                    self.varTable.setItem(row, 0, item)

                    item = QTableWidgetItem('')
                    self.varTable.setItem(row, 1, item)

                    item = QTableWidgetItem('')
                    self.varTable.setItem(row, 2, item)

        self.varTable.blockSignals(False)

    def setTableMag(self):
        self.magTable.blockSignals(True)
        layer_idx = self.layerBox.currentIndex()

        # set the magnetic direction combobox to the correct magnetization direction
        mag_idx = 0
        dir = self.magDirection[layer_idx]
        if dir=='x':
            mag_idx = 0
        elif dir =='y':
            mag_idx = 1
        elif dir =='z':
            mag_idx = 2

        self.magDirBox.setCurrentIndex(mag_idx)

        labels = []
        density = []
        sf = []
        # Loops through all of the elements
        for ele_idx in range(len(self.structTableInfo[layer_idx])):
            element = self.structTableInfo[layer_idx][ele_idx][0]

            names = self.magData[element][layer_idx][0]
            D = self.magData[element][layer_idx][1]
            S = self.magData[element][layer_idx][2]

            num = len(names)
            if num != 0:
                for i in range(num):
                    labels.append(names[i])
                    if len(D) != 0:
                        density.append(D[i])
                    else:
                        density.append('')
                    if len(S) != 0:
                        sf.append(S[i])
                    else:
                        sf.append('')
            else:
                labels.append(element)
                if len(D) != 0:
                    density.append(D[0])
                else:
                    density.append('')
                if len(S) != 0:
                    sf.append(S[0])
                else:
                    sf.append('')

        row_num = len(labels)
        self.magTable.setRowCount(row_num)
        self.magTable.setVerticalHeaderLabels(
            labels)
        self.magTable.resizeColumnsToContents()  # resizes columns to fit

        for row in range(row_num):
            name = self.magTable.verticalHeaderItem(row).text()  # for color purposes
            mydensity = density[row]
            mysf = sf[row]
            if mydensity == '':
                item = QTableWidgetItem('')
                self.magTable.setItem(row,0,item)
            else:
                item = QTableWidgetItem(str(mydensity))
                self.magTable.setItem(row, 0, item)

                for fit in self.parameterFit:
                    if len(fit) == 3 and fit[0] != 'SCATTERING FACTOR':
                        if layer_idx == fit[0] and fit[1] == 'MAGNETIC' and fit[2] == name:
                            self.magTable.item(row,0).setBackground(QtGui.QColor(150, 255, 150))
                    elif len(fit) == 4 and fit[1] == 'MAGNETIC':
                        if layer_idx == fit[0] and fit[1] == 'MAGNETIC' and fit[3] == name:
                            self.magTable.item(row,0).setBackground(QtGui.QColor(150, 255, 150))
            if mysf == '':
                item = QTableWidgetItem('')
                self.magTable.setItem(row, 1, item)
            else:
                item = QTableWidgetItem(mysf)
                self.magTable.setItem(row, 1, item)
                for fit in self.parameterFit:
                    if len(fit) == 3:
                        scattering_factor = self.magTable.item(row,1).text()
                        if fit[0] == 'SCATTERING FACTOR' and fit[1] == 'MAGNETIC' and fit[2] == scattering_factor:
                            self.magTable.item(row, 1).setBackground(QtGui.QColor(150, 255, 150))

        self.magTable.blockSignals(False)


    def _addLayer(self):

        addLayerApp = compoundInput()
        addLayerApp.show()
        addLayerApp.exec_()
        userinput = addLayerApp.val
        addLayerApp.close()

        # check to see if we have a new element
        num_layers = len(self.structTableInfo)
        my_elements = []
        # checks to see if we have a new element and adds it to the element variation list
        for i in range(len(userinput)):
            # includes new scattering factors into the energy shift
            if userinput[i][5] not in self.sample.eShift.keys():
                self.sample.eShift[userinput[i][5]] = 0
                name = 'ff-' + userinput[i][5]
                self.eShift[name] = 0
            my_elements.append(userinput[i][0])
            if userinput[i][0] not in list(self.varData.keys()):
                self.varData[userinput[i][0]] = [[['',''],['',''],['','']] for j in range(num_layers)]
                self.magData[userinput[i][0]] = [[[userinput[i][0]],[''],['']] for j in range(num_layers)]


        if len(userinput) != 0:
            num = self.layerBox.count()
            idx = self.layerBox.currentIndex()
            if num == 0:
                self.varData = dict()
                self.layerBox.addItem('Substrate')
                self.structTableInfo.insert(0, userinput)
                for i in range(len(userinput)):
                    self.varData[userinput[i][0]] = [['',''],['',''],['','']]
                    self.magData[userinput[i][0]] = [[userinput[i][0]],[''],['']]

            else:
                self.layerBox.addItem('Layer ' + str(num))
                self.structTableInfo.insert(idx+1, userinput)

                for key in list(self.varData.keys()):

                    isVar = False
                    isMag = False
                    data = []
                    data_mag = []
                    # gets the correct elements
                    if key in my_elements:
                        for info_idx in range(len(self.varData[key])):
                            info = self.varData[key][info_idx]
                            info_mag = self.magData[key][info_idx]
                            if not isVar:
                                if info[0][0] != '':

                                    data = [info[0], [1,0], info[2]]
                                    isVar = True
                                else:
                                    data = [['',''],['',''],['','']]

                            if not isMag:
                                if info_mag[1][0] != '' or info_mag[2][0] != '':
                                    n = len(info_mag[0])
                                    data_mag = [info_mag[0],[0 for i in range(n)], info_mag[2] ]
                                    isMag = True
                                else:
                                    data_mag = [[key],[''],['']]
                    else:
                        data = [['',''],['',''],['','']]
                        data_mag = [key,[''],['']]

                    self.varData[key].insert(idx+1, data)
                    self.magData[key].insert(idx+1,data_mag)
                    self.magDirection.insert(idx+1,'z')


    def _removeLayer(self):
        num = self.layerBox.count()  # get the number of layers in the material

        idx = self.layerBox.currentIndex()  # determine which layer has been selected to remove
        if num != 0:
            self.layerBox.removeItem(num-1)

        self.structTableInfo.pop(idx)  # removes the information about that layer

        # removes the element variation data
        for key in list(self.varData.keys()):
            self.varData[key].pop(idx)
            self.magData[key].pop(idx)

        self.setTable()  # sets the table for layer that replaces the removed layer

    def _copyLayer(self):
        num = self.layerBox.count()
        if num == 0:
            self.layerBox.addItem('Substrate')
        else:
            self.layerBox.addItem('Layer ' + str(num))

        idx = self.layerBox.currentIndex()

        newLayer = copy.deepcopy(self.structTableInfo[idx])
        self.structTableInfo.insert(idx+1,newLayer)


        # goes through all the keys
        for key in list(self.varData.keys()):

            info = copy.deepcopy(self.varData[key][idx])
            info_mag = copy.deepcopy(self.magData[key][idx])

            self.varData[key].insert(idx+1,info)
            self.magData[key].insert(idx+1,info_mag)

        new_dir = copy.deepcopy(self.magDirection[idx])
        self.magDirection.insert(idx+1,new_dir)

    def _structural(self):
        self.structButton.setStyleSheet('background: blue; color: white')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.magButton.setStyleSheet('background: lightGrey')
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.sampleInfoLayout.setCurrentIndex(0)

    def _elementVariation(self):
        self.structButton.setStyleSheet('background: lightGrey')
        self.polyButton.setStyleSheet('background: blue; color: white')
        self.magButton.setStyleSheet('background: lightGrey')
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.sampleInfoLayout.setCurrentIndex(1)

    def _magnetic(self):
        self.structButton.setStyleSheet('background: lightGrey')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.magButton.setStyleSheet('background: blue; color: white')
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.sampleInfoLayout.setCurrentIndex(2)

    def _densityprofile(self):

        step_size = float(self.step_size.text())
        self.sample = self._createSample()

        thickness, density, density_magnetic = self.sample.density_profile(step=step_size)

        self.densityWidget.clear()
        self._plotDensityProfile(thickness, density, density_magnetic)


    def _plotDensityProfile(self,thickness, density, density_magnetic):


        num = len(density)
        num = num + len(density_magnetic)

        val = list(density.values())
        mag_val = list(density_magnetic.values())
        check = []
        for key in list(density.keys()):
            if key[-1].isdigit():
                check.append(True)
            else:
                check.append(False)

        for idx in range(len(val)):
            if check[idx]:
                self.densityWidget.plot(thickness,val[idx], pen=pg.mkPen((idx,num),width=2), name=list(density.keys())[idx])
            else:
                self.densityWidget.plot(thickness, val[idx], pen=pg.mkPen((idx,num),width=2), name=list(density.keys())[idx])

        for idx in range(len(mag_val)):
            myname = 'Mag: ' + list(density_magnetic.keys())[idx]
            self.densityWidget.plot(thickness, -mag_val[idx], pen=pg.mkPen((num-idx,num),width=2, style=Qt.DashLine), name=myname)
        self.densityWidget.setLabel('left', "Density (mol/cm^3)")
        self.densityWidget.setLabel('bottom', "Thickness (Å)")
    def _createSample(self):
        # This function takes the information from the tables and converts it into a usable form

        m = len(self.structTableInfo)  # determines how many layers in the sample
        sample = ms.slab(m)  # initializes the slab class

        # loops through each layer and sets the appropriate parameters
        for idx in range(m):
            formula = ''  # used to determine the chemical formula
            thickness = []  # thickness of the layer
            density = []  # density of the layer
            roughness = []  # roughness of the layer
            linked_roughness = []  # linked roughness of the layer

            layer = self.structTableInfo[idx]  # gets the layer information

            for ele in range(len(layer)):
                element = layer[ele]  # element data

                # recreates the chemical formula string
                formula = formula + element[0]
                if element[6] != '1':
                    formula = formula + element[6]

                thickness.append(float(element[1]))  # gets thickness data
                density.append(float(element[2]))  # gets density data
                roughness.append(float(element[3]))  # gets roughness data
                if element[4] != '' and element[4] != 'False' and element[4] != False:
                    linked_roughness.append(float(element[4]))
                else:
                    linked_roughness.append(False)

                # still need to take into account sf that are different than element
            sample.addlayer(idx,formula,thickness,density=density, roughness=roughness, linked_roughness=linked_roughness)

        for key, value in self.eShift.items():
            if str(key).startswith('ff-'):
                sample.eShift[str(key)[3:]] = value
            elif str(key).startswith('ffm-'):
                sample.mag_eShift[str(key)[4:]] = value

        for idx in range(m):
            layer = self.structTableInfo[idx]  # gets the layer information
            for ele in range(len(layer)):
                ele_name = layer[ele][0]

                poly = self.varData[ele_name][idx]  # retrieves the element variation data for particular layer

                names = poly[0]
                ratio = poly[1]
                scattering_factor = poly[2]

                if len(names) > 1:
                    if names[0] != '':
                        ratio = [float(ratio[i]) for i in range(len(ratio))]
                        sample.polymorphous(idx,ele_name,names,ratio,sf=scattering_factor)

        for idx in range(m):
            layer = self.structTableInfo[idx]  # gets the layer information
            for ele in range(len(layer)):
                ele_name = layer[ele][0]

                mag = self.magData[ele_name][idx]


                names = mag[0]
                ratio = mag[1]
                scattering_factor = mag[2]

                if ratio[0] != '' and names[0] != '' and scattering_factor[0] != '':
                    ratio = [float(ratio[i]) for i in range(len(ratio))]
                    sample.magnetization(idx,names,ratio,scattering_factor)


        return sample


    def getData(self):

        # this function retrieves all of the important data from the sample
        for j in range(len(self.sample.structure)):
            layer = self.sample.structure[j]
            elekeys = list(layer.keys())

            for ele in elekeys:

                if len(layer[ele].polymorph) != 0:
                    mag_density = ['' for i in range(len(layer[ele].polymorph))]
                    mag_sf = ['' for i in range(len(layer[ele].polymorph))]
                    self.varData[ele][j] = [layer[ele].polymorph, list(layer[ele].poly_ratio),
                                            layer[ele].scattering_factor]

                    if len(layer[ele].mag_density) != 0:
                        mag_density = list(layer[ele].mag_density)
                    if len(layer[ele].mag_scattering_factor) != 0:
                        mag_sf = layer[ele].mag_scattering_factor

                    self.magData[ele][j] = [layer[ele].polymorph, mag_density, mag_sf]


                else:
                    mag_density = ['']
                    mag_sf = ['']
                    if len(layer[ele].mag_density) != 0:
                        mag_density = list(layer[ele].mag_density)
                    if len(layer[ele].mag_scattering_factor) != 0:
                        mag_sf = layer[ele].mag_scattering_factor

                        # implementation will change
                        for scat in mag_sf:
                            name = 'ffm-'+scat
                            if name not in list(self.eShift.keys()):
                                self.eShift[name] = 0

                    self.magData[ele][j] = [[ele],mag_density, mag_sf]

            # gets the magnetic direction for that particular layer
            gamma = layer[ele].gamma
            phi = layer[ele].phi
            if gamma == 90 and phi == 90:
                self.magDirection[j] = 'y'
            elif gamma == 0 and phi == 90:
                self.magDirection[j] = 'x'
            elif gamma == 0 and phi == 0:
                self.magDirection[j] = 'z'


class reflectivityWidget(QWidget):
    def __init__(self, sWidget, data, data_dict, sim_dict):
        super().__init__()

        self.rom = [True, False, False]
        self.romSim = [True, False, False]

        self.scanType = True  # if true do reflectivity scan

        self.bs = dict()  # background shift
        self.sf = dict()  # scaling factor

        self.sfBsFitParams = []  # variable used to keep track of the fitting parameters
        self.currentVal = []

        self.sWidget = sWidget
        self.sample = sWidget.sample
        self.data = data
        self.data_dict = data_dict

        self.fit = []  # scans to be fit
        self.bounds = []  # bounds
        self.weights = []  # weights of the bounds

        self.isChangeTable = False
        self.axis_state = True
        self.scan_state = True
        self.previousIdx = 0

        # Adding the plotting Widget
        self.spectrumWidget = pg.PlotWidget()
        self.spectrumWidget.setBackground('w')

        self.spectrumWidget.addLegend()
        # This will be used to determine which scan to view
        whichScanLabel = QLabel('Scans:')
        whichScanLabel.setFixedWidth(30)
        self.whichScan = QComboBox()
        self.whichScan.setFixedWidth(250)
        for scan in data:
            self.whichScan.addItem(scan[2])

        self.whichScan.currentIndexChanged.connect(self.myPlotting)

        self.whichScan.setCurrentIndex(0)
        whichScanLayout = QHBoxLayout()

        whichScanLayout.addWidget(whichScanLabel)
        whichScanLayout.addWidget(self.whichScan)

        self.plot_scans()

        self.fitButton = QPushButton('Fit Scan')
        self.fitButton.clicked.connect(self._scanSelection)


        # create the widgets for doing a simulation only
        simulationLabel = QLabel('Simulations')

        simEnergyLayout = QHBoxLayout()
        simEnergyLabel = QLabel('Energy (eV): ')
        simEnergyLabel.setFixedWidth(80)
        self.simEnergyEdit = QLineEdit()
        self.simEnergyEdit.setText('500')
        self.simEnergyEdit.editingFinished.connect(self.setEnergy)
        self.simEnergyEdit.setFixedWidth(200)
        simEnergyLayout.addWidget(simEnergyLabel)
        simEnergyLayout.addWidget(self.simEnergyEdit)

        simAngleLayout = QHBoxLayout()
        simAngleLabel = QLabel('Angle (degrees): ')
        simAngleLabel.setFixedWidth(80)
        self.simAngleEdit = QLineEdit()
        self.simAngleEdit.setText('5')
        self.simAngleEdit.editingFinished.connect(self.setAngle)
        self.simAngleEdit.setFixedWidth(200)
        simAngleLayout.addWidget(simAngleLabel)
        simAngleLayout.addWidget(self.simAngleEdit)

        energyRangeLayout = QHBoxLayout()
        energyRangeLabel = QLabel('Energy Range: ')
        energyRangeLabel.setFixedWidth(100)
        self.simLowEnergy = QLineEdit()  # energy ranges
        self.simLowEnergy.setText('450')
        self.simLowEnergy.editingFinished.connect(self.setLowerEnergy)
        self.simLowEnergy.setFixedWidth(90)
        self.simUpEnergy = QLineEdit()
        self.simUpEnergy.setText('550')
        self.simUpEnergy.editingFinished.connect(self.setUpperEnergy)
        self.simUpEnergy.setFixedWidth(90)
        energyRangeLayout.addWidget(energyRangeLabel)
        energyRangeLayout.addWidget(self.simLowEnergy)
        energyRangeLayout.addWidget(self.simUpEnergy)

        self.polarBox = QComboBox()
        self.polarBox.addItems(['s-polarized', 'p-polarized', 'left circular', 'right circular',
                                'linear asymmetry', 'circular asymmetry'])
        self.polarBox.currentIndexChanged.connect(self.mySimPlotting)

        self.reflectButton = QPushButton('Reflectivity Scan')
        self.reflectButton.clicked.connect(self.changeToReflScan)
        self.reflectButton.setStyleSheet('background: cyan')
        self.reflectButton.setFixedWidth(300)
        self.energyButton = QPushButton('Energy Scan')
        self.energyButton.clicked.connect(self.changeToEnergyScan)
        self.energyButton.setStyleSheet('background: grey')
        self.energyButton.setFixedWidth(300)

        self.rButtonSim = QPushButton('Reflectometry')
        self.rButtonSim.setStyleSheet('background: grey')
        self.rButtonSim.clicked.connect(self.rPlotSim)
        self.opButtonSim = QPushButton('Optical')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.clicked.connect(self.opPlotSim)
        self.opmButtonSim = QPushButton('Magneto-Optical')
        self.opmButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.clicked.connect(self.opmPlotSim)

        simButtonLayout = QHBoxLayout()
        simButtonLayout.addWidget(self.rButtonSim)
        simButtonLayout.addWidget(self.opButtonSim)
        simButtonLayout.addWidget(self.opmButtonSim)
        # continue on with stuff
        self.boundWeightTable = QTableWidget()

        self.selectedScans = QComboBox()

        self.selectedScans.activated.connect(self.changeColorFit)
        self.selectedScans.activated.connect(self.readTable)
        self.whichScan.activated.connect(self.changeColorScan)
        self.selectedScans.activated.connect(self.myPlotting)
        self.selectedScans.activated.connect(self.setTable)
        self.selectedScans.activated.connect(self.changeFitColor)

        self.addBoundaryButton = QPushButton('Add Boundary')
        self.addBoundaryButton.setFixedWidth(250)
        self.removeBoundaryButton = QPushButton('Remove Boundary')
        self.removeBoundaryButton.setMaximumWidth(250)
        self.removeScan = QPushButton('Remove Scan')
        self.removeScan.setMaximumWidth(250)
        self.removeScan.clicked.connect(self._removeScanSelection)
        self.removeScan.setStyleSheet("background-color : cyan")

        self.qz = QRadioButton('qz (A)', self)
        self.qz.setChecked(True)
        self.qz.toggled.connect(self.updateAxis)
        self.angle = QRadioButton('Theta (degrees)', self)
        self.angle.toggled.connect(self.updateAxis)
        axis_label = QLabel('Reflectivity Axis: ')
        hbox = QHBoxLayout()

        self.rButton = QPushButton('Reflectometry')
        self.rButton.setStyleSheet('background: cyan')
        self.rButton.clicked.connect(self.rPlot)
        self.opButton = QPushButton('Optical')
        self.opButton.setStyleSheet('background: grey')
        self.opButton.clicked.connect(self.opPlot)
        self.opmButton = QPushButton('Magneto-Optical')
        self.opmButton.setStyleSheet('background: grey')
        self.opmButton.clicked.connect(self.opmPlot)

        self.stepWidget = QLineEdit()
        self.stepWidget.textChanged.connect(self.changeStepSize)
        self.stepWidget.setText(self.sWidget._step_size)
        self.stepWidget.setMaximumWidth(175)
        stepLabel = QLabel('Step Size (Å):')
        stepLabel.setMaximumWidth(65)
        stepLayout = QHBoxLayout()
        stepLayout.addWidget(stepLabel)
        stepLayout.addWidget(self.stepWidget)

        hbox.addWidget(axis_label)
        hbox.addWidget(self.qz)
        hbox.addWidget(self.angle)

        selectedScansLabel = QLabel('Fit Scans: ')
        selectedScansLabel.setMaximumWidth(60)
        selectedScansLayout = QHBoxLayout()
        selectedScansLayout.addWidget(selectedScansLabel)
        selectedScansLayout.addWidget(self.selectedScans)

        buttonLayout = QHBoxLayout()
        buttonLayout.addWidget(self.rButton)
        buttonLayout.addWidget(self.opButton)
        buttonLayout.addWidget(self.opmButton)

        scanSelectionLayout = QVBoxLayout()
        scanSelectionLayout.addLayout(whichScanLayout)
        scanSelectionLayout.addWidget(self.fitButton)
        scanSelectionLayout.addLayout(buttonLayout)
        scanSelectionLayout.addStretch(1)
        scanSelectionLayout.addLayout(hbox)
        scanSelectionLayout.addLayout(stepLayout)
        scanSelectionLayout.addStretch(1)
        scanSelectionLayout.addWidget(simulationLabel)
        scanSelectionLayout.addLayout(simButtonLayout)
        scanSelectionLayout.addLayout(simAngleLayout)
        scanSelectionLayout.addLayout(simEnergyLayout)
        scanSelectionLayout.addLayout(energyRangeLayout)
        scanSelectionLayout.addWidget(self.polarBox)
        scanSelectionLayout.addWidget(self.reflectButton)
        scanSelectionLayout.addWidget(self.energyButton)


        #scanSelectionLayout.addStretch(1)
        #scanSelectionLayout.addLayout(selectedScansLayout)


        toplayout = QHBoxLayout()
        toplayout.addWidget(self.spectrumWidget)
        toplayout.addLayout(scanSelectionLayout)


        boundWidget = QWidget()
        boundLayout = QVBoxLayout()
        boundLayout.addLayout(selectedScansLayout)
        self.addBoundaryButton.clicked.connect(self.addBoundWeight)
        boundLayout.addWidget(self.addBoundaryButton)
        self.removeBoundaryButton.clicked.connect(self.removeBoundWeight)
        boundLayout.addWidget(self.removeBoundaryButton)
        boundLayout.addWidget(self.removeScan)

        boundWidget.setLayout(boundLayout)
        allScansLayout = QHBoxLayout()
        allScanLabel = QLabel('All Scans: ')
        allScanLabel.setFixedWidth(65)
        self.allScan = QCheckBox()
        self.allScan.setChecked(1)
        self.allScan.stateChanged.connect(self.allScanStateChanged)
        allScansLayout.addWidget(allScanLabel)
        allScansLayout.addWidget(self.allScan)

        sfbsLayout = QVBoxLayout()

        bsLayout = QHBoxLayout()
        bsLabel = QLabel('Background Shift:')
        bsLabel.setFixedWidth(90)
        self.backgroundShift = QLineEdit()
        self.scalingFactor = QLineEdit()
        self.backgroundShift.textChanged.connect(self.changeSFandBS)
        self.backgroundShift.setFixedWidth(100)
        self.backgroundShift.setFixedHeight(25)
        bsLayout.addWidget(bsLabel)
        bsLayout.addWidget(self.backgroundShift)

        sfLayout = QHBoxLayout()
        sfLabel = QLabel('Scaling Factor:')
        sfLabel.setFixedWidth(90)
        self.scalingFactor.textChanged.connect(self.changeSFandBS)
        self.scalingFactor.setFixedWidth(100)
        self.scalingFactor.setFixedHeight(25)
        sfLayout.addWidget(sfLabel)
        sfLayout.addWidget(self.scalingFactor)

        self.bsFit = QPushButton('Fit')
        self.bsFit.clicked.connect(self.fitBackgroundShift)
        self.bsUnfit = QPushButton('Remove Fit')
        self.bsUnfit.clicked.connect(self.unfitBackgroundShift)
        bsFitLayout = QHBoxLayout()
        bsFitLayout.addWidget(self.bsFit)
        bsFitLayout.addWidget(self.bsUnfit)

        self.sfFit = QPushButton('Fit')
        self.sfFit.clicked.connect(self.fitScalingFactor)
        self.sfUnfit = QPushButton('Remove Fit')
        self.sfUnfit.clicked.connect(self.unfitScalingFactor)
        sfFitLayout = QHBoxLayout()
        sfFitLayout.addWidget(self.sfFit)
        sfFitLayout.addWidget(self.sfUnfit)

        sfbsLayout.addLayout(allScansLayout)
        sfbsLayout.addLayout(bsLayout)
        sfbsLayout.addLayout(bsFitLayout)
        sfbsLayout.addLayout(sfLayout)
        sfbsLayout.addLayout(sfFitLayout)

        semiBotLayout = QHBoxLayout()
        semiBotLayout.addLayout(sfbsLayout)
        semiBotLayout.addWidget(self.boundWeightTable)
        self.boundWeightTable.itemChanged.connect(self.value_changed)

        bottomlayout = QHBoxLayout()
        bottomlayout.addLayout(semiBotLayout)
        bottomlayout.addWidget(boundWidget)

        pagelayout = QVBoxLayout()
        pagelayout.addLayout(toplayout)
        pagelayout.addLayout(bottomlayout)

        self.backgroundShift.editingFinished.connect(self.bsChange)
        self.scalingFactor.editingFinished.connect(self.sfChange)
        self.setLayout(pagelayout)

    def changeToReflScan(self):
        self.energyButton.setStyleSheet('background: grey')
        self.reflectButton.setStyleSheet('background: cyan')
        self.scanType = True
        self.mySimPlotting()

    def changeToEnergyScan(self):
        self.energyButton.setStyleSheet('background: cyan')
        self.reflectButton.setStyleSheet('background: grey')
        self.scanType = False
        self.mySimPlotting()

    def rPlotSim(self):
        self.romSim = [True, False, False]
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: cyan')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        self.spectrumWidget.clear()
        pol = self.polarBox.currentText()
        step_size = float(self.sWidget._step_size)
        polName = 'S'
        if pol == 's-polarized':
            polName = 'S'
        elif pol == 'p-polarized':
            polName = 'P'
        elif pol == 'left circular':
            polName = 'LC'
        elif pol == 'right circular':
            polName = 'RC'
        elif pol == 'linear asymmetry':
            polName = 'AL'
        else:
            polName = 'AC'

        # number of points
        n = 1001
        energy = float(self.simEnergyEdit.text())
        if self.scanType:
            Theta = np.linspace(0.1,89.1,num=n)
            qz = np.sin(Theta*np.pi/180)*(energy * 0.001013546143)
            qz, R = self.sample.reflectivity(energy, qz, s_min=step_size)
            R = R[polName]
            if self.axis_state:  # momentum transfer
                self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 1),width=2), name='Simulation')
                self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
            else:  # angle
                self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0,1), width=2), name='Simulation')
                self.spectrumWidget.setLabel('bottom', "Angle (degrees)")

            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            if polName == 'S' or polName == 'P' or polName == 'LC' or polName == 'RC':
                self.spectrumWidget.setLogMode(False, True)
            else:
                self.spectrumWidget.setLogMode(False, False)

        else:
            angle = float(self.simAngleEdit.text())
            lw = float(self.simLowEnergy.text())
            up = float(self.simUpEnergy.text())
            E = np.linspace(lw,up, num=n)
            E, R = self.sample.energy_scan(angle, E, s_min=step_size)
            R = R[polName]

            self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 1),width=2), name='Simulation')
            self.spectrumWidget.setLabel('bottom', "Energy (eV)")
            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            self.spectrumWidget.setLogMode(False, False)


    def opPlotSim(self):
        self.spectrumWidget.clear()
        self.romSim = [False, True, False]
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: cyan')
        self.opmButtonSim.setStyleSheet('background: grey')

        E = float(self.simEnergyEdit.text())

        step_size = float(self.sWidget._step_size)
        thickness, density, density_magnetic = self.sample.density_profile(step=step_size)

        sf = dict()  # form factors of non-magnetic components
        sfm = dict()  # form factors of magnetic components

        # Non-Magnetic Scattering Factor
        for e in self.sample.find_sf[0].keys():
            name = 'ff-' + self.sample.find_sf[0][e]
            dE = float(self.sWidget.eShift[name])
            sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E + dE, False)
        # Magnetic Scattering Factor
        for em in self.sample.find_sf[1].keys():
            name = 'ffm-' + self.sample.find_sf[1][em]
            dE = float(self.sWidget.eShift[name])
            sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E + dE, True)

        delta, beta = ms.index_of_refraction(density, sf, E)  # calculates dielectric constant for structural component

        self.spectrumWidget.plot(thickness, delta, pen=pg.mkPen((0, 2), width=2), name='delta')
        self.spectrumWidget.plot(thickness, beta, pen=pg.mkPen((1, 2), width=2), name='beta')
        self.spectrumWidget.setLabel('left', "Reflectivity, R")
        self.spectrumWidget.setLabel('bottom', "Thickness, Å")
        self.spectrumWidget.setLogMode(False, False)

    def opmPlotSim(self):

        self.spectrumWidget.clear()
        self.romSim = [False, False, True]
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: cyan')

        E = float(self.simEnergyEdit.text())

        step_size = float(self.sWidget._step_size)
        thickness, density, density_magnetic = self.sample.density_profile(step=step_size)

        sf = dict()  # form factors of non-magnetic components
        sfm = dict()  # form factors of magnetic components

        # Non-Magnetic Scattering Factor
        for e in self.sample.find_sf[0].keys():
            name = 'ff-' + self.sample.find_sf[0][e]
            dE = float(self.sWidget.eShift[name])
            sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E + dE, False)
        # Magnetic Scattering Factor
        for em in self.sample.find_sf[1].keys():
            name = 'ffm-' + self.sample.find_sf[1][em]
            dE = float(self.sWidget.eShift[name])
            sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E + dE, True)

        delta_m, beta_m = ms.magnetic_optical_constant(density_magnetic, sfm, E)

        if len(density_magnetic) == 0:
            self.spectrumWidget.plot(thickness, np.zeros(len(thickness)), pen=pg.mkPen((0, 2), width=2), name='delta_m')
            self.spectrumWidget.plot(thickness, np.zeros(len(thickness)), pen=pg.mkPen((1, 2), width=2), name='beta_m')
        else:
            self.spectrumWidget.plot(thickness, delta_m, pen=pg.mkPen((0, 2), width=2), name='delta_m')
            self.spectrumWidget.plot(thickness, beta_m, pen=pg.mkPen((1, 2), width=2), name='beta_m')

        self.spectrumWidget.setLabel('left', "Reflectivity, R")
        self.spectrumWidget.setLabel('bottom', "Thickness, Å")
        self.spectrumWidget.setLogMode(False, False)

    def setAngle(self):
        angle = self.simAngleEdit.text()

        if angle != '' and angle != '-':
            angle = float(angle)
            if angle < 0 or angle == 0:
                angle = 0.1
                self.simAngleEdit.setText(str(angle))
                self.mySimPlotting()
            elif angle > 90 or angle == 90:
                angle = 89.9
                self.simAngleEdit.setText(str(angle))
                self.mySimPlotting()


    def setEnergy(self):
        energy = self.simEnergyEdit.text()

        if energy != '':
            energy = float(energy)
            lower = str(energy - 50)
            upper = str(energy + 50)
            self.simLowEnergy.setText(lower)
            self.simUpEnergy.setText(upper)
            self.mySimPlotting()

    def setLowerEnergy(self):
        energy = float(self.simEnergyEdit.text())
        lower = float(self.simLowEnergy.text())

        if energy < lower:
            lower = energy - 50

        self.simLowEnergy.setText(str(lower))
        self.mySimPlotting()

    def setUpperEnergy(self):
        energy = float(self.simEnergyEdit.text())
        upper = float(self.simUpEnergy.text())

        if energy > upper:
            upper = energy + 50

        self.simUpEnergy.setText(str(upper))
        self.mySimPlotting()
    def value_changed(self):
        if not self.isChangeTable:
            idx = self.selectedScans.currentIndex()
            row = self.boundWeightTable.currentRow()
            column = self.boundWeightTable.currentColumn()
            item = self.boundWeightTable.currentItem().text()
            if row == 0 or row == 1: # upper or lower bounds
                self.bounds[idx][column][row] = item
            elif row == 2:  # weights
                self.weights[idx][column] = item

    def bsChange(self):
        name = self.selectedScans.currentText()
        idx = self.selectedScans.currentIndex()
        var = self.backgroundShift.text()
        copy_of_fit = copy.deepcopy(self.sfBsFitParams)
        old_var = ''
        if name != '':
            # first change the value
            if self.allScan.checkState() == 0:
                old_var = self.bs[name]
                self.bs[name] = var
            else:
                old_var = self.bs[name]
                for key in list(self.bs.keys()):
                    self.bs[key] = var

            for fit in copy_of_fit:
                if fit[0] == 'BACKGROUND SHIFT':
                    upper = str(float(var) + 5e-8)
                    lower = str(float(var) - 5e-8)
                    if self.allScan.checkState() != 0:
                        idx = self.sfBsFitParams.index(['BACKGROUND SHIFT', 'ALL SCANS'])
                        self.currentVal[idx] = [var, [lower, upper]]
                    else:
                        idx = self.sfBsFitParams.index(['BACKGROUND SHIFT', name])
                        self.currentVal[idx] = [var, [lower, upper]]


    def sfChange(self):
        name = self.selectedScans.currentText()
        idx = self.selectedScans.currentIndex()
        var = self.scalingFactor.text()
        copy_of_fit = copy.deepcopy(self.sfBsFitParams)

        if name != '':
            # first change the value
            if self.allScan.checkState() == 0:
                self.sf[name] = var
            else:
                for key in list(self.sf.keys()):
                    self.sf[key] = var

            for fit in copy_of_fit:
                if fit[0] == 'SCALING FACTOR':
                    upper = str(float(var) + 0.2)
                    lower = str(float(var) - 0.2)
                    if self.allScan.checkState() != 0:
                        idx = self.sfBsFitParams.index(['SCALING FACTOR', 'ALL SCANS'])
                        self.currentVal[idx] = [var, [lower, upper]]
                    else:
                        idx = self.sfBsFitParams.index(['SCALING FACTOR', name])
                        self.currentVal[idx] = [var, [lower, upper]]

    def changeFitColor(self):
        name = self.selectedScans.currentText()
        self.backgroundShift.setStyleSheet('background: white')
        self.scalingFactor.setStyleSheet('background: white')
        for fit in self.sfBsFitParams:

            if fit[0] == 'BACKGROUND SHIFT' and fit[1] == 'ALL SCANS':
                self.backgroundShift.setStyleSheet('background: red')
            elif fit[0] == 'BACKGROUND SHIFT' and fit[1] == name:
                self.backgroundShift.setStyleSheet('background: red')


            if fit[0] == 'SCALING FACTOR' and fit[1] == 'ALL SCANS':
                self.scalingFactor.setStyleSheet('background: red')
            elif fit[0] == 'SCALING FACTOR' and fit[1] == name:
                self.scalingFactor.setStyleSheet('background: red')


    def allScanStateChanged(self):
        # erase the fitting parameters previously input by the user
        self.sfBsFitParams = []
        self.changeFitColor()

    def fitBackgroundShift(self):
        state = self.allScan.checkState()
        name = self.selectedScans.currentText()
        value = self.backgroundShift.text()
        if name != '':
            if state == 2:
                fit = ['BACKGROUND SHIFT', 'ALL SCANS']
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 5e-8
                    upper = float(value) + 5e-8
                    self.currentVal.append([float(value), [lower, upper]])
            else:
                fit = ['BACKGROUND SHIFT', name]
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 5e-8
                    upper = float(value) + 5e-8
                    self.currentVal.append([float(value), [lower, upper]])


        self.changeFitColor()

    def fitScalingFactor(self):
        state = self.allScan.checkState()
        name = self.selectedScans.currentText()
        value = self.scalingFactor.text()
        if name != '':
            if state == 2:
                fit = ['SCALING FACTOR', 'ALL SCANS']
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 0.2
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 0.2
                    self.currentVal.append([float(value), [lower, upper]])
            else:
                fit = ['SCALING FACTOR', name]
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 0.2
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 0.2
                    self.currentVal.append([float(value), [lower, upper]])

        self.changeFitColor()

    def unfitBackgroundShift(self):
        state = self.allScan.checkState()
        name = self.selectedScans.currentText()
        if name != '':
            if state == 2:
                if ['BACKGROUND SHIFT', 'ALL SCANS'] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['BACKGROUND SHIFT', 'ALL SCANS'])
                    self.sfBsFitParams.remove(['BACKGROUND SHIFT', 'ALL SCANS'])
                    self.currentVal.pop(idx)
            else:
                if ['BACKGROUND SHIFT', name] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['BACKGROUND SHIFT', name])
                    self.sfBsFitParams.remove(['BACKGROUND SHIFT', name])
                    self.currentVal.pop(idx)
        self.changeFitColor()

    def unfitScalingFactor(self):
        state = self.allScan.checkState()
        name = self.selectedScans.currentText()
        if name != '':
            if state == 2:
                if ['SCALING FACTOR', 'ALL SCANS'] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['SCALING FACTOR', 'ALL SCANS'])
                    self.sfBsFitParams.remove(['SCALING FACTOR', 'ALL SCANS'])
                    self.currentVal.pop(idx)
            else:
                if ['SCALING FACTOR', name] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['SCALING FACTOR', name])
                    self.sfBsFitParams.remove(['SCALING FACTOR', name])
                    self.currentVal.pop(idx)

        self.changeFitColor()

    def changeSFandBS(self):
        idx = self.selectedScans.currentIndex()
        name = self.selectedScans.currentText()
        bs = self.backgroundShift.text()
        sf = self.scalingFactor.text()
        if bs != '' and sf != '':
            if self.allScan.checkState() == 0:  # case where all scans have different bs and sf
                self.bs[name] = bs
                self.sf[name] = sf
            else: # case where all scans have same bs and sf
                for key in list(self.bs.keys()):
                    self.bs[key] = bs
                    self.sf[key] = sf


    def changeStepSize(self):
        self.sWidget._step_size = self.stepWidget.text()

    def mySimPlotting(self):
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        # self.sample = self.sWidget.sample
        idx = self.romSim.index(True)
        if idx == 0:
            self.rPlotSim()
        elif idx == 1:
            self.opPlotSim()
        elif idx == 2:
            self.opmPlotSim()

    def myPlotting(self):

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        #self.sample = self.sWidget.sample
        idx = self.rom.index(True)
        if idx == 0:
            self.rPlot()
        elif idx == 1:
            self.opPlot()
        elif idx == 2:
            self.opmPlot()

    def rPlot(self):
        self.rom = [True, False, False]

        self.rButton.setStyleSheet('background: cyan')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        if self.scan_state:
            self.plot_scans()
        else:
            self.plot_selected_scans()
            self.setTable()

    def opPlot(self):
        self.rom = [False,True, False]
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: cyan')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        name = ''
        self.spectrumWidget.clear()

        if self.scan_state:
            name = self.whichScan.currentText()
        else:
            name = self.selectedScans.currentText()
        if name != '':
            E = self.data_dict[name]['Energy']

            step_size = float(self.sWidget._step_size)
            thickness, density, density_magnetic = self.sample.density_profile(step=step_size)  # Computes the density profile

            sf = dict()  # form factors of non-magnetic components
            sfm = dict()  # form factors of magnetic components

            # Non-Magnetic Scattering Factor
            for e in self.sample.find_sf[0].keys():
                name = 'ff-' + self.sample.find_sf[0][e]
                dE = float(self.sWidget.eShift[name])
                sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E+dE, False)
            # Magnetic Scattering Factor
            for em in self.sample.find_sf[1].keys():
                name = 'ffm-' + self.sample.find_sf[1][em]
                dE = float(self.sWidget.eShift[name])
                sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E+dE, True)

            delta, beta = ms.index_of_refraction(density, sf, E)  # calculates dielectric constant for structural component

            self.spectrumWidget.plot(thickness, delta, pen=pg.mkPen((0, 2), width=2), name='delta')
            self.spectrumWidget.plot(thickness, beta, pen=pg.mkPen((1, 2), width=2), name='beta')
            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            self.spectrumWidget.setLabel('bottom', "Thickness, Å")
            self.spectrumWidget.setLogMode(False, False)
            #delta_m, beta_m = ms.magnetic_optical_constant(density_magnetic, sfm, E)  # calculates dielectric constant for magnetic component

    def opmPlot(self):
        self.rom = [False, False, True]

        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: cyan')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        name = ''
        self.spectrumWidget.clear()

        if self.scan_state:
            name = self.whichScan.currentText()
        else:
            name = self.selectedScans.currentText()

        if name != '':
            E = self.data_dict[name]['Energy']

            step_size = float(self.sWidget._step_size)
            thickness, density, density_magnetic = self.sample.density_profile(step=step_size)  # Computes the density profile

            sf = dict()  # form factors of non-magnetic components
            sfm = dict()  # form factors of magnetic components

            # Non-Magnetic Scattering Factor

            for em in self.sample.find_sf[1].keys():
                sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E, True)

            delta_m, beta_m = ms.magnetic_optical_constant(density_magnetic, sfm, E)  # calculates dielectric constant for magnetic component
            self.spectrumWidget.plot(thickness, delta_m, pen=pg.mkPen((0, 2), width=2), name='delta_m')
            self.spectrumWidget.plot(thickness, beta_m, pen=pg.mkPen((1, 2), width=2), name='beta_m')
            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            self.spectrumWidget.setLabel('bottom', "Thickness, Å")
            self.spectrumWidget.setLogMode(False, False)

    def updateAxis(self):
        self.readTable()
        rbtn = self.sender()

        if rbtn.isChecked() == True:
            if rbtn.text() == 'qz (A)':
                self.axis_state = True
            else:
                self.axis_state = False

        self.myPlotting()
        self.mySimPlotting()


    def changeColorScan(self):
        self.selectedScans.setStyleSheet('background: white; selection-background-color: grey')
        self.whichScan.setStyleSheet('background: red; selection-background-color: red')
        self.scan_state = True

    def changeColorFit(self):
        self.selectedScans.setStyleSheet('background: red; selection-background-color: red')
        self.whichScan.setStyleSheet('background: white; selection-background-color: grey')
        self.scan_state = False

    def setTable(self):

        self.isChangeTable = True
        idx = self.selectedScans.currentIndex()
        name = self.selectedScans.currentText()
        self.scalingFactor.setText(self.sf[name]) # setting the appropriate scaling factor

        self.backgroundShift.setText(self.bs[name])  # setting the appropriate background shift
        scan = self.fit[idx]

        if scan != '':

            E = self.data_dict[scan]['Energy']
            mykeys = list(self.data_dict[scan].keys())

            bound = self.bounds[idx]
            weight = self.weights[idx]
            col = len(bound)
            row = 3

            self.boundWeightTable.setRowCount(row)
            self.boundWeightTable.setColumnCount(col)

            self.boundWeightTable.setVerticalHeaderLabels(['Lower Bound', 'Upper Bound', 'Weight'])

            for i in range(row):
                for j in range(col):
                    if i == 0:
                        myitem = ''
                        if 'Angle' not in mykeys:
                            if not self.axis_state:
                                if len(bound[j][0]) != 0:
                                    myitem = str(np.arcsin(float(bound[j][0])/(E * 0.001013546143))* 180 / np.pi)[:7]
                            else:
                                myitem = copy.copy(bound[j][0])
                        else:
                            myitem = copy.copy(bound[j][0])


                        item = QTableWidgetItem(myitem)
                        self.boundWeightTable.setItem(i,j,item)

                        #self.boundWeightTable.setItem(i, j, item)
                    elif i == 1:
                        myitem = ''
                        if 'Angle' not in mykeys:
                            if not self.axis_state:
                                if len(bound[j][1]) != 0:
                                    myitem = str(np.arcsin(float(bound[j][1])/(E * 0.001013546143)) * 180 / np.pi)[:7]
                            else:
                                myitem = copy.copy(bound[j][1])
                        else:
                            myitem = copy.copy(bound[j][1])

                        item = QTableWidgetItem(myitem)
                        self.boundWeightTable.setItem(i, j, item)
                    elif i == 2:

                        item = QTableWidgetItem(weight[j])
                        self.boundWeightTable.setItem(i, j, item)
        self.isChangeTable = False

    def _scanSelection(self):
        idx = self.whichScan.currentIndex()
        name = self.whichScan.currentText()

        if name not in self.fit:
            if 'Angle' in list(self.data_dict[name].keys()):
                lower = str(self.data_dict[name]['Data'][3][0])[0:7]
                upper = str(self.data_dict[name]['Data'][3][-1])[0:7]
            else:
                lower = str(self.data_dict[name]['Data'][0][0])[0:7]
                upper = str(self.data_dict[name]['Data'][0][-1])[0:7]

            self.fit.append(name)
            self.bounds.append([[lower,upper]])
            self.weights.append(['1'])
            self.selectedScans.addItem(name)

            # sets the bs and sf values
            if self.allScan.checkState() == 2: # all scans state checked
                if len(self.bs) == 0:
                    self.bs[name] = str(self.data_dict[name]['Background Shift'])
                    self.sf[name] = str(self.data_dict[name]['Scaling Factor'])
                else:
                    key = list(self.bs.keys())[0]
                    self.bs[name] = str(self.bs[key])
                    self.sf[name] = str(self.sf[key])
            else:
                self.bs[name] = str(self.data_dict[name]['Background Shift'])
                self.sf[name] = str(self.data_dict[name]['Scaling Factor'])


    def _removeScanSelection(self):
        idx = self.selectedScans.currentIndex()
        name = self.selectedScans.currentText()
        # takes care of proper indexing
        if idx == self.previousIdx and self.previousIdx != 0:
            self.previousIdx = self.previousIdx - 1

        self.selectedScans.removeItem(idx)  # selected scans case where all scans have same bs and sf
        self.fit.pop(idx)
        self.bounds.pop(idx)
        self.weights.pop(idx)

        del self.bs[name]
        del self.sf[name]

        self.setTable()  # makes sure that the table is switched
        self.myPlotting()

    def addBoundWeight(self):
        col = self.boundWeightTable.columnCount()
        idx = self.selectedScans.currentIndex()
        n = len(self.bounds[idx])
        upper = self.bounds[idx][n-1][1]
        self.bounds[idx][n-1][1] = ''
        self.bounds[idx].append(['',upper])
        self.weights[idx].append('1')
        self.boundWeightTable.setColumnCount(col+1)

        self.setTable()

    def removeBoundWeight(self):
        col = self.boundWeightTable.columnCount()
        idx = self.selectedScans.currentIndex()  # gets the selected scan

        if col != 1:
            n = len(self.bounds[idx])  # get the number of boundaries
            upper = self.bounds[idx][n - 1][1]  # gets the proper upper boundary
            self.bounds[idx][n - 2][1] = upper
            self.bounds[idx].pop()
            self.weights[idx].pop()
            self.boundWeightTable.setColumnCount(col-1)


        self.setTable()

    def readTable(self):
        idx = self.selectedScans.currentIndex()
        name = self.selectedScans.currentText()

        row = self.boundWeightTable.rowCount()
        column = self.boundWeightTable.columnCount()
        if name != '':
            E = self.data_dict[name]['Energy']
            for i in range(row):
                for j in range(column):

                    item = self.boundWeightTable.item(i,j).text()
                    if i == 0:
                        if 'Angle' not in list(self.data_dict[name].keys()):
                            if not(self.axis_state):  # we are in angle state
                                if len(item) != 0:
                                    item = str(np.sin(float(item)*np.pi/180)*(E * 0.001013546143))
                        if len(item) != 0:

                            if len(item) < 8:
                                self.bounds[self.previousIdx][j][0] = item
                            else:
                                self.bounds[self.previousIdx][j][0] = item[:7]
                        else:
                            self.bounds[self.previousIdx][j][0] = ''
                    elif i == 1:
                        if 'Angle' not in list(self.data_dict[name].keys()):
                            if not (self.axis_state):  # we are in angle state
                                if len(item) != 0:
                                    item = str(np.sin(float(item)*np.pi/180)*(E * 0.001013546143))

                        if len(item) != 0:
                            if len(item) < 8:
                                self.bounds[self.previousIdx][j][1] = item
                            else:
                                self.bounds[self.previousIdx][j][1] = item[:7]
                        else:
                            self.bounds[self.previousIdx][j][1] = ''
                    elif i == 2:

                        self.weights[self.previousIdx][j] = item


            self.previousIdx = idx

    def plot_scans(self):


        if len(self.data) != 0:
            #self.sample = self.sWidget.sample

            self.spectrumWidget.clear()
            idx = self.whichScan.currentIndex()
            name = self.whichScan.currentText()
            if name != '':
                background_shift = self.data_dict[name]['Background Shift']
                scaling_factor = self.data_dict[name]['Scaling Factor']
                dat = self.data_dict[name]['Data']
                pol = self.data_dict[name]['Polarization']
                scan_type = self.data[idx][1]
                step_size = float(self.sWidget._step_size)

                if scan_type == 'Reflectivity':
                    qz = dat[0]

                    R = dat[2]
                    E = self.data_dict[name]['Energy']

                    qz, Rsim = self.sample.reflectivity(E,qz, s_min=step_size, bShift=background_shift, sFactor=scaling_factor)
                    Theta = np.arcsin(qz / (E * 0.001013546143))*180/np.pi
                    Rsim = Rsim[pol]
                    if pol == 'S' or pol =='P' or pol =='LC' or pol == 'RC':

                        if self.axis_state:
                            self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                            self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                        else:
                            self.spectrumWidget.plot(Theta, R,pen=pg.mkPen((0,2), width=2), name='Data')
                            self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                        self.spectrumWidget.setLabel('left', "Reflectivity, R")
                        self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                        self.spectrumWidget.setLogMode(False,True)
                    elif pol == 'AL' or pol =='AC':
                        rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                        if self.axis_state:

                            self.spectrumWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                            self.spectrumWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2), name='Simulation')
                        else:
                            self.spectrumWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                            self.spectrumWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2), name='Simulation')

                        self.spectrumWidget.setLogMode(False, False)
                        self.spectrumWidget.setLabel('left', "Reflectivity, R")
                        self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                elif scan_type == 'Energy':
                    E = dat[3]
                    R = dat[2]
                    Theta = self.data_dict[name]['Angle']
                    E, Rsim = self.sample.energy_scan(Theta,E, s_min=step_size, bShift=background_shift, sFactor=scaling_factor)
                    Rsim = Rsim[pol]
                    self.spectrumWidget.setLogMode(False, False)
                    self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                    self.spectrumWidget.plot(E, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Energy, E (eV)")

    def plot_selected_scans(self):
        step_size = float(self.sWidget._step_size)
        self.sample = self.sWidget.sample
        self.spectrumWidget.clear()
        name = self.selectedScans.currentText()
        b_idx = self.selectedScans.currentIndex()


        if name != '':
            bound = self.bounds[b_idx]
            lower = float(bound[0][0])
            upper = float(bound[-1][-1])
            background_shift = self.data_dict[name]['Background Shift']
            scaling_factor = self.data_dict[name]['Scaling Factor']
            idx=0
            notDone = True
            while notDone and idx==len(self.data)-1:
                temp_name = self.data[idx][2]
                if temp_name == name:
                    notDone=False
                else:
                    idx = idx + 1


            dat = self.data_dict[name]['Data']
            pol = self.data_dict[name]['Polarization']

            if 'Angle' not in list(self.data_dict[name].keys()):
                qz = dat[0]
                R = dat[2]
                E = self.data_dict[name]['Energy']
                qz, Rsim = self.sample.reflectivity(E,qz, s_min=step_size, bShift=background_shift, sFactor=scaling_factor)
                Theta = np.arcsin(qz/(E*0.001013546143))*180/np.pi

                Rsim = Rsim[pol]
                if pol == 'S' or pol =='P' or pol =='LC' or pol == 'RC':

                    if self.axis_state:
                        self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                    else:
                        self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                        lower = np.arcsin(lower/(E*0.001013546143))*180/np.pi
                        upper = np.arcsin(upper/(E*0.001013546143))*180/np.pi
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.spectrumWidget.setLogMode(False,True)

                    self.spectrumWidget.setXRange(lower,upper)
                elif pol == 'AL' or pol =='AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    if self.axis_state:
                        self.spectrumWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2), name='Simulation')
                    else:
                        self.spectrumWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2), name='Simulation')

                    self.spectrumWidget.setLogMode(False, False)
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.spectrumWidget.setXRange(lower, upper)
            else:
                E = dat[3]
                R = dat[2]
                Theta = self.data_dict[name]['Angle']
                E, Rsim = self.sample.energy_scan(Theta,E, s_min=step_size, bShift=background_shift, sFactor=scaling_factor)
                Rsim = Rsim[pol]
                self.spectrumWidget.setLogMode(False, False)
                self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                self.spectrumWidget.plot(E, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                self.spectrumWidget.setLabel('left', "Reflectivity, R")
                self.spectrumWidget.setLabel('bottom', "Energy, E (eV)")

class Worker(QObject):
    finished = pyqtSignal()
    progress = pyqtSignal(int)


    def __init__(self, function):
        super().__init__()
        self.function = function

    def run(self):
        x, fun = self.function._optimizer()
        self.function.x = x
        self.function.fun = fun
        self.finished.emit()




class callback():
    def __init__(self):
        self.Finish = False

    def stop_evolution(self, x, convergence):

        if stop:
            return True
        else:
            return False

    def stop_simplicial(self, x):

        if stop:
            return True
        else:
            return False
    def stop_annealing(self, x, f, contect):

        if stop:
            return True
        else:
            return False

class GlobalOptimizationWidget(QWidget):
    def __init__(self, sWidget, rWidget, rApp):
        super().__init__()

        self.sWidget = sWidget
        self.rWidget = rWidget
        self.rApp = rApp

        self.sample = copy.deepcopy(self.sWidget.sample)
        self.temp_sample = copy.deepcopy(self.sample)
        self.sampleBounds = []
        self.sfBounds = []
        self.otherBounds = []

        self.x = []
        self.fun = 0
        self.callback = callback()
        self.progressFinished = True
        self.objective = 'Chi-Square'
        self.shape_weight = 0

        # plotting layout ---------------------------------------------------------------------------------------------
        plotLayout = QHBoxLayout()

        self.goParameters = {'differential evolution': ['currenttobest1bin',2,15, 1e-6, 0,0.5,1, 0.7, True,'latinhypercube','immediate'],
                             'simplicial homology': ['None', 1, 'simplicial'],
                             'dual annealing': [150, 5230.0,2e-5,2.62,5.0,10000000.0,True]}

        self.paramChange = True
        self.parameters = []
        selectedScansLayout = QHBoxLayout()
        self.selectedScans = QComboBox()  # shows the scans selected for the data fitting process
        self.selectedScans.currentTextChanged.connect(self.plot_scan)
        selectedScansLabel = QLabel('Fit Scans:')
        selectedScansLayout.addWidget(selectedScansLabel)
        selectedScansLayout.addWidget(self.selectedScans)

        # Adding the plotting Widget
        self.plotWidget = pg.PlotWidget()
        self.plotWidget.setBackground('w')

        self.plotWidget.addLegend()


        # Global optimization parameters and fitting
        buttonLayout = QVBoxLayout()


        self.runButton = QPushButton('Run Optimization')
        self.runButton.pressed.connect(self._run_global_optimization)
        self.runButton.setStyleSheet('background: green')

        self.stopButton = QPushButton('Stop Optimization')
        self.stopButton.clicked.connect(self._stop_optimization)
        self.optButton = QPushButton('Update Sample')
        self.optButton.clicked.connect(self._save_optimization)
        self.optButton.setStyleSheet('background: cyan')

        buttonLayout.addWidget(self.selectedScans)
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(self.runButton)
        buttonLayout.addWidget(self.stopButton)
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(self.optButton)

        # Addinng Widgets to plotting layout
        plotLayout.addWidget(self.plotWidget)
        plotLayout.addLayout(buttonLayout)




        self.fittingParamTable = QTableWidget()
        self.fittingParamTable.setColumnCount(5)
        self.fittingParamTable.setHorizontalHeaderLabels(['Name', 'Current Value', 'Lower Boundary', 'Upper Boundary', 'New'])
        tableLayout = QVBoxLayout()
        tableLayout.addWidget(self.fittingParamTable)

        # Include the different input parameters for the global optimization and their algorithms
        self.goParamWidget = QWidget()

        # Adding objective function parameters

        # start with the layout of the "main" window
        algorithmLayout = QVBoxLayout()
        algorithmLabel = QLabel('Algorithm Selection')
        self.algorithmSelect = QComboBox()
        self.algorithmSelect.addItems(['differential evolution', 'simplicial homology', 'dual annealing'])
        self.algorithmSelect.currentIndexChanged.connect(self.change_algorithm)
        algorithmLayout.addWidget(algorithmLabel)
        algorithmLayout.addWidget(self.algorithmSelect)
        algorithmLayout.addStretch(1)

        # parameter Labels
        objectiveFunction = QLabel('Objective Function Parameters:')
        self.chi = QRadioButton('Chi-Square', self)
        self.chi.setChecked(True)
        self.chi.toggled.connect(self._changeObjectiveFunction)
        self.L1 = QRadioButton('L1-Norm', self)
        self.L1.toggled.connect(self._changeObjectiveFunction)
        self.L2 = QRadioButton('L2-Norm', self)
        self.L2.toggled.connect(self._changeObjectiveFunction)

        totLabel = QLabel('Total Variation: ')
        totLabel.setFixedWidth(80)
        totLayout = QHBoxLayout()
        self.totalVarWeight = QLineEdit('0')
        self.totalVarWeight.textChanged.connect(self._changeShapeWeight)
        self.totalVarWeight.setFixedWidth(50)
        totLayout.addWidget(totLabel)
        totLayout.addWidget(self.totalVarWeight)
        totLayout.addStretch(1)


        vbox = QVBoxLayout()
        vbox.addWidget(objectiveFunction)
        vbox.addWidget(self.chi)
        vbox.addWidget(self.L1)
        vbox.addWidget(self.L2)
        vbox.addLayout(totLayout)

        algorithmLayout.addLayout(vbox)
        algorithmWidget = QWidget()
        algorithmWidget.setStyleSheet("border: 1px solid black;")
        algorithmWidget.setLayout(algorithmLayout)

        self.goStackLayout = QStackedLayout()

        # differential evolution
        self.evolutionWidget = QWidget()
        evolutionLayout = QVBoxLayout()

        eStrategyLayout = QHBoxLayout()
        self.eStrategy = QComboBox()
        self.eStrategy.addItems(['best1bin','best1exp', 'rand1exp','randtobest1exp', 'best2exp', 'rand2exp',
                                  'randtobest1bin','currenttobest1bin', 'best2bin', 'rand2bin', 'rand1bin'])
        self.eStrategy.currentIndexChanged.connect(self.getGOParameters)
        stratLabel = QLabel('Strategy: ')
        stratLabel.setFixedWidth(70)
        eStrategyLayout.addWidget(stratLabel)
        eStrategyLayout.addWidget(self.eStrategy)
        evolutionLayout.addLayout(eStrategyLayout)

        eMaxiterLayout = QHBoxLayout()
        self.eMaxiter = QLineEdit()
        self.eMaxiter.textChanged.connect(self.getGOParameters)
        eMaxiterLabel = QLabel('maxIter')
        eMaxiterLabel.setFixedWidth(70)
        eMaxiterLayout.addWidget(eMaxiterLabel)
        eMaxiterLayout.addWidget(self.eMaxiter)
        evolutionLayout.addLayout(eMaxiterLayout)

        ePopsizeLayout = QHBoxLayout()
        self.ePopsize = QLineEdit()
        self.ePopsize.textChanged.connect(self.getGOParameters)
        popsizeLabel = QLabel('popsize: ')
        popsizeLabel.setFixedWidth(70)
        ePopsizeLayout.addWidget(popsizeLabel)
        ePopsizeLayout.addWidget(self.ePopsize)
        evolutionLayout.addLayout(ePopsizeLayout)


        eTolLayout = QHBoxLayout()
        self.eTol = QLineEdit()
        eTolLabel = QLabel('tol: ')
        self.eTol.textChanged.connect(self.getGOParameters)
        eTolLabel.setFixedWidth(70)
        eTolLayout.addWidget(eTolLabel)
        eTolLayout.addWidget(self.eTol)
        evolutionLayout.addLayout(eTolLayout)

        eAtolLayout = QHBoxLayout()
        self.eAtol = QLineEdit()
        self.eAtol.textChanged.connect(self.getGOParameters)
        eAtolLabel = QLabel('atol: ')
        eAtolLabel.setFixedWidth(70)
        eAtolLayout.addWidget(eAtolLabel)
        eAtolLayout.addWidget(self.eAtol)
        evolutionLayout.addLayout(eAtolLayout)

        eMinMutationLayout = QHBoxLayout()
        self.eMinMutation = QLineEdit()
        self.eMinMutation.textChanged.connect(self.getGOParameters)
        eMinMutationLabel = QLabel('min. mutation: ')
        eMinMutationLabel.setFixedWidth(70)
        eMinMutationLayout.addWidget(eMinMutationLabel)
        eMinMutationLayout.addWidget(self.eMinMutation)
        evolutionLayout.addLayout(eMinMutationLayout)

        eMaxMutationLayout = QHBoxLayout()
        self.eMaxMutation = QLineEdit()
        self.eMaxMutation.textChanged.connect(self.getGOParameters)
        eMaxMutationLabel = QLabel('max. mutation: ')
        eMaxMutationLabel.setFixedWidth(70)
        eMaxMutationLayout.addWidget(eMaxMutationLabel)
        eMaxMutationLayout.addWidget(self.eMaxMutation)
        evolutionLayout.addLayout(eMaxMutationLayout)

        eRecombLayout = QHBoxLayout()
        self.eRecomb = QLineEdit()
        self.eRecomb.textChanged.connect(self.getGOParameters)
        recombLabel = QLabel('recombination: ')
        recombLabel.setFixedWidth(70)
        eRecombLayout.addWidget(recombLabel)
        eRecombLayout.addWidget(self.eRecomb)
        evolutionLayout.addLayout(eRecombLayout)

        ePolishLayout = QHBoxLayout()
        self.ePolish = QCheckBox()
        self.ePolish.stateChanged.connect(self.getGOParameters)
        polishLabel = QLabel('polish')
        polishLabel.setFixedWidth(70)
        ePolishLayout.addWidget(polishLabel)
        ePolishLayout.addWidget(self.ePolish)
        evolutionLayout.addLayout(ePolishLayout)

        eInitLayout = QHBoxLayout()
        self.eInit = QComboBox()
        self.eInit.addItems(['latinhypercube','sobol','halton','random'])
        self.eInit.currentIndexChanged.connect(self.getGOParameters)
        initLabel = QLabel('init: ')
        initLabel.setFixedWidth(70)
        eInitLayout.addWidget(initLabel)
        eInitLayout.addWidget(self.eInit)
        evolutionLayout.addLayout(eInitLayout)

        eUpdatingLayout = QHBoxLayout()
        self.eUpdating = QComboBox()
        self.eUpdating.addItems(['immediate', 'deferred'])
        self.eUpdating.currentIndexChanged.connect(self.getGOParameters)
        updateLabel = QLabel('updating: ')
        updateLabel.setFixedWidth(70)
        eUpdatingLayout.addWidget(updateLabel)
        eUpdatingLayout.addWidget(self.eUpdating)
        evolutionLayout.addLayout(eUpdatingLayout)

        self.evolutionWidget.setLayout(evolutionLayout)

        # shgo algorithm
        shgoLayout = QVBoxLayout()
        self.shgoWidget = QWidget()

        shgoNLayout = QHBoxLayout()
        nLabel = QLabel('n: ')
        nLabel.setFixedWidth(70)
        self.shgoN = QLineEdit()
        self.shgoN.textChanged.connect(self.getGOParameters)
        shgoNLayout.addWidget(nLabel)
        shgoNLayout.addWidget(self.shgoN)
        shgoLayout.addLayout(shgoNLayout)

        shgoIterLayout = QHBoxLayout()
        iterLabel = QLabel('iter: ')
        iterLabel.setFixedWidth(70)
        self.shgoIter = QLineEdit()
        self.shgoIter.textChanged.connect(self.getGOParameters)
        shgoIterLayout.addWidget(iterLabel)
        shgoIterLayout.addWidget(self.shgoIter)
        shgoLayout.addLayout(shgoIterLayout)

        shgoSamplingLayout = QHBoxLayout()
        samplingLabel = QLabel('sampling: ')
        samplingLabel.setFixedWidth(70)
        self.shgoSampling = QComboBox()
        self.shgoSampling.addItems(['simplicial','halton','sobol'])
        self.shgoSampling.currentIndexChanged.connect(self.getGOParameters)
        shgoSamplingLayout.addWidget(samplingLabel)
        shgoSamplingLayout.addWidget(self.shgoSampling)
        shgoLayout.addLayout(shgoSamplingLayout)

        self.shgoWidget.setLayout(shgoLayout)

        # dual annealing parameters

        self.dualWidget = QWidget()
        dualLayout = QVBoxLayout()

        dualMaxiterLayout = QHBoxLayout()
        dualMaxiterLabel = QLabel('maxiter: ')
        dualMaxiterLabel.setFixedWidth(70)
        self.dualMaxiter = QLineEdit()
        self.dualMaxiter.textChanged.connect(self.getGOParameters)
        dualMaxiterLayout.addWidget(dualMaxiterLabel)
        dualMaxiterLayout.addWidget(self.dualMaxiter)
        dualLayout.addLayout(dualMaxiterLayout)

        dualInitTempLayout = QHBoxLayout()
        dualInitTempLabel = QLabel('initial temp: ')
        dualInitTempLabel.setFixedWidth(70)
        self.dualInitTemp = QLineEdit()
        self.dualInitTemp.textChanged.connect(self.getGOParameters)
        dualInitTempLayout.addWidget(dualInitTempLabel)
        dualInitTempLayout.addWidget(self.dualInitTemp)
        dualLayout.addLayout(dualInitTempLayout)

        dualRestartTempLayout = QHBoxLayout()
        dualRestartTempLabel = QLabel('restart temp: ')
        dualRestartTempLabel.setFixedWidth(70)
        self.dualRestartTemp = QLineEdit()
        self.dualRestartTemp.textChanged.connect(self.getGOParameters)
        dualRestartTempLayout.addWidget(dualRestartTempLabel)
        dualRestartTempLayout.addWidget(self.dualRestartTemp)
        dualLayout.addLayout(dualRestartTempLayout)

        dualVisitLayout = QHBoxLayout()
        dualVisitLabel = QLabel('visit: ')
        dualVisitLabel.setFixedWidth(70)
        self.dualVisit = QLineEdit()
        self.dualVisit.textChanged.connect(self.getGOParameters)
        dualVisitLayout.addWidget(dualVisitLabel)
        dualVisitLayout.addWidget(self.dualVisit)
        dualLayout.addLayout(dualVisitLayout)

        dualAcceptLayout = QHBoxLayout()
        dualAcceptLabel = QLabel('accept: ')
        dualAcceptLabel.setFixedWidth(70)
        self.dualAccept = QLineEdit()
        self.dualAccept.textChanged.connect(self.getGOParameters)
        dualAcceptLayout.addWidget(dualAcceptLabel)
        dualAcceptLayout.addWidget(self.dualAccept)
        dualLayout.addLayout(dualAcceptLayout)

        dualMaxfunLayout = QHBoxLayout()
        dualMaxfunLabel = QLabel('maxfun: ')
        dualMaxfunLabel.setFixedWidth(70)
        self.dualMaxfun = QLineEdit()
        self.dualMaxfun.textChanged.connect(self.getGOParameters)
        dualMaxfunLayout.addWidget(dualMaxfunLabel)
        dualMaxfunLayout.addWidget(self.dualMaxfun)
        dualLayout.addLayout(dualMaxfunLayout)

        dualLocalLayout = QHBoxLayout()
        dualLocalLabel = QLabel('local search: ')
        dualLocalLabel.setFixedWidth(70)
        self.dualLocal = QCheckBox()
        self.dualLocal.stateChanged.connect(self.getGOParameters)
        dualLocalLayout.addWidget(dualLocalLabel)
        dualLocalLayout.addWidget(self.dualLocal)
        dualLayout.addLayout(dualLocalLayout)


        self.dualWidget.setLayout(dualLayout)

        # adding the algorithm widgets to stacked layout
        self.goStackLayout.addWidget(self.evolutionWidget)
        self.goStackLayout.addWidget(self.shgoWidget)
        self.goStackLayout.addWidget(self.dualWidget)

        goLayout = QHBoxLayout()
        goLayout.addWidget(algorithmWidget)
        goLayout.addLayout(self.goStackLayout)

        self.goParamWidget.setLayout(goLayout)

        bottomLayout = QHBoxLayout()
        bottomLayout.addLayout(tableLayout)
        bottomLayout.addWidget(self.goParamWidget)



        pagelayout = QVBoxLayout()
        pagelayout.addLayout(plotLayout)
        pagelayout.addLayout(bottomLayout)
        self.setGOParameters()
        self.setLayout(pagelayout)
        self.setTableFit()

    def _changeObjectiveFunction(self):
        rbtn = self.sender()

        if rbtn.isChecked() == True:
            self.objective = rbtn.text()
    def _changeShapeWeight(self):
        value = self.totalVarWeight.text()
        if value != '':
            self.shape_weight = float(self.totalVarWeight.text())


    def _stop_optimization(self):
        global stop
        stop = True

    def changeFitParameters(self):
        # This function simply takes all the fitting parameters and saves the new fit
        for idx,fit in enumerate(self.parameters):
            if type(fit[0]) != str: # structural, polymorphous, magnetic
                layer = fit[0]
                my_type = fit[1]
                if my_type == 'STRUCTURAL':
                    mode = fit[2]
                    if mode == 'COMPOUND':
                        char = fit[3]
                        ele_idx = fit[4]  # keeps track of the element index
                        if char == 'THICKNESS':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][1])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][1] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][1] = float(self.sWidget.structTableInfo[layer][i][1]) - diff

                        elif char == 'DENSITY':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][2])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                s = float(self.sWidget.structTableInfo[layer][i][6])
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][2] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][2] = float(self.sWidget.structTableInfo[layer][i][2]) - s*diff

                        elif char == 'ROUGHNESS':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][3])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][3] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][3] = float(self.sWidget.structTableInfo[layer][i][3]) - diff

                        elif char == 'LINKED ROUGHNESS':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][4])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][4] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][4] = float(self.sWidget.structTableInfo[layer][i][4]) - diff

                    elif mode == 'ELEMENT':  # element mode
                        element = fit[3]
                        char = fit[4]

                        ele_idx = 0
                        for i in range(len(self.sWidget.structTableInfo[layer])):
                            if self.sWidget.structTableInfo[layer][i][0] == element:
                                ele_idx = i

                        if char == 'THICKNESS':
                            self.sWidget.structTableInfo[layer][ele_idx][1] = self.x[idx]
                        elif char == 'DENSITY':
                            self.sWidget.structTableInfo[layer][ele_idx][2] = self.x[idx]
                        elif char == 'ROUGHNESS':
                            self.sWidget.structTableInfo[layer][ele_idx][3] = self.x[idx]
                        elif char == 'LINKED ROUGHNESS':
                            self.sWidget.structTableInfo[layer][ele_idx][4] = self.x[idx]
                elif my_type == 'POLYMORPHOUS':
                    element = fit[2]
                    poly = fit[3]
                    j = list(self.sWidget.varData[element][layer][0]).index(poly)
                    self.sWidget.varData[element][layer][1][j] = self.x[idx]

                    # will need to change for more than 2 element variations
                    if j == 1:
                        self.sWidget.varData[element][layer][1][0] = 1-float(self.x[idx])
                    elif j == 0:
                        self.sWidget.varData[element][layer][1][1] = 1-float(self.x[idx])

                elif my_type == 'MAGNETIC':
                    if len(fit) == 3:
                        element = fit[2]
                        self.sWidget.magData[element][layer][1][0] = self.x[idx]
                    elif len(fit) == 4:
                        element = fit[2]
                        poly = fit[3]
                        j = list(self.sWidget.magData[element][layer][0]).index(poly)
                        self.sWidget.magData[element][layer][1][j] = self.x[idx]

            else: # scattering factor, background shift, scaling factor
                if fit[0] == 'SCATTERING FACTOR':
                    my_type = fit[1]
                    if my_type == 'STRUCTURAL':
                        sf = fit[2]
                        name = 'ff-' + sf
                        self.sWidget.eShift[name] = self.x[idx]
                    elif my_type == 'MAGNETIC':
                        sf = fit[2]
                        name = 'ffm-' + sf
                        self.sWidget.eShift[name] = self.x[idx]

                elif fit[0] == 'BACKGROUND SHIFT':
                    scans = self.rWidget.fit
                    if fit[1] == 'ALL SCANS':
                        for scan in scans:
                            self.rWidget.bs[scan] = "{:e}".format(self.x[idx])
                    else:
                        scan = fit[1]
                        self.rWidget.bs[scan] = "{:e}".format(self.x[idx])

                elif fit[0] == 'SCALING FACTOR':
                    scans = self.rWidget.fit
                    if fit[1] == 'ALL SCANS':
                        for scan in scans:
                            self.rWidget.sf[scan] = str(self.x[idx])
                    else:
                        scan = fit[1]
                        self.rWidget.sf[scan] = str(self.x[idx])


    def _save_optimization(self):


        row = 0
        # first we need to change the boundaries
        for idx in range(len(self.sWidget.currentVal)):
            fit = self.parameters[idx]
            self.sWidget.currentVal[idx][0] = str(self.x[row])
            if fit[1] == 'STRUCTURAL':
                if fit[2] == 'COMPOUND':
                    if fit[3] == "THICKNESS":
                        lower = self.x[row] - 5
                        upper = self.x[row] + 5
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                    elif fit[3] == "DENSITY":
                        lower = self.x[row] - 0.01
                        upper = self.x[row] + 0.01
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                    elif fit[3] == "ROUGHNESS":
                        lower = self.x[row] - 1
                        upper = self.x[row] + 1
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                    elif fit[3] == "LINKED ROUGHNESS":
                        lower = self.x[row] - 1
                        upper = self.x[row] + 1
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                elif fit[2] == 'ELEMENT':
                    if fit[4] == "THICKNESS":
                        lower = self.x[row] - 5
                        upper = self.x[row] + 5
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                    elif fit[4] == "DENSITY":
                        lower = self.x[row] - 0.01
                        upper = self.x[row] + 0.01
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                    elif fit[4] == "ROUGHNESS":
                        lower = self.x[row] - 1
                        upper = self.x[row] + 1
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
                    elif fit[4] == "LINKED ROUGHNESS":
                        lower = self.x[row] - 1
                        upper = self.x[row] + 1
                        if lower < 0:
                            lower = 0

                        self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]

            elif fit[1] == 'POLYMORPHOUS':
                lower = self.x[row] - 0.2
                upper = self.x[row] + 0.2
                if lower < 0:
                    lower = 0
                if upper > 1:
                    upper = 1
                self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]
            elif fit[1] == 'MAGNETIC':
                lower = self.x[row] - 0.01
                upper = self.x[row] + 0.01
                if lower < 0:
                    lower = 0
                self.sWidget.currentVal[idx][1] = [str(lower), str(upper)]

            row = row + 1

        for idx in range(len(self.rWidget.currentVal)):
            fit = self.parameters[idx]
            self.rWidget.currentVal[idx][0] = str(self.x[row])
            if fit[0] == "BACKGROUND SHIFT":
                lower = str(self.x[row] - 5e-8)
                upper = str(self.x[row] + 5e-8)
                self.rWidget.currentVal[idx][1] = [lower, upper]
            elif fit[0] == "SCALING FACTOR":
                lower = str(self.x[row] - 0.2)
                upper = str(self.x[row] + 0.2)
                self.rWidget.currentVal[idx][1] = [lower, upper]
            elif fit[0] == "SCATTERING FACTOR":
                lower = str(self.x[row] - 0.5)
                upper = str(self.x[row] + 0.5)
                self.rWidget.currentVal[idx][1] = [lower, upper]
            row = row + 1

        self.changeFitParameters()

        self.sWidget.sample, self.rWidget.bs, self.rWidget.sf = go.changeSampleParams(self.x, self.parameters, self.sWidget.sample, self.rWidget.bs, self.rWidget.sf)

        # update the background shift and scaling factors
        name = self.rWidget.selectedScans.currentText()
        self.rWidget.backgroundShift.blockSignals(True)
        self.rWidget.scalingFactor.blockSignals(True)
        self.rWidget.backgroundShift.setText(self.rWidget.bs[name])
        self.rWidget.scalingFactor.setText(self.rWidget.sf[name])
        self.rWidget.backgroundShift.blockSignals(False)
        self.rWidget.scalingFactor.blockSignals(False)

        self.setTableFit()
        self.sWidget.setTable()
        self.sWidget.setTableVar()
        self.sWidget.setTableMag()


    def plot_scan(self):

        self.plotWidget.clear()
        name = self.selectedScans.currentText()

        if name != '':
            dat = self.rWidget.data_dict[name]['Data']
            pol = self.rWidget.data_dict[name]['Polarization']


            idx = 0
            notDone = True
            while notDone and idx == len(self.rWidget.data) - 1:
                temp_name = self.rWidget.data[idx][2]
                if temp_name == name:
                    notDone = False
                else:
                    idx = idx + 1

            scan_type = self.rWidget.data[idx][1]
            step_size = float(self.sWidget._step_size)

            sample1 = copy.deepcopy(self.sample)
            sample2 = self.sample  # just to make python happy
            isGO = False

            backS = copy.deepcopy(self.rWidget.bs)
            scaleF = copy.deepcopy(self.rWidget.sf)
            if len(self.x) != 0:
                sample2, backS, scaleS = go.changeSampleParams(self.x,self.parameters, copy.deepcopy(sample1), backS, scaleF)
                isGO = True

            scaling_factor_old = float(self.rWidget.sf[name])
            background_shift_old = float(self.rWidget.bs[name])

            scaling_factor = float(scaleF[name])
            background_shift = float(backS[name])

            if scan_type == 'Reflectivity':
                qz = dat[0]


                R = dat[2]
                E = self.rWidget.data_dict[name]['Energy']
                qz, Rsim = sample1.reflectivity(E, qz, s_min=step_size,sFactor=scaling_factor_old, bShift=background_shift_old)
                if isGO:
                    qz,Rgo = sample2.reflectivity(E, qz, s_min=step_size, sFactor=scaling_factor, bShift=background_shift)
                Theta = np.arcsin(qz / (E * 0.001013546143)) * 180 / np.pi
                Rsim = Rsim[pol]

                if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(qz, Rsim, pen=pg.mkPen((1, 3), width=2), name='Simulation')
                        if isGO:

                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size,sFactor=scaling_factor, bShift=background_shift)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(qz, Rgo, pen=pg.mkPen((2, 3), width=2), name='Optimized')

                    else:

                        self.plotWidget.plot(Theta, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 3), width=2), name='Simulation')
                        if isGO:
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size,sFactor=scaling_factor, bShift=background_shift)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(Theta, Rgo, pen=pg.mkPen((2, 3), width=2), name='Optimized')

                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.plotWidget.setLogMode(False, True)
                elif pol == 'AL' or pol == 'AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 3), width=2), name='Simulation')
                        if isGO:
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size,bShift=background_shift, sFactor=scaling_factor)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(qz[rm_idx], Rgo[rm_idx], pen=pg.mkPen((2, 3), width=2), name='Optimized')
                    else:
                        self.plotWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 3), width=2), name='Simulation')
                        if isGO:
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size, sFactor=scaling_factor, bShift=background_shift)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(Theta[rm_idx], Rgo[rm_idx], pen=pg.mkPen((2, 3), width=2), name='Optimized')

                    self.plotWidget.setLogMode(False, False)
                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
            elif scan_type == 'Energy':
                E = dat[3]
                R = dat[2]
                Theta = self.rWidget.data_dict[name]['Angle']
                E, Rsim = sample1.energy_scan(Theta, E, s_min=step_size, sFactor=scaling_factor_old, bShift=background_shift_old)
                if isGO:
                    qz,Rgo = sample2.energy_scan(Theta, E, s_min=step_size, sFactor=scaling_factor, bShift=background_shift)
                    Rgo = Rgo[pol]

                Rsim = Rsim[pol]
                self.plotWidget.setLogMode(False, False)
                self.plotWidget.plot(E, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                self.plotWidget.plot(E, Rsim, pen=pg.mkPen((1, 3), width=2), name='Simulation')
                if isGO:
                    qz, Rgo = sample2.energy_scan(Theta, E, s_min=step_size, sFactor=scaling_factor, bShift=background_shift)
                    Rgo = Rgo[pol]
                    self.plotWidget.plot(E,Rgo, pen=pg.mkPen((2, 3), width=2), name='Optimized')
                self.plotWidget.setLabel('left', "Reflectivity, R")
                self.plotWidget.setLabel('bottom', "Energy, E (eV)")


    def run_first(self):
        global stop
        stop = False
        self.sample = copy.deepcopy(self.sWidget._createSample())
        self.temp_sample = copy.deepcopy(self.sample)  # makes sure that
        self.worker = Worker(self)
        self.runButton.setStyleSheet('background: red')
        self.runButton.blockSignals(True)

    def optimizationFinished(self):
        global stop
        stop = True
        self.sample = copy.deepcopy(self.temp_sample)

        # purpose of this is to reset the structure from anything the user did before optimization finished
        self.sWidget._setStructFromSample(self.sample)
        self.sWidget._setVarMagFromSample(self.sample)
        self.sWidget.getData()
        self.sWidget.setTable()
        self.sWidget.setTableMag()
        self.sWidget.setTableVar()

        # make sure that I all the other parameters are returned back to original value after the global optimization

        self.worker = Worker(self)
        self.plot_scan()
        self.setTableFit()
        self.runButton.setStyleSheet('background: green')
        self.runButton.blockSignals(False)

    def _run_global_optimization(self):

        # runs the optimizer method to perform the global optimization
        self.thread = QThread()  # initialize the thread
        self.worker = Worker(self)
        self.worker.moveToThread(self.thread)

        self.thread.started.connect(self.run_first)
        self.thread.started.connect(self.worker.run)
        self.worker.finished.connect(self.thread.quit)
        self.worker.finished.connect(self.optimizationFinished)
        self.worker.finished.connect(self.worker.deleteLater)
        self.thread.finished.connect(self.thread.deleteLater)

        self.thread.start()


    def _optimizer(self):
        # getting the scans and putting them in their proper format
        # putting the parameters and their boundaries in the proper format!
        parameters = copy.deepcopy(self.sWidget.parameterFit)
        for fit in self.rWidget.sfBsFitParams:
            parameters.append(fit)

        self.parameters = parameters  # needed for creating new sample
        lw = []
        up = []

        for b in self.sWidget.currentVal:
            lw.append(float(b[1][0]))
            up.append(float(b[1][1]))
        for b in self.rWidget.currentVal:
            lw.append(float(b[1][0]))
            up.append(float(b[1][1]))

        bounds = list(zip(lw, up))

        scans = copy.deepcopy(self.rWidget.fit)

        # determines the boundaries of the scans
        sBounds = []
        for bound in self.rWidget.bounds:
            temp = []
            for b in bound:
                temp.append((float(b[0]), float(b[1])))
            sBounds.append(temp)

        sWeights = []
        for weight in self.rWidget.weights:
            temp = []
            for w in weight:
                temp.append(float(w))
            sWeights.append(temp)

        x = []
        fun = 0
        data_dict = self.rWidget.data_dict
        data = self.rWidget.data
        sample = copy.deepcopy(self.sample)

        backS = copy.deepcopy(self.rWidget.bs)
        scaleF = copy.deepcopy(self.rWidget.sf)

        idx = self.algorithmSelect.currentIndex()

        if len(parameters) != 0 and len(scans) != 0:
            if idx == 0:
                x, fun = go.differential_evolution(sample, data,data_dict, scans, backS, scaleF, parameters, bounds, sBounds, sWeights,
                                                   self.goParameters['differential evolution'],self.callback,
                                                   self.objective, self.shape_weight)
            elif idx == 1:
                x, fun = go.shgo(sample,data,data_dict, scans, backS, scaleF, parameters, bounds, sBounds, sWeights,
                                 self.goParameters['simplicial homology'], self.callback,
                                 self.objective, self.shape_weight)
            elif idx == 2:
                x, fun = go.dual_annealing(sample, data, data_dict, scans, backS, scaleF, parameters, bounds, sBounds, sWeights,
                                           self.goParameters['dual annealing'], self.callback,
                                           self.objective, self.shape_weight)
        else:
            print('Try try again')
        print('Done')

        return x, fun

    def getGOParameters(self):


        idx = self.algorithmSelect.currentIndex()
        if idx == 0:
            self.goParameters['differential evolution'][0] = self.eStrategy.currentText()
            self.goParameters['differential evolution'][1] = self.eMaxiter.text()
            self.goParameters['differential evolution'][2] = self.ePopsize.text()
            self.goParameters['differential evolution'][3] = self.eTol.text()
            self.goParameters['differential evolution'][4] = self.eAtol.text()
            self.goParameters['differential evolution'][5] = self.eMinMutation.text()
            self.goParameters['differential evolution'][6] = self.eMaxMutation.text()
            self.goParameters['differential evolution'][7] = self.eRecomb.text()
            if self.ePolish.checkState == 0:
                self.goParameters['differential evolution'][8] = 'True'
            else:
                self.goParameters['differential evolution'][8] = 'False'
            self.goParameters['differential evolution'][9] = self.eInit.currentText()
            self.goParameters['differential evolution'][10] = self.eUpdating.currentText()

        elif idx == 1:  # simplicial homology
            self.goParameters['simplicial homology'][0] = self.shgoN.text()
            self.goParameters['simplicial homology'][1] = self.shgoIter.text()
            self.goParameters['simplicial homology'][2] = self.shgoSampling.currentText()
        elif idx == 2:  # dual annealing
            self.goParameters['dual annealing'][0] = self.dualMaxiter.text()
            self.goParameters['dual annealing'][1] = self.dualInitTemp.text()
            self.goParameters['dual annealing'][2] = self.dualRestartTemp.text()
            self.goParameters['dual annealing'][3] = self.dualVisit.text()
            self.goParameters['dual annealing'][4] = self.dualAccept.text()
            self.goParameters['dual annealing'][5] = self.dualMaxfun.text()
            if self.dualLocal.checkState() == 0:
                self.goParameters['dual annealing'][6] = 'False'
            else:
                self.goParameters['dual annealing'][6] = 'True'

        self.setGOParameters()

    def setGOParameters(self):
        self.eStrategy.blockSignals(True)
        self.eMaxiter.blockSignals(True)
        self.ePopsize.blockSignals(True)
        self.eTol.blockSignals(True)
        self.eAtol.blockSignals(True)
        self.eMinMutation.blockSignals(True)
        self.eMaxMutation.blockSignals(True)
        self.eRecomb.blockSignals(True)
        self.ePolish.blockSignals(True)
        self.eInit.blockSignals(True)
        self.eUpdating.blockSignals(True)
        self.shgoN.blockSignals(True)
        self.shgoIter.blockSignals(True)
        self.shgoSampling.blockSignals(True)
        self.dualMaxiter.blockSignals(True)
        self.dualInitTemp.blockSignals(True)
        self.dualRestartTemp.blockSignals(True)
        self.dualVisit.blockSignals(True)
        self.dualAccept.blockSignals(True)
        self.dualMaxfun.blockSignals(True)
        self.dualLocal.blockSignals(True)


        idx = self.algorithmSelect.currentIndex()

        self.eStrategy.setCurrentText(self.goParameters['differential evolution'][0])
        self.eMaxiter.setText(str(self.goParameters['differential evolution'][1]))
        self.ePopsize.setText(str(self.goParameters['differential evolution'][2]))
        self.eTol.setText(str(self.goParameters['differential evolution'][3]))
        self.eAtol.setText(str(self.goParameters['differential evolution'][4]))
        self.eMinMutation.setText(str(self.goParameters['differential evolution'][5]))
        self.eMaxMutation.setText(str(self.goParameters['differential evolution'][6]))
        self.eRecomb.setText(str(self.goParameters['differential evolution'][7]))
        if str(self.goParameters['differential evolution'][8]) == 'True':
            self.ePolish.setCheckState(2)
        elif str(self.goParameters['differential evolution'][8]) == 'False':
            self.ePolish.setCheckState(0)
        self.eInit.setCurrentText(str(self.goParameters['differential evolution'][9]))
        self.eUpdating.setCurrentText(str(self.goParameters['differential evolution'][10]))


        self.shgoN.setText(str(self.goParameters['simplicial homology'][0]))
        self.shgoIter.setText(str(self.goParameters['simplicial homology'][1]))
        self.shgoSampling.setCurrentText(str(self.goParameters['simplicial homology'][2]))

        self.dualMaxiter.setText(str(self.goParameters['dual annealing'][0]))
        self.dualInitTemp.setText(str(self.goParameters['dual annealing'][1]))
        self.dualRestartTemp.setText(str(self.goParameters['dual annealing'][2]))
        self.dualVisit.setText(str(self.goParameters['dual annealing'][3]))
        self.dualAccept.setText(str(self.goParameters['dual annealing'][4]))
        self.dualMaxfun.setText(str(self.goParameters['dual annealing'][5]))
        if self.goParameters['dual annealing'][6] == 'False':
            self.dualLocal.setChecked(2)
        elif self.goParameters['dual annealing'][6] == 'True':
            self.dualLocal.setChecked(0)

        self.eStrategy.blockSignals(False)
        self.eMaxiter.blockSignals(False)
        self.ePopsize.blockSignals(False)
        self.eTol.blockSignals(False)
        self.eAtol.blockSignals(False)
        self.eMinMutation.blockSignals(False)
        self.eMaxMutation.blockSignals(False)
        self.eRecomb.blockSignals(False)
        self.ePolish.blockSignals(False)
        self.eInit.blockSignals(False)
        self.eUpdating.blockSignals(False)
        self.shgoN.blockSignals(False)
        self.shgoIter.blockSignals(False)
        self.shgoSampling.blockSignals(False)
        self.dualMaxiter.blockSignals(False)
        self.dualInitTemp.blockSignals(False)
        self.dualRestartTemp.blockSignals(False)
        self.dualVisit.blockSignals(False)
        self.dualAccept.blockSignals(False)
        self.dualMaxfun.blockSignals(False)
        self.dualLocal.blockSignals(False)

    def change_algorithm(self):
        # set the proper algorithm widget
        idx = self.algorithmSelect.currentIndex()
        self.goStackLayout.setCurrentIndex(idx)

    def updateScreen(self):

        # shows only the selected scans for fitting
        self.selectedScans.clear()
        for text in self.rWidget.fit:
            self.selectedScans.addItem(text)

    def setTableFit(self):

        # initialize the boundaries
        rows = len(self.sWidget.parameterFit) + len(self.rWidget.sfBsFitParams)  # total number of rows required
        self.fittingParamTable.setRowCount(rows)

        # create the names and set number of rows
        row = 0
        for idx,param in enumerate(self.sWidget.parameterFit):
            name = self.getName(param)
            value = str(self.sWidget.currentVal[idx][0])
            lower = str(self.sWidget.currentVal[idx][1][0])
            upper = str(self.sWidget.currentVal[idx][1][1])
            item1 = QTableWidgetItem(name)
            item2 = QTableWidgetItem(value)

            item3 = QTableWidgetItem(lower)
            item4 = QTableWidgetItem(upper)

            self.fittingParamTable.setItem(row,0, item1)
            self.fittingParamTable.setItem(row, 1, item2)
            self.fittingParamTable.setItem(row, 2, item3)
            self.fittingParamTable.setItem(row, 3, item4)

            row = row + 1

        for idx, param in enumerate(self.rWidget.sfBsFitParams):
            name = self.getName(param)
            value = str(self.rWidget.currentVal[idx][0])
            lower = str(self.rWidget.currentVal[idx][1][0])
            upper = str(self.rWidget.currentVal[idx][1][1])

            item1 = QTableWidgetItem(name)
            item2 = QTableWidgetItem(str(value))

            item3 = QTableWidgetItem(lower)
            item4 = QTableWidgetItem(upper)

            self.fittingParamTable.setItem(row, 0, item1)
            self.fittingParamTable.setItem(row, 1, item2)
            self.fittingParamTable.setItem(row, 2, item3)
            self.fittingParamTable.setItem(row, 3, item4)

            row = row + 1

        for idx in range(len(self.x)):
            item = QTableWidgetItem(str(self.x[idx]))
            self.fittingParamTable.setItem(idx, 4, item)


    def getName(self, p):
        name = ''
        n = len(p)
        shift = 0
        if n != 0:
            if type(p[0]) == int:  # sample parameters
                layer = p[0]
                param_type = p[1]

                if param_type == 'STRUCTURAL':  # structural case
                    mode = p[2]
                    if mode == 'ELEMENT':
                        element = p[3]
                        char = p[4]
                        if char == 'THICKNESS':
                            name = element + '-' + 'th. ' + str(layer)
                        elif char == 'DENSITY':
                            name = element + '-' + 'dens. ' + str(layer)
                        elif char == 'ROUGHNESS':
                            name = element + '-' + 'rough. ' + str(layer)
                        elif char == 'LINKED ROUGHNESS':
                            name = element + '-' + 'Lrough. ' + str(layer)
                    elif mode == 'COMPOUND':
                        char = p[3]
                        compound = ''

                        # gets all the elements in the layer
                        for e in self.sWidget.structTableInfo[layer]:
                            compound = compound + e[0]

                        if char == 'THICKNESS':
                            name = compound + '-th. ' + str(layer)
                        elif char == 'DENSITY':
                            name = compound + '-dens. ' + str(layer)
                        elif char == 'ROUGHNESS':
                            name = compound + '-rough. ' + str(layer)
                        elif char == 'LINKED ROUGHNESS':
                            name = compound + '-Lrough. ' + str(layer)

                elif param_type == 'POLYMORPHOUS':
                    var = p[-1]
                    name = var + ' -ratio ' + str(layer)

                elif param_type == 'MAGNETIC':
                    var = p[-1]
                    name = var + ' -mdens. ' + str(layer)

            elif p[0] == 'SCATTERING FACTOR':  # scattering factor case
                name = name + 'ff'
                scattering_factor = p[2]
                param_type = p[1]

                if param_type == 'MAGNETIC':
                    name = name + 'm-' + scattering_factor
                else:
                    name = name + '-' + scattering_factor

            elif p[0] == 'BACKGROUND SHIFT':  # background shift case
                if p[1] == 'ALL SCANS':
                    name = 'bShift-All'
                else:
                    name = 'bShift-'+p[1]
            elif p[0] == 'SCALING FACTOR':  # scaling factor
                if p[1] == 'ALL SCANS':
                    name = 'sFactor-All'
                else:
                    name = 'sFactor-' + p[1]

        return name

class ReflectometryApp(QMainWindow):
    def __init__(self):
        super().__init__()
        cwd = os.getcwd()
        self.fname = cwd + '\demo.h5'
        self.data = []
        self.data_dict = dict()
        self.sim_dict = dict()

        self.sample = ms.slab(1)  # app is initialized and no project is selected
        self.sample.addlayer(0,'SrTiO3',50)

        # set the title
        self.setWindowTitle('Reflectometry of Quantum Materials')

        # set the geometry of the window
        self.setGeometry(180,60,1400,800)

        pagelayout = QVBoxLayout()
        buttonlayout = QHBoxLayout()
        self.stackedlayout = QStackedLayout()

        pagelayout.addLayout(buttonlayout)
        pagelayout.addLayout(self.stackedlayout)

        label1 = QLabel('Label 1')
        label2 = QLabel('Label 2')
        label3 = QLabel('Label 3')

        self.sampleButton = QPushButton('Sample')
        self.reflButton = QPushButton('Reflectivity')
        self.goButton = QPushButton('Global Optimization')

        self._sampleWidget = sampleWidget(self.sample)  # initialize the sample widget
        self._reflectivityWidget = reflectivityWidget(self._sampleWidget, self.data, self.data_dict, self.sim_dict)
        self._goWidget = GlobalOptimizationWidget(self._sampleWidget, self._reflectivityWidget, self)

        self.sampleButton.setStyleSheet("background-color : magenta")
        self.sampleButton.clicked.connect(self.activate_tab_1)
        buttonlayout.addWidget(self.sampleButton)
        self.stackedlayout.addWidget(self._sampleWidget)

        self.reflButton.setStyleSheet("background-color : pink")
        self.reflButton.clicked.connect(self.activate_tab_2)
        buttonlayout.addWidget(self.reflButton)
        self.stackedlayout.addWidget(self._reflectivityWidget)

        self.goButton.setStyleSheet("background-color : pink")
        self.goButton.clicked.connect(self.activate_tab_3)
        buttonlayout.addWidget(self.goButton)
        self.stackedlayout.addWidget(self._goWidget)

        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)

        # starting my menu bar
        menuBar = self.menuBar()

        # saving and loading area
        fileMenu = QMenu("&File", self)

        self.newFile = QAction("&New", self)
        self.newFile.triggered.connect(self._newFile)
        fileMenu.addAction(self.newFile)

        self.loadFile = QAction("&Load", self)
        self.loadFile.triggered.connect(self._loadFile)
        fileMenu.addAction(self.loadFile)

        self.saveFile = QAction("&Save", self)
        self.saveFile.triggered.connect(self._saveFile)
        fileMenu.addAction(self.saveFile)

        self.saveAsFile = QAction("&Save As", self)
        self.saveAsFile.triggered.connect(self._saveAsFile)
        fileMenu.addAction(self.saveAsFile)
        fileMenu.addSeparator()
        self.saveSampleFile = QAction("&Save Sample", self)
        self.saveSampleFile.triggered.connect(self._saveSample)
        fileMenu.addAction(self.saveSampleFile)

        self.saveSimulationFile = QAction("&Save Simulation", self)
        self.saveSimulationFile.triggered.connect(self._saveSimulation)
        fileMenu.addAction(self.saveSimulationFile)
        fileMenu.addSeparator()
        self.importDataset = QAction("&Import Dataset", self)
        self.importDataset.triggered.connect(self._importDataSet)
        fileMenu.addAction(self.importDataset)

        self.saveSummary = QAction("&Save Summary", self)
        self.saveSummary.triggered.connect(self._summary)
        fileMenu.addAction(self.saveSummary)
        fileMenu.addSeparator()
        self.exitFile = QAction("&Exit", self)
        self.exitFile.triggered.connect(self._exitApplication)
        fileMenu.addAction(self.exitFile)

        # Tools menu
        menuBar.addMenu(fileMenu)
        toolsMenu = menuBar.addMenu("&Tools")
        self.script = QAction('Script', self)
        self.script.triggered.connect(self._script)
        toolsMenu.addAction(self.script)

        self.showFormFactor = QAction('Form Factors', self)
        self.showFormFactor.triggered.connect(self._showFormFactor)
        toolsMenu.addAction(self.showFormFactor)

        helpMenu = menuBar.addMenu("&Help")
        self.license = QAction('License', self)
        self.license.triggered.connect(self._license)
        helpMenu.addAction(self.license)

        self.about = QAction('About',self)
        self.about.triggered.connect(self._help)
        helpMenu.addAction(self.about)

    def _newFile(self):

        # create a new file with the inputted
        filename, _ = QFileDialog.getSaveFileName()
        fname = filename.split('/')[-1]

        # checks to make sure filename is in the correct format
        cont = True
        if filename == '' or fname == '':
            cont = False
        elif fname.endswith('.h5'):
            self.fname = filename  # change the file name that we will be using
        elif '.' not in fname:
            self.fname = filename + '.h5'
        else:
            cont = False


        if cont:  # create the new file
            ds.newFileHDF5(self.fname)

            sample = ms.slab(2)
            sample.addlayer(0, 'SrTiO3', 50)
            sample.addlayer(1, 'LaMnO3', 10)
            self.sample = sample
            self._sampleWidget.sample = sample

            self.data = list()
            self.data_dict = dict()
            self.sim_dict = dict()

            # loading in the background shift and scaling factor
            self._reflectivityWidget.bs = dict()
            self._reflectivityWidget.sf = dict()

            self._sampleWidget._setStructFromSample(sample)
            self._sampleWidget._setVarMagFromSample(sample)

            self._sampleWidget.eShift = dict()

            # now it's time to load the other information
            self._reflectivityWidget.selectedScans.clear()
            self._reflectivityWidget.whichScan.clear()

            self._reflectivityWidget.data = self.data
            self._reflectivityWidget.data_dict = self.data_dict

            for scan in self.data:
                self._reflectivityWidget.whichScan.addItem(scan[2])

            for key in list(sample.eShift.keys()):
                name = 'ff-' + key
                self._sampleWidget.eShift[name] = sample.eShift[key]

            for key in list(sample.mag_eShift.keys()):
                name = 'ffm-' + key
                self._sampleWidget.eShift[name] = sample.mag_eShift[key]

            self._sampleWidget.setTableEShift()

            layerList = []
            for i in range(len(sample.structure)):
                if i == 0:
                    layerList.append('Substrate')
                else:
                    layerList.append('Layer ' + str(i))

            self._sampleWidget.layerBox.clear()
            # change this for an arbitrary sample model
            self._sampleWidget.layerBox.addItems(layerList)

            self._goWidget.goParameters = {'differential evolution': ['currenttobest1bin',2,15, 1e-6, 0,0.5,1, 0.7, True,'latinhypercube','immediate'],
                             'simplicial homology': ['None', 1, 'simplicial'],
                             'dual annealing': [150, 5230.0,2e-5,2.62,5.0,10000000.0,True]}

            self._goWidget.setGOParameters()
            self._goWidget.setTableFit()


            # for now let's clear all the fitting parameters
            self._reflectivityWidget.sfBsFitParams = list()
            self._reflectivityWidget.currentVal = list()
            self._reflectivityWidget.rom = [True, False, False]
            self._reflectivityWidget.fit = list()
            self._reflectivityWidget.bounds = list()
            self._reflectivityWidget.weights = list()

            self._sampleWidget.parameterFit = list()
            self._sampleWidget.currentVal = list()

            self._goWidget.x = list()
            self._goWidget.fun = 0

            self._goWidget.parameters = []
            for param in self._sampleWidget.parameterFit:
                self._goWidget.parameters.append(param)
            for param in self._reflectivityWidget.sfbsFitParams:
                self._goWidget.parameters.append(param)

            # reset all of the tables!!!
            # - otherwise we have everything for the load function
        else:
            print('Unable to create the file')

    def _loadFile(self):

        self.fname, _ = QFileDialog.getOpenFileName(self, 'Open File')
        fname = self.fname.split('/')[-1]

        # when loading files I need to be able to scan the entire
        if fname.endswith('.h5') or fname.endswith('.all'):
            if fname.endswith('.h5') and fname != 'demo.h5':
                self.sample = ds.ReadSampleHDF5(self.fname)
                self._sampleWidget.sample = self.sample

                self.data, self.data_dict, self.sim_dict = ds.ReadDataHDF5(self.fname)

                # loading in the background shift and scaling factor
                self._reflectivityWidget.bs = dict()
                self._reflectivityWidget.sf = dict()


                self._sampleWidget._setStructFromSample(self.sample)
                self._sampleWidget._setVarMagFromSample(self.sample)


                # now it's time to load the other information
                self._reflectivityWidget.selectedScans.clear()
                self._reflectivityWidget.whichScan.clear()

                self._reflectivityWidget.data = self.data
                self._reflectivityWidget.data_dict = self.data_dict

                # make sure we are not plotting when we do not need to
                self._reflectivityWidget.whichScan.blockSignals(True)
                for scan in self.data:
                    self._reflectivityWidget.whichScan.addItem(scan[2])
                self._reflectivityWidget.whichScan.blockSignals(False)

                for key in list(self.sample.eShift.keys()):
                    name = 'ff-'+key
                    self._sampleWidget.eShift[name] = self.sample.eShift[key]

                for key in list(self.sample.mag_eShift.keys()):
                    name = 'ffm-'+key

                    self._sampleWidget.eShift[name] = self.sample.mag_eShift[key]

                self._sampleWidget.setTableEShift()

                layerList = []
                for i in range(len(self.sample.structure)):
                    if i == 0:
                        layerList.append('Substrate')
                    else:
                        layerList.append('Layer ' + str(i))

                self._sampleWidget.layerBox.clear()
                # change this for an arbitrary sample model
                self._sampleWidget.layerBox.addItems(layerList)

                self._goWidget.goParameters = ds.ReadAlgorithmHDF5(self.fname)
                self._goWidget.setGOParameters()
                self._goWidget.setTableFit()

                fitParams = ds.ReadFitHDF5(self.fname)

                # for now let's clear all the fitting parameters
                self._reflectivityWidget.sfBsFitParams = fitParams[0]
                self._reflectivityWidget.currentVal = fitParams[1]
                self._reflectivityWidget.rom = [True, False, False]
                self._reflectivityWidget.fit = fitParams[4]
                self._reflectivityWidget.bounds = fitParams[5]
                self._reflectivityWidget.weights = fitParams[6]

                self._sampleWidget.parameterFit = fitParams[2]
                self._sampleWidget.currentVal = fitParams[3]

                self._goWidget.x = fitParams[7]
                self._goWidget.fun = fitParams[8]

                self._goWidget.parameters = []
                for param in self._sampleWidget.parameterFit:
                    self._goWidget.parameters.append(param)
                for param in self._reflectivityWidget.sfBsFitParams:
                    self._goWidget.parameters.append(param)

                # reset all of the tables!!!
                # - otherwise we have everything for the load function
            elif fname.endswith('.all'):
                print('Currently no implementation')

    def _saveFile(self):
        # work on saving the current file
        # saving function is used to save entire workspace
        filename = self.fname
        fname = filename.split('/')[-1]

        if fname != 'demo.h5':
            self.sample = self._sampleWidget._createSample()
            self._sampleWidget.sample = self.sample
            self._reflectivityWidget.sample = self.sample

            # save the sample information to the file
            #ds.WriteSampleHDF5(self.fname, self.sample)

            data_dict = self.data_dict

            fitParams = [self._reflectivityWidget.sfBsFitParams,self._reflectivityWidget.currentVal,
                         self._sampleWidget.parameterFit, self._sampleWidget.currentVal,
                         self._reflectivityWidget.fit, self._reflectivityWidget.bounds,
                         self._reflectivityWidget.weights, self._goWidget.x, self._goWidget.fun]

            optParams = self._goWidget.goParameters

            ds.saveFileHDF5(fname, self.sample, data_dict,  fitParams, optParams)

    def _saveAsFile(self):
        # create a new file with the inputted
        filename, _ = QFileDialog.getSaveFileName()
        fname = filename.split('/')[-1]

        # checks to make sure filename is in the correct format
        cont = True
        if filename == '' or fname == '':
            cont = False
        elif fname.endswith('.h5'):
            self.fname = filename  # change the file name that we will be using
        elif '.' not in fname:
            self.fname = filename + '.h5'
        else:
            cont = False

        if cont and fname != 'demo.h5':  # create the new file
            data_dict = self.data_dict
            sim_dict = self.sim_dict
            fitParams = [self._reflectivityWidget.sfBsFitParams, self._reflectivityWidget.currentVal,
                         self._sampleWidget.parameterFit, self._sampleWidget.currentVal,
                         self._reflectivityWidget.fit, self._reflectivityWidget.bounds,
                         self._reflectivityWidget.weights, self._goWidget.x, self._goWidget.fun]

            optParams = self._goWidget.goParameters

            self.sample = self._sampleWidget._createSample()
            self._sampleWidget.sample = self.sample
            self._reflectivityWidget.sample = self.sample

            ds.saveAsFileHDF5(fname, self.sample,data_dict, sim_dict, fitParams, optParams)

    def _saveSimulation(self):
        sim_dict = copy.deepcopy(self.data_dict)

        fname = self.fname

        if len(sim_dict) != 0:
            loadingApp = LoadingScreen(self.sample, sim_dict)
            loadingApp.show()
            loadingApp.exec_()
            sim_dict = loadingApp.sim_dict
            loadingApp.close()

            ds.saveSimulationHDF5(fname, sim_dict)




    def _saveSample(self):
        self.sample = self._sampleWidget._createSample()
        self._sampleWidget.sample = self.sample
        self._reflectivityWidget.sample = self.sample

        # save the sample information to the file
        ds.WriteSampleHDF5(self.fname, self.sample)

    def _importDataSet(self):
        print('import dataset')

    def _summary(self):
        sample = self._sampleWidget._createSample()
        # create a new file with the inputted (save as a textfile!)
        filename, _ = QFileDialog.getSaveFileName()
        fname = filename.split('/')[-1]
        cont = True
        if fname.endswith('.txt'):
            self.fname = filename  # change the file name that we will be using
        elif fname == '':
            cont = False
        elif '.' not in fname:
            filename = filename + '.txt'

        if cont:
            with open(filename, 'w') as file:
                file.write("# Structure \n")  # header defining that the sample information is starting
                n = len(sample.structure)

                # General information for the sample model
                file.write("numberlayers = %s \n" % str(n))
                file.write("polyelements = %s \n" % str(sample.poly_elements))
                file.write("magelements = %s \n" % str(sample.mag_elements))
                file.write("layermagnetized = %s \n" % sample.layer_magnetized)
                file.write("energyShift = %s \n" % sample.eShift)
                file.write("magEnergyShift = %s \n\n" % sample.mag_eShift)

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

                file.write("# Fit Parameters \n\n")

                parameters = []
                for param in self._sampleWidget.parameterFit:
                    parameters.append(param)
                for param in self._reflectivityWidget.sfBsFitParams:
                    parameters.append(param)

                current_val = []
                for val in self._reflectivityWidget.currentVal:
                    current_val.append(val)
                for val in self._sampleWidget.currentVal:
                    current_val.append(val)

                file.write("scans = %s \n" % self._reflectivityWidget.fit)
                file.write("scanBounds = %s \n" % self._reflectivityWidget.bounds)
                file.write("scanWeights = %s \n\n" % self._reflectivityWidget.weights)

                file.write("fitParams = %s \n" % parameters)
                file.write("values = %s \n" % current_val)
                file.write("results = %s \n" % self._goWidget.x)
                file.write("chiSquare = %s \n\n" % self._goWidget.fun)

                file.write('# Optimization Algorithms \n\n')

                globOpt = self._goWidget.goParameters
                file.write("algorithm = differential_evolution \n")
                file.write("strategy = %s \n" % globOpt['differential evolution'][0])
                file.write("maxIter = %s \n" % globOpt['differential evolution'][1])
                file.write("popsize = %s \n" % globOpt['differential evolution'][2])
                file.write("tol = %s \n" % globOpt['differential evolution'][3])
                file.write("atol = %s \n" % globOpt['differential evolution'][4])
                file.write("minMutation = %s \n" % globOpt['differential evolution'][5])
                file.write("maxMutation = %s \n" % globOpt['differential evolution'][6])
                file.write("recombination = %s \n" % globOpt['differential evolution'][7])
                file.write("polish = %s \n" % globOpt['differential evolution'][8])
                file.write("init = %s \n" % globOpt['differential evolution'][9])
                file.write("updating = %s \n\n" % globOpt['differential evolution'][10])

                file.write("algorithm = simplicial_homology \n")
                file.write("n = %s \n" % globOpt['simplicial homology'][0])
                file.write("iter = %s \n" % globOpt['simplicial homology'][1])
                file.write("sampling = %s \n\n" % globOpt['simplicial homology'][2])

                file.write("algorithm = dual annealing \n")
                file.write("maxiter = %s \n" % globOpt['dual annealing'][0])
                file.write("initialTemp = %s \n" % globOpt['dual annealing'][1])
                file.write("restartTemp = %s \n" % globOpt['dual annealing'][2])
                file.write("visit = %s \n" % globOpt['dual annealing'][3])
                file.write("accept = %s \n" % globOpt['dual annealing'][4])
                file.write("maxfun = %s \n" % globOpt['dual annealing'][5])
                file.write("localSearch = %s \n\n" % globOpt['dual annealing'][6])
                file.close()

    def _exitApplication(self):
        # exit the program
        sys.exit()

    def _script(self):
        pass

    def _showFormFactor(self):
        pass

    def _license(self):
        pass

    def _help(self):
        pass
    def activate_tab_1(self):
        self.sample = copy.deepcopy(self._sampleWidget._createSample())
        self._reflectivityWidget.sample = self.sample
        self._goWidget.sample = self.sample

        self._sampleWidget.step_size.setText(self._sampleWidget._step_size)
        self.sampleButton.setStyleSheet("background-color : magenta")
        self.reflButton.setStyleSheet("background-color : pink")
        self.goButton.setStyleSheet("background-color : pink")
        self.stackedlayout.setCurrentIndex(0)
    def activate_tab_2(self):
        self.sample = copy.deepcopy(self._sampleWidget._createSample())
        self._reflectivityWidget.sample = self.sample
        self._goWidget.sample = self.sample

        self._reflectivityWidget.sample = self.sample
        self._reflectivityWidget.myPlotting()
        self._reflectivityWidget.stepWidget.setText(self._sampleWidget._step_size)
        self.sampleButton.setStyleSheet("background-color : pink")
        self.reflButton.setStyleSheet("background-color : magenta")
        self.goButton.setStyleSheet("background-color : pink")
        self.stackedlayout.setCurrentIndex(1)

    def activate_tab_3(self):
        self.sample = copy.deepcopy(self._sampleWidget._createSample())
        self._reflectivityWidget.sample = self.sample
        self._goWidget.sample = self.sample

        self.sampleButton.setStyleSheet("background-color : pink")
        self.reflButton.setStyleSheet("background-color : pink")
        self.goButton.setStyleSheet("background-color : magenta")
        self._goWidget.setTableFit()
        self._goWidget.updateScreen()

        self.stackedlayout.setCurrentIndex(2)


class LoadingScreen(QDialog):
    def __init__(self, sample, sim_dict):
        super().__init__()
        self.setWindowTitle('Saving Simulation')
        self.sample = sample
        self.sim_dict = sim_dict
        self.n = len(list(self.sim_dict.keys()))
        layout = QHBoxLayout()
        button = QPushButton('Start Save')
        button.clicked.connect(self.run)
        self.progress = QProgressBar(self)
        self.progress.setRange(0,self.n-1)
        self.progress.setVisible(True)
        layout.addWidget(button)
        layout.addWidget(self.progress)
        self.setLayout(layout)

    def run(self):

        my_keys = list(self.sim_dict.keys())
        n = self.n

        for idx in range(len(my_keys)):

            if idx % 2 == 0:
                self.progress.setValue(idx)

            key = my_keys[idx]
            pol = self.sim_dict[key]['Polarization']

            if 'Angle' in self.sim_dict[key].keys():  # energy scan
                E = self.sim_dict[key]['Data'][3]  # get energy axis
                Theta = self.sim_dict[key]['Angle']  # get angle
                E, R = self.sample.energy_scan(Theta, E)  # calculate energy scan
                R = R[pol]  # polarization
                self.sim_dict[key]['Data'][2] = list(R)
            else:  # reflectivity scan
                qz = self.sim_dict[key]['Data'][0]
                energy = self.sim_dict[key]['Energy']
                qz, R = self.sample.reflectivity(energy, qz)
                R = R[pol]
                self.sim_dict[key]['Data'][2] = list(R)


        self.accept()










if __name__ == '__main__':


    fname = 'Pim10uc.h5'

    app = QApplication(sys.argv)
    demo = ReflectometryApp()

    demo.show()

    sys.exit(app.exec_())