from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.figure as Figure
import sys
import material_structure as ms
import os
import pyqtgraph as pg
import data_structure as ds
import copy

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

        # gets the roughness

        # gets the linked roughness





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

        self.varTable = QTableWidget()
        self.elementBox = QComboBox()
        self.variationElements = sample.poly_elements

        self.magData = {ele: [[[''], [''], ['']] for i in range(len(sample.structure))] for ele in
                        sample.myelements}
        self.magDirection = ['z' for i in range(len(sample.structure))]
        self.magDirBox = QComboBox()
        self.getData() # gets the element variation and magnetic information

        self.magTable = QTableWidget()


        self.change_elements = False
        self.element_index = 0

        self.previousLayer = 0
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

        # changes the table on the screen when new layer selected
        self.layerBox.currentIndexChanged.connect(self.changetable)

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
        for i in range(len(self.sample.structure)):
            self.setStructFromSample(i)

        # setTable
        self.setTable(0)


        # Element Variation Stuff
        self.elementVariation = variationWidget(self,sample)
        self.elementMagnetic = magneticWidget(self,sample)

        self.sampleInfoLayout.addWidget(self.structTable)
        self.sampleInfoLayout.addWidget(self.elementVariation)
        self.sampleInfoLayout.addWidget(self.elementMagnetic)

        self.structTable.viewport().installEventFilter(self)
        self.varTable.viewport().installEventFilter(self)
        self.magTable.viewport().installEventFilter(self)

        selectlayout = QVBoxLayout()

        # buttons for choosing which parameters to choose
        structButton = QPushButton('Structure')
        structButton.clicked.connect(self._structural)
        structButton.clicked.connect(self.setTableVar)
        structButton.clicked.connect(self.setTableMag)
        selectlayout.addWidget(structButton)

        polyButton = QPushButton('Element Variation')
        polyButton.clicked.connect(self._elementVariation)
        polyButton.clicked.connect(self.setTableVar)
        selectlayout.addWidget(polyButton)

        magButton = QPushButton('Magnetic')
        magButton.clicked.connect(self._magnetic)
        magButton.clicked.connect(self.setTableMag)
        selectlayout.addWidget(magButton)

        self.structBool = True
        self.polyBool = False
        self.magBool = False

        dpButton = QPushButton('Density Profile')
        dpButton.clicked.connect(self.changetable)
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


        self.setLayout(mylayout)

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

        return False

    def mag_handler(self):
        idx = self.sampleInfoLayout.currentIndex()  # keeps track of which parameters are to be fit
        my_layer = self.layerBox.currentIndex()  # layer index

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

        # I now have the layer, element and the name of the element
        if column == 0:
            pass
        elif column == 1:
            scattering_factor = self.magTable.item(row, column).text()
            self.parameterFit.append(['Magnetic Scattering Factor'])

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
                    if fit[1] == 'Polymorphous' and fit[0] == my_layer:
                        name = self.varTable.item(row, 0).text()
                        if element == fit[2] and name == fit[3]:
                            alreadySelected = True
                elif column == 2:
                    scattering_factor = self.varTable.item(row,column).text()
                    if len(fit) == 2 and scattering_factor == fit[1]:
                        alreadySelected = True

            if not alreadySelected:
                if column == 1:
                    name = self.varTable.item(row, 0).text()
                    ratio = self.varTable.item(row, 1).text()

                    if ratio != '':
                        self.parameterFit.append([my_layer, 'Polymorphous', element, name])
                elif column == 2:
                    scattering_factor = self.varTable.item(row,2).text()
                    if scattering_factor != '':
                        self.parameterFit.append(['Scattering Factor', scattering_factor])

        elif action == _remove_fit:

            for fit in copy_fit_list:
                if column == 1:
                    if len(fit) == 4 and fit[1] == 'Polymorphous':
                        name = self.varTable.item(row,0).text()
                        if my_layer == fit[0] and element == fit[2] and name == fit[3]:
                            self.parameterFit.remove(fit)
                elif column == 2:
                    if len(fit) == 2:
                        scattering_factor = self.varTable.item(row,2).text()
                        if scattering_factor == fit[1]:
                            self.parameterFit.remove(fit)
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
        elif column == 0:
            _element_fit.setDisabled(True)
            _compound_fit.setDisabled(True)
            _remove_fit.setDisabled(True)

        action = menu.exec_(QtGui.QCursor.pos())

        # Element Mode
        if action == _element_fit:


            element = self.structTable.item(row, 0).text()
            alreadySelected = False

            # Check to make sure parameter is not already selected
            for fit in copy_fit_list:
                # check if layer and parameter in compound mode
                n = len(fit)
                if n == 2:  # scattering factor
                    if column == 5:  # trying to fit scattering factor
                        item = self.structTable.currentItem().text()
                        if item == fit[1]:
                            alreadySelected = True

                elif n == 3:  # compound mode
                    layer = fit[0]
                    param = fit[2]
                    param_num = 1
                    if param == 'Thickness':
                        param_num = 1
                    elif param == 'Density':
                        param_num = 2
                    elif param == 'Roughness':
                        param_num = 3
                    elif param == 'Linked Roughness':
                        param_num = 4

                    if layer == my_layer and column == param_num:
                        self.parameterFit.remove(fit)

                elif n == 4:  # element mode
                    layer = fit[0]
                    ele = fit[2]
                    param = fit[3]
                    my_ele = self.structTableInfo[idx][row][0]

                    param_num = 1
                    if param == 'Thickness':
                        param_num = 1
                    elif param == 'Density':
                        param_num = 2
                    elif param == 'Roughness':
                        param_num = 3
                    elif param == 'Linked Roughness':
                        param_num = 4

                    if layer == my_layer and column == param_num and ele == my_ele:
                        alreadySelected = True

            if not alreadySelected:
                if column == 1:  # thickness
                    self.parameterFit.append([my_layer, 'Element', element, 'Thickness'])
                elif column == 2:  # density
                    self.parameterFit.append([my_layer, 'Element', element, 'Density'])
                elif column == 3:  # roughness
                    self.parameterFit.append([my_layer, 'Element', element, 'Roughness'])
                elif column == 4:  # linked roughness
                    self.parameterFit.append([my_layer, 'Element', element, 'Linked Roughness'])
                elif column == 5:  # scattering factor
                    scattering_factor = self.structTable.item(row, 5).text()
                    if scattering_factor[0] != '[':
                        self.parameterFit.append(['Scattering Factor', scattering_factor])

        elif action == _compound_fit:

            alreadySelected = False
            for fit in copy_fit_list:
                n = len(fit)
                if n == 3:  # compound check
                    layer = fit[0]
                    param = fit[2]

                    param_n = 1
                    if param == 'Thickness':
                        param_n = 1
                    elif param == 'Density':
                        param_n = 2
                    elif param == 'Roughness':
                        param_n = 3
                    elif param == 'Linked Roughness':
                        param_n = 4


                    if layer == my_layer and param_n == column:
                        alreadySelected = True

                elif n == 4:  # element check
                    layer = fit[0]
                    param = fit[3]
                    param_n = 1
                    if param == 'Thickness':
                        param_n = 1
                    elif param == 'Density':
                        param_n = 2
                    elif param == 'Roughness':
                        param_n = 3
                    elif param == 'Linked Roughness':
                        param_n = 4


                    if param_n == column and layer == my_layer:
                        self.parameterFit.remove(fit)

                if not alreadySelected:

                    if column == 1:  # thickness
                        my_fit = [my_layer, 'Compound', 'Thickness']
                        self.parameterFit.append(my_fit)
                    elif column == 2:  # density
                        my_fit = [my_layer, 'Compound', 'Density']
                        self.parameterFit.append(my_fit)
                    elif column == 3:  # roughness
                        my_fit = [my_layer, 'Compound', 'Roughness']
                        self.parameterFit.append(my_fit)
                    elif column == 4:  # linked roughness
                        my_fit = [my_layer, 'Compound', 'Linked Roughness']
                        self.parameterFit.append(my_fit)

        elif action == _remove_fit:
            element = self.structTableInfo[my_layer][row][0]
            scattering_factor = self.structTableInfo[my_layer][row][5]
            for fit in copy_fit_list:
                n = len(fit)
                if column == 1:
                    mode = fit[1]
                    if mode == "Element":
                        ele = fit[2]
                        if ele == element and fit[3] == 'Thickness':
                            self.parameterFit.remove(fit)
                    elif mode == 'Compound':
                        if fit[2] == 'Thickness':
                            self.parameterFit.remove(fit)
                elif column == 2:
                    mode = fit[1]
                    if mode == "Element":
                        ele = fit[2]
                        if ele == element and fit[3] == 'Density':
                            self.parameterFit.remove(fit)
                    elif mode == 'Compound':
                        if fit[2] == 'Density':
                            self.parameterFit.remove(fit)
                elif column == 3:
                    mode = fit[1]
                    if mode == "Element":
                        ele = fit[2]
                        if ele == element and fit[3] == 'Roughness':
                            self.parameterFit.remove(fit)
                    elif mode == 'Compound':
                        if fit[2] == 'Roughness':
                            self.parameterFit.remove(fit)
                elif column == 4:
                    mode = fit[1]
                    if mode == "Element":
                        ele = fit[2]
                        if ele == element and fit[3] == 'Linked Roughness':
                            self.parameterFit.remove(fit)
                    elif mode == 'Compound':
                        if fit[2] == 'Linked Roughness':
                            self.parameterFit.remove(fit)
                elif column == 5 and n == 2:
                    if scattering_factor == fit[1]:
                        self.parameterFit.remove(fit)

        self.setTable(my_layer)

    def changeStepSize(self):
        self._step_size = self.step_size.text()

    def setStructFromSample(self, idx):
        # this function will take the upon initialization and load in all parameters
        structInfo = self.sample.structure
        numLay = len(structInfo)
        tempArray = np.empty((3,7), dtype=object)
        for row in range(3):
            for col in range(7):
                if col == 0:
                    element = list(structInfo[idx].keys())[row]
                    tempArray[row,col] = str(element)
                elif col == 1:
                    element = list(structInfo[idx].keys())[row]
                    thickness = structInfo[idx][element].thickness
                    tempArray[row,col] = str(thickness)
                elif col == 2:
                    element = list(structInfo[idx].keys())[row]
                    density = structInfo[idx][element].density
                    tempArray[row,col] = str(density)
                elif col == 3:
                    element = list(structInfo[idx].keys())[row]
                    roughness = structInfo[idx][element].roughness
                    tempArray[row,col] = str(roughness)
                elif col == 4:
                    element = list(structInfo[idx].keys())[row]
                    linked_roughness = structInfo[idx][element].linked_roughness
                    tempArray[row,col] = str(linked_roughness)
                elif col == 5:
                    element = list(structInfo[idx].keys())[row]
                    scattering_factor = structInfo[idx][element].scattering_factor
                    tempArray[row,col] = str(scattering_factor)
                elif col == 6:
                    element = list(structInfo[idx].keys())[row]
                    stoichiometry = structInfo[idx][element].stoichiometry
                    tempArray[row,col] = str(stoichiometry)


        self.structTableInfo.append(tempArray)

    def setTable(self, idx):
        tableInfo = self.structTableInfo[idx]
        num_rows = self.structTable.rowCount()

        for col in range(6):
            for row in range(num_rows):

                item = QTableWidgetItem(str(tableInfo[row][col]))
                self.structTable.setItem(row,col, item)


        for fit in self.parameterFit:
            layer = fit[0]
            n = len(fit)
            if layer == idx:  # not scattering factor parameters
                if n == 3: # compound mode
                    param = fit[2]
                    param_n = 0

                    if param == 'Thickness':
                        param_n = 1
                    elif param == 'Density':
                        param_n = 2
                    elif param == 'Roughness':
                        param_n = 3
                    elif param == 'Linked Roughness':
                        param_n = 4

                    for row in range(num_rows):
                        if param_n != 0:
                            self.structTable.item(row, param_n).setBackground(QtGui.QColor(150, 150, 255))
                elif n == 4:  # element mode
                    ele = fit[2]
                    param = fit[3]

                    param_n = 0
                    if param == 'Thickness':
                        param_n = 1
                    elif param == 'Density':
                        param_n = 2
                    elif param == 'Roughness':
                        param_n = 3
                    elif param == 'Linked Roughness':
                        param_n = 4

                    for row in range(num_rows):
                        my_ele = self.structTableInfo[idx][row][0]
                        if my_ele == ele:
                            if param_n != 0:
                                self.structTable.item(row, param_n).setBackground(QtGui.QColor(150, 255, 150))
            if layer == 'Scattering Factor':
                for row in range(num_rows):
                    if fit[1] == self.structTableInfo[idx][row][5]:
                        self.structTable.item(row,5).setBackground(QtGui.QColor(150, 255, 150))


    def setTableVar(self):
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
                            if len(fit) == 2 and fit[1] == info[2][row]:
                                self.varTable.item(row,2).setBackground(QtGui.QColor(150, 255, 150))
                            elif len(fit) == 4:
                                if fit[0] == layer_idx and fit[1] == "Polymorphous" and fit[2] == ele and fit[3] == info[0][row]:
                                    self.varTable.item(row,1).setBackground(QtGui.QColor(150, 255, 150))
            else:
                for row in range(self.varTable.rowCount()):
                    item = QTableWidgetItem('')
                    self.varTable.setItem(row, 0, item)

                    item = QTableWidgetItem('')
                    self.varTable.setItem(row, 1, item)

                    item = QTableWidgetItem('')
                    self.varTable.setItem(row, 2, item)





    def setTableMag(self):

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

            mydensity = density[row]
            mysf = sf[row]
            if mydensity == '':
                item = QTableWidgetItem('')
                self.magTable.setItem(row,0,item)
            else:
                item = QTableWidgetItem(str(mydensity))
                self.magTable.setItem(row, 0, item)

            if mysf == '':
                item = QTableWidgetItem('')
                self.magTable.setItem(row, 1, item)
            else:
                item = QTableWidgetItem(mysf)
                self.magTable.setItem(row, 1, item)


    def changetable(self):

        idx = self.layerBox.currentIndex()

        if self.structBool:
            for row in range(len(self.structTableInfo[self.previousLayer])):
                for col in range(len(self.structTableInfo[self.previousLayer][0])-1):
                    item = self.structTable.item(row,col).text()
                    self.structTableInfo[self.previousLayer][row][col] = str(item)

            self.previousLayer = idx


        if self.polyBool:
            # get information from variation table
            if not self.change_elements:
                self.element_index = self.elementBox.currentIndex()

            ele = self.structTableInfo[self.previousLayer][self.element_index][0]

            name = ['' for i in range(self.varTable.rowCount())]  # initialize for no input case
            ratio = ['' for i in range(self.varTable.rowCount())]
            scat = ['' for i in range(self.varTable.rowCount())]

            for row in range(self.varTable.rowCount()):
                for col in range(self.varTable.columnCount()):
                    if self.varTable.item(row, col) == None:
                        pass
                    else:
                        item = self.varTable.item(row, col).text()
                        if col == 0:
                            name[row] = item
                        elif col == 1:
                            ratio[row] = item
                        elif col == 2:
                            scat[row] = item

            # Makes sure that if the name or scattering factor is changed, that we change it throughout
            for i in range(len(self.structTableInfo)):
                for j in range(len(self.structTableInfo[i])):

                    new_ele = self.structTableInfo[i][j][0]
                    if new_ele == ele and i != self.previousLayer:

                        self.varData[ele][i][0] = name
                        self.varData[ele][i][2] = scat

                        if self.magData[ele][i][0][0] != '' and len(self.magData[ele][i][0]) != 1:
                            self.magData[ele][i][0] = name  # makes sure that magnetic components has correct names



            self.varData[ele][self.previousLayer] = [name, ratio, scat]
            if self.magData[ele][self.previousLayer][0][0] != '' and len(self.magData[ele][self.previousLayer][0]) != 1:
                self.magData[ele][self.previousLayer][0] = name  # makes sure that magnetic components has correct names

            self.previousLayer = idx

        if self.magBool:
            layer = self.structTableInfo[self.previousLayer]

            elements = []

            for ele_idx in range(len(layer)):
                ele = layer[ele_idx][0]
                elements.append(ele)


            e = 0  # element index
            v = 0 # element variation index
            for row in range(self.magTable.rowCount()):

                element = elements[e]  # gets the proper element


                names = self.magData[element][self.previousLayer][0]
                num_v = len(names)

                for col in range(self.magTable.columnCount()):
                    item = self.magTable.item(row,col).text()  # gets the current item

                    if num_v == 1:
                        if col == 0:
                            self.magData[element][self.previousLayer][1][0] = item
                        elif col == 1:
                            self.magData[element][self.previousLayer][2][0] = item
                            e = e + 1

                    else:

                        if col == 0:
                            self.magData[element][self.previousLayer][1][v] = item
                        elif col == 1:
                            self.magData[element][self.previousLayer][2][v] = item
                            v = v + 1
                            if v > num_v-1:
                                v = 0
                                e = e + 1

                        # Makes sure that if the name or scattering factor is changed, that we change it throughout
            for i in range(len(self.structTableInfo[self.previousLayer])):
                element = self.structTableInfo[self.previousLayer][i][0]
                new_sf = self.magData[element][self.previousLayer][2]
                for lay_idx in range(len(self.magData[element])):

                    if self.magData[element][lay_idx][0][0] != '':
                        self.magData[element][lay_idx][2] = new_sf
                        if self.magData[element][lay_idx][1][0] == '' and new_sf[0] != '':
                            self.magData[element][lay_idx][1][0] = 0
                        elif new_sf[0] == '':
                            self.magData[element][lay_idx][1][0] = ''

            self.previousLayer = idx

        self.setTable(idx)



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
            my_elements.append(userinput[i][0])
            if userinput[i][0] not in list(self.varData.keys()):
                self.varData[userinput[i][0]] = [[['',''],['',''],['','']] for i in range(num_layers)]
                self.magData[userinput[i][0]] = [[[userinput[i][0]],[''],['']] for i in range(num_layers)]


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

        self.setTable(idx)  # sets the table for layer that replaces the removed layer

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
        self.changetable()

        self.structBool = True
        self.polyBool = False
        self.magBool = False

        self.changetable()

        self.sampleInfoLayout.setCurrentIndex(0)

    def _elementVariation(self):
        self.changetable()

        self.structBool = False
        self.polyBool = True
        self.magBool = False

        self.changetable()

        self.sampleInfoLayout.setCurrentIndex(1)

    def _magnetic(self):
        self.changetable()

        self.structBool = False
        self.polyBool = False
        self.magBool = True

        self.changetable()

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

    def _createSample(self):

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
                if element[4].isdigit():
                    linked_roughness.append(float(element[4]))
                else:
                    linked_roughness.append(False)

                # still need to take into account sf that are different than element
            sample.addlayer(idx,formula,thickness,density=density, roughness=roughness, linked_roughness=linked_roughness)

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
        self.bs = []  # background shift
        self.sf = []  # scaling factor

        self.sWidget = sWidget
        self.sample = sWidget.sample
        self.data = data
        self.data_dict = data_dict

        self.fit = []  # scans to be fit
        self.bounds = []  # bounds
        self.weights = []  # weights of the bounds
        self.axis_state = True
        self.scan_state = True
        self.previousIdx = 0

        # Adding the plotting Widget
        self.spectrumWidget = pg.PlotWidget()
        self.spectrumWidget.setBackground('w')

        self.spectrumWidget.addLegend()
        # This will be used to determine which scan to view
        whichScanLabel = QLabel('Scans:')
        self.whichScan = QComboBox()
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

        self.boundWeightTable = QTableWidget()

        self.selectedScans = QComboBox()

        self.selectedScans.activated.connect(self.changeColorFit)
        self.selectedScans.activated.connect(self.readTable)
        self.whichScan.activated.connect(self.changeColorScan)
        self.selectedScans.activated.connect(self.myPlotting)
        self.selectedScans.activated.connect(self.setTable)

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
        self.opButton = QPushButton('Optical Profile')
        self.opButton.setStyleSheet('background: grey')
        self.opButton.clicked.connect(self.opPlot)
        self.opmButton = QPushButton('Magneto-Optical Profile')
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

        scanSelectionLayout = QVBoxLayout()
        scanSelectionLayout.addLayout(whichScanLayout)
        scanSelectionLayout.addWidget(self.fitButton)
        scanSelectionLayout.addStretch(1)
        scanSelectionLayout.addLayout(hbox)
        scanSelectionLayout.addLayout(stepLayout)
        scanSelectionLayout.addStretch(1)
        scanSelectionLayout.addWidget(self.rButton)
        scanSelectionLayout.addWidget(self.opButton)
        scanSelectionLayout.addWidget(self.opmButton)
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

        sfbsLayout.addLayout(allScansLayout)
        sfbsLayout.addLayout(bsLayout)
        sfbsLayout.addLayout(sfLayout)

        semiBotLayout = QHBoxLayout()
        semiBotLayout.addLayout(sfbsLayout)
        semiBotLayout.addWidget(self.boundWeightTable)


        bottomlayout = QHBoxLayout()
        bottomlayout.addLayout(semiBotLayout)
        bottomlayout.addWidget(boundWidget)

        pagelayout = QVBoxLayout()
        pagelayout.addLayout(toplayout)
        pagelayout.addLayout(bottomlayout)

        self.setLayout(pagelayout)

    def changeSFandBS(self):
        idx = self.selectedScans.currentIndex()
        state = self.allScan.checkState()

        bs = self.backgroundShift.text()
        sf = self.scalingFactor.text()
        if bs != '' and sf != '':
            if self.allScan.checkState() == 0:
                self.bs[idx] = bs
                self.sf[idx] = sf
            else:
                for i in range(len(self.bs)):
                    self.bs[i] = bs
                    self.sf[i] = sf


    def changeStepSize(self):
        self.sWidget._step_size = self.stepWidget.text()

    def myPlotting(self):
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
        name = ''
        self.spectrumWidget.clear()

        if self.scan_state:
            name = self.whichScan.currentText()
        else:
            name = self.selectedScans.currentText()

        E = self.data_dict[name]['Energy']

        step_size = float(self.sWidget._step_size)
        thickness, density, density_magnetic = self.sample.density_profile(step=step_size)  # Computes the density profile

        sf = dict()  # form factors of non-magnetic components
        sfm = dict()  # form factors of magnetic components

        # Non-Magnetic Scattering Factor
        for e in self.sample.find_sf[0].keys():
            sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E, False)
        # Magnetic Scattering Factor
        for em in self.sample.find_sf[1].keys():
            sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E, True)

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
        name = ''
        self.spectrumWidget.clear()

        if self.scan_state:
            name = self.whichScan.currentText()
        else:
            name = self.selectedScans.currentText()

        E = self.data_dict[name]['Energy']

        step_size = float(self.sWidget._step_size)
        thickness, density, density_magnetic = self.sample.density_profile(step=step_size)  # Computes the density profile

        sf = dict()  # form factors of non-magnetic components
        sfm = dict()  # form factors of magnetic components

        # Non-Magnetic Scattering Factor
        for e in self.sample.find_sf[0].keys():
            sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E, False)
        # Magnetic Scattering Factor
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


    def changeColorScan(self):
        self.selectedScans.setStyleSheet('background: white; selection-background-color: grey')
        self.whichScan.setStyleSheet('background: red; selection-background-color: red')
        self.scan_state = True

    def changeColorFit(self):
        self.selectedScans.setStyleSheet('background: red; selection-background-color: red')
        self.whichScan.setStyleSheet('background: white; selection-background-color: grey')
        self.scan_state = False

    def setTable(self):

        idx = self.selectedScans.currentIndex()
        self.scalingFactor.setText(self.sf[idx]) # setting the appropriate scaling factor

        self.backgroundShift.setText(self.bs[idx])  # setting the appropriate background shift
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

            idx = self.selectedScans.currentIndex()
            if self.allScan.checkState() == 0:
                self.bs.append('0')  # background shift
                self.sf.append('1')  # scaling factor
            else:
                if len(self.bs) != 0 and len(self.bs) != 1:
                    if idx == 0:
                        self.bs.append(self.bs[idx+1])
                        self.sf.append(self.sf[idx+1])
                    else:
                        self.bs.append(self.bs[idx - 1])
                        self.sf.append(self.sf[idx - 1])
                elif len(self.bs) == 1:
                    self.bs.append(self.bs[0])
                    self.sf.append(self.sf[0])
                else:
                    self.bs.append('0')  # background shift
                    self.sf.append('1')  # scaling factor

    def _removeScanSelection(self):
        idx = self.selectedScans.currentIndex()

        # takes care of proper indexing
        if idx == self.previousIdx and self.previousIdx != 0:
            self.previousIdx = self.previousIdx - 1

        self.selectedScans.removeItem(idx)
        self.fit.pop(idx)
        self.bounds.pop(idx)
        self.weights.pop(idx)

        self.bs.pop(idx)  # background shift
        self.sf.pop(idx)  # scaling factor

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
        idx = self.selectedScans.currentIndex()

        if col != 1:
            n = len(self.bounds[idx])
            upper = self.bounds[idx][n - 1][1]
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

        self.sample = self.sWidget.sample
        self.spectrumWidget.clear()
        idx = self.whichScan.currentIndex()
        name = self.whichScan.currentText()
        dat = self.data_dict[name]['Data']
        pol = self.data_dict[name]['Polarization']
        scan_type = self.data[idx][1]
        step_size = float(self.sWidget._step_size)

        if scan_type == 'Reflectivity':
            qz = dat[0]

            R = dat[2]
            E = self.data_dict[name]['Energy']
            qz, Rsim = self.sample.reflectivity(E,qz, s_min=step_size)
            Theta = np.arcsin(qz / (E * 0.001013546143))*180/np.pi
            Rsim = Rsim[pol]
            if pol == 'S' or pol =='P' or pol =='LC' or pol == 'RC':

                if self.axis_state:
                    self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                    self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                else:
                    self.spectrumWidget.plot(Theta,R,pen=pg.mkPen((0,2), width=2), name='Data')
                    self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                self.spectrumWidget.setLabel('left', "Reflectivity, R")
                self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                self.spectrumWidget.setLogMode(False,True)
            elif pol == 'AL' or pol =='AC':
                if self.axis_state:
                    self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                    self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                else:
                    self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                    self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                self.spectrumWidget.setLogMode(False, False)
                self.spectrumWidget.setLabel('left', "Reflectivity, R")
                self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
        elif scan_type == 'Energy':
            E = dat[3]
            R = dat[2]
            Theta = self.data_dict[name]['Angle']
            E, Rsim = self.sample.energy_scan(Theta,E, s_min=step_size)
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
                qz, Rsim = self.sample.reflectivity(E,qz, s_min=step_size)
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
                    if self.axis_state:
                        self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                    else:
                        self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                    self.spectrumWidget.setLogMode(False, False)
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.spectrumWidget.setXRange(lower, upper)
            else:
                E = dat[3]
                R = dat[2]
                Theta = self.data_dict[name]['Angle']
                E, Rsim = self.sample.energy_scan(Theta,E, s_min=step_size)
                Rsim = Rsim[pol]
                self.spectrumWidget.setLogMode(False, False)
                self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                self.spectrumWidget.plot(E, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                self.spectrumWidget.setLabel('left', "Reflectivity, R")
                self.spectrumWidget.setLabel('bottom', "Energy, E (eV)")

class GlobalOptimization(QMainWindow):
    def __init__(self):
        super().__init__()

class ReflectometryApp(QMainWindow):
    def __init__(self, fname):
        super().__init__()
        sample = ds.ReadSampleHDF5(fname)  # temporary way of loading in data
        data, data_dict, sim_dict = ds.ReadDataHDF5(fname)
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

        self._sampleWidget = sampleWidget(sample)  # initialize the sample widget
        self._reflectivityWidget = reflectivityWidget(self._sampleWidget, data, data_dict, sim_dict)

        self.sampleButton.setStyleSheet("background-color : pink")
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
        self.stackedlayout.addWidget(label3)

        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)

    def activate_tab_1(self):
        self._sampleWidget.step_size.setText(self._sampleWidget._step_size)
        self.stackedlayout.setCurrentIndex(0)
    def activate_tab_2(self):
        self._reflectivityWidget.stepWidget.setText(self._sampleWidget._step_size)
        self.stackedlayout.setCurrentIndex(1)
    def activate_tab_3(self):
        self.stackedlayout.setCurrentIndex(2)





if __name__ == '__main__':


    fname = 'Pim10uc.h5'

    app = QApplication(sys.argv)
    demo = ReflectometryApp(fname)

    demo.show()

    sys.exit(app.exec_())