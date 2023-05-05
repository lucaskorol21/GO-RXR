"""
Library: GUI_GO
Version: 0.2
Author: Lucas Korol
Institution: University of Saskatchewan
Last Updated: March 28nd, 2023
Python: version 3.7

Purpose: This python file contains the graphical user interface for the GO-RXR software package.

---------------------------------------------------------------------------------------------------------------------------------------
Imported Libraries

material_model (version 0.1) - Part of the 'name' software package and is used to retrieve the form factors and
                               calculate the optical constants

material_structure (version 0.1) - Part of the 'name' software package and is used to retrieve to calculate the reflectivity spectra

global_optimization (version 0.1) - Part of the 'name' software package and is used to perform the global optimization

data_structure (version 0.2) - Part of the 'name' software package that is used to saving and loading the workspace

numpy (version 1.21.4) - used for array manipulation

scipy (version 1.7.1) - used for data smoothing

PyQt5 (version 5.15.7) - This library is used to contruct the application

pyqtgraph (version 0.12.4) - This library is used for the plotting widget

pyinstaller (version 5.9.0) - This library is used to create a standalone executable file


--------------------------------------------------------------------------------------------------------------------------------------
Note: There are a few immediate changes that could be implemented for future versions which are listed below:

1. Magnetization direction. Currently, we can only calculate the reflectivity spectra for a magnetization direction
   in the x, y, and z-directions. In the case where we can have an arbitarty direction that are defined by phi and theta
   this should be altered. Instead of using a QComboBox in magneticWidget for the magnetization direction I would create
   two QLineEdits for phi and theta respectively. I would make sure to include checks that the angles are physically
   appropriate. Note that if this is done changes would need to be made in the magneticWidget, sampleWidget, and the
   data_structure python file (save the magnetization direction as phi and theta instead of a string x,y,z). A good plan
   to make sure you have changed everything accordingly would be to search magDirection in the python file.

2. Layer magnetization. Currently, the user has the ability to give different magnetization directions per layer.
   The current reflectivity calculation implementation is unable to handle multiple magnetization directions, so
   it would make sense to remove this capability. If I have time I will get to this.

3. Inclusion of other global optimization algorithms. I would like to include the direct algorithm into the the list of
   algorithms to use. The issue is that we required python 3.8 and above to use it, but some of the underlining
   code that is used for the reflectivity calculations does not allow for the use of python 3.8 and above. I've already
   included a lot of the code to include the direct algorithm, or any other algorithm. I would suggest searching 'direct'
   in all of the python files (or another algorithm) to see where to make the appropriate changes.

4. Data smoothing. It may be worth including other data smoothing methods. For example, there has been discussion of
   using neural networks to perform some of the smoothing. The neural networks have been found to remove the noise
   while maintaining the shape, even for very noisy signals.

Warning: Any changes to the data type would need to be carefully considered.

Instructions to create executable using pytinstaller:

1. Install pyinstaller in the python environement.
2. Run 'pyinstaller GUI.GO.py' in the terminal.
3. A new file named GUI_GO.spec will appear in the project file. Open this file.
4. In this file there is a filed called datas = []. Replace the empty array with:
    datas = [('.\\global_optimization.py', '.'), ('.\\data_structure.py','.'),('.\\material_structure.py','.'),
    ('.\\material_model.py','.'), ('.\\default_script.txt','.'), ('.\\form_factor.pkl','.'),
    ('.\\form_factor_magnetic.pkl','.'), ('.\\Perovskite_Density.txt','.'), ('.\\Atomic_Mass.txt','.'),
    ('.\\demo.h5','.')]
5. This includes all the python files, text files, and all other data used to execute GUI_GO.py. If newer versions
   of GO-RXR depend on more files than they must be included into this array. Follow the same naming convention as shown.
6. Once satisfied run 'pyinstaller GUI_GO.spec' in the terminal. This will ensure that all desired data is included into
   the data file.
7. The executable file will be found in the 'dist' directory in the project workspace. There should now be a directory
   in this workspace with the name GUI_GO.
8. You can copy this file an input it into another directory. From here I would suggest creating a file zip file that
   can be distributed.
9. A one file that can compile the executable into a single file instead of a directory. However, this can take much longer
   to initialize as with the one file all the libraries and files need to be unpacked before the software can be run.
"""

import ast
from scipy import interpolate
from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import *
from PyQt5.QtGui import *
import traceback
import numpy as np
import time
import sys
import material_structure as ms
import os
import pyqtgraph as pg
import data_structure as ds
import copy
import global_optimization as go

import material_model as mm
from scipy import signal
import h5py
import multiprocessing as mp
import pickle
from scipy.interpolate import interp1d, UnivariateSpline
from scipy.fft import fft, fftfreq, fftshift, ifft, ifftshift


global stop
stop = False

# keep track of the current state of the global optimization
global x_vars
x_vars = []


def stringcheck(string):
    """
    Purpose: Checks to make sure sample properties are in the correct format
    :param string: Parameters value
    :return: Boolean determining if value is in the correct format
    """
    num = 0  # used to check how many '.' are found in the string
    correctFormat = True  # is value in the correct format?
    if len(string) == 0:  # string is empty
        correctFormat = False
    else:
        first = True  # checks the first character in the string
        for char in string:
            if first:  # makes sure first character is a digit
                if char == '.':
                    correctFormat = False
                elif not char.isdigit():
                    correctFormat = False
                first = False

            elif not char.isdigit():  # checks if all other characters are digits
                if char == '.':
                    num = num + 1
                    if num > 1:
                        correctFormat = False
                else:
                    correctFormat = False

    return correctFormat


class compoundInput(QDialog):
    """
    Purpose: Creates a widget when user wants to add another layer. This widget allows the user to select the chemical
             formula, thickness, density (g/cm^3), roughness, and linked roughness for the layer.
    """

    def __init__(self):
        super().__init__()
        self.val = []  # class value used to store layer properties
        pagelayout = QVBoxLayout()  # page layout
        infolayout = QGridLayout()  # compound information layout

        # Chemical formula of added layer
        formula = QLabel('Formula: ')
        self.formula = QLineEdit()
        self.formula.editingFinished.connect(self.formulaDone)

        # thickness of layer
        thickness = QLabel('Thickness (A): ')
        self.thickness = QLineEdit()
        self.thickness.setText('10')

        # Density (g/cm^3)
        density = QLabel('Density (g/cm^3): ')
        self.density = QLineEdit()

        # roughness of layer
        roughness = QLabel('Roughness (A): ')
        self.roughness = QLineEdit()
        self.roughness.setText('2')

        # include linked roughness if applicable to layer
        linkedroughnesslayout = QHBoxLayout()  # layout for linked roughness
        linkedroughness = QLabel('Linked Roughness (A): ')
        self.linkedroughness = QLineEdit()
        self.linkedroughness.setHidden(True)

        # checkbox to determine if user wants to add a linked roughness to layer
        self.checkbox = QCheckBox()
        self.checkbox.stateChanged.connect(self.linkedroughnessbox)
        self.checkboxstate = 0

        linkedroughnesslayout.addWidget(self.checkbox)
        linkedroughnesslayout.addWidget(linkedroughness)

        # add labels to the info layout
        infolayout.addWidget(formula, 0, 0)
        infolayout.addWidget(thickness, 1, 0)
        infolayout.addWidget(density, 2, 0)
        infolayout.addWidget(roughness, 3, 0)
        infolayout.addLayout(linkedroughnesslayout, 4, 0)

        # add the text edit widgets to the info layout
        infolayout.addWidget(self.formula, 0, 1)
        infolayout.addWidget(self.thickness, 1, 1)
        infolayout.addWidget(self.density, 2, 1)
        infolayout.addWidget(self.roughness, 3, 1)
        infolayout.addWidget(self.linkedroughness, 4, 1)

        # add layer properties to the current sample model
        enterButton = QPushButton('Enter')
        enterButton.clicked.connect(self.inputComplete)

        self.errorMessage = QLabel('')  # let user know if they have made a mistake in their input

        # add all layouts to the page layout
        pagelayout.addLayout(infolayout)
        pagelayout.addWidget(enterButton)
        pagelayout.addWidget(self.errorMessage)
        self.setLayout(pagelayout)

    def formulaDone(self):
        """
        Purpose: searches database for material density
        :return: density (g/cm^3)
        """
        cwd = os.getcwd()
        filename = 'Perovskite_Density.txt'

        found = False  # boolean user that determines if material found in database
        with open(filename) as file:
            for line in file:
                myformula = line.split()[0]  # retrieves the chemical formula
                mydensity = line.split()[1]  # retrieves the density
                if not found:
                    if self.formula.text() == myformula:  # material in database
                        self.density.setText(mydensity)  # set the density in the widget
                        found = True
                    else:  # material not found in the database
                        self.density.clear()  # reset the density in the widget

    def linkedroughnessbox(self):
        """
        Purpose: Hide or make visible the linked roughness lineEdit widget
        """
        self.checkboxstate = self.checkbox.checkState()  # retrieve the checkbox state
        if self.checkboxstate > 0:
            self.linkedroughness.setHidden(False)
        else:
            self.linkedroughness.setHidden(True)

    def inputComplete(self):
        """
        Purpose: Checks to make sure parameters are in correct format before sending back
                 to main widget to add to current sample
        :return: List of the layer parameters
        """

        finished = True
        # gets the elements and their stoichiometry
        myElements = ms.find_stoichiometry(self.formula.text())  # gets the elements and their stoichiometry

        # gets the density
        myThickness = self.thickness.text()
        thicknessCorrect = stringcheck(myThickness)  # checks thickness format

        myDensity = self.density.text()
        densityCorrect = stringcheck(myDensity)  # checks density format

        # gets the density
        myRoughness = self.roughness.text()
        roughnessCorrect = stringcheck(myRoughness)  # checks roughness format

        # gets the linked roughness
        myLinkedroughness = self.linkedroughness.text()

        # sends error message to user if formate is incorrect
        linkedroughnessCorrect = True
        if myLinkedroughness != '':
            linkedroughnessCorrect = stringcheck(myLinkedroughness)

        if not (thicknessCorrect) or not (densityCorrect) or not (roughnessCorrect) or not (linkedroughnessCorrect):
            if not (thicknessCorrect):
                self.errorMessage.setText('Please check thickness!')
            elif not (densityCorrect):
                self.errorMessage.setText('Please check density!')
            elif not (roughnessCorrect):
                self.errorMessage.setText('Please check roughness!')
            elif not (linkedroughnessCorrect):
                self.errorMessage.setText('Please check linked roughness!')
        else:  # transform sample parameters into correct format
            molar_mass = 0  # molar mass
            elements = list(myElements[0].keys())  # list of elements in layer
            # gets the molar mass of the compound
            for ele in elements:
                stoich = myElements[0][ele].stoichiometry

                if ele[-1].isdigit():
                    ele = ele.rstrip(ele[-1])

                molar_mass = molar_mass + ms.atomic_mass(ele) * stoich

            tempArray = []

            # put layer info in correct format depending for both linked roughness cases
            if myLinkedroughness == '':  # no linked roughness case
                # pre-set form factor to element name
                for ele in elements:
                    stoich = myElements[0][ele].stoichiometry
                    density = float(myDensity)
                    if ele[-1].isdigit():
                        ff_ele = ele.rstrip(ele[-1])
                    else:
                        ff_ele = ele

                    tempArray.append(
                        [ele, myThickness, str(density * float(stoich) / molar_mass), myRoughness, False, ff_ele,
                         stoich])
            else:  # linked roughness case
                for ele in elements:
                    stoich = myElements[0][ele].stoichiometry
                    density = float(myDensity)
                    # pre-set form factor to element symbol
                    if ele[-1].isdigit():
                        ff_ele = ele.rstrip(ele[-1])
                    else:
                        ff_ele = ele

                    tempArray.append(
                        [ele, myThickness, str(density * float(stoich) / molar_mass), myRoughness, myLinkedroughness,
                         ff_ele, stoich])

            self.val = np.array(tempArray)  # reset layer values
            self.accept()  # close the widget


class variationWidget(QDialog):
    """
    Purpose: This widget provides the user a workspace to handle elemental variations.It allows for the user to include
             the different oxidation states of a material. The only properites that the user may change is the number
             of element variations, their ratio, and their form factors.
    """

    def __init__(self, mainWidget, sample):
        """
        Purpose: Initialize the widget
        :param mainWidget: This is the sampleWidget. This is used as an input value so the information stored in sample
                           widget can be used in the variationWidget.
        :param sample: This is the most recent material model (slab class)
        """
        super().__init__()

        pagelayout = QHBoxLayout()  # page layout

        self.elelayout = QVBoxLayout()  # This is the element layout
        self.mainWidget = mainWidget  # referring to the structuralWidget
        self.mainWidget.layerBox.currentIndexChanged.connect(self.changeElements)  # change elements when layer changed

        # add element variation
        addButton = QPushButton('Add')
        addButton.clicked.connect(self.addVarEle)
        # remove element variation
        deleteButton = QPushButton('Delete')
        deleteButton.clicked.connect(self.deleteVarEle)

        self.sample = sample  # reset sample information

        #self.radiobutton = QRadioButton()
        idx = self.mainWidget.layerBox.currentIndex()  # retrieves the current layer index

        # adding elements in the current layer to the combobox
        for j in range(len(list(self.sample.structure[idx].keys()))):
            ele = list(self.sample.structure[idx].keys())[j]
            self.mainWidget.elementBox.addItem(ele)

        self.mainWidget.elementBox.currentIndexChanged.connect(self.mainWidget.setTableVar)  # element selection changed
        self.elelayout.addWidget(self.mainWidget.elementBox)  # add combobox to layout

        self.mainWidget.setTableVar()  # set the element variation table

        # add the buttons to the element layout
        self.elelayout.addWidget(addButton)
        self.elelayout.addWidget(deleteButton)

        # setting the headers for the element variation table
        self.mainWidget.varTable.setRowCount(2)
        self.mainWidget.varTable.setColumnCount(3)
        self.mainWidget.varTable.setHorizontalHeaderLabels(
            ['Name', 'Ratio', 'Form Factor'])

        aFont = QtGui.QFont()
        aFont.setBold(True)
        self.mainWidget.varTable.horizontalHeader().setFont(aFont)

        pagelayout.addLayout(self.elelayout)  # add element layout to the page layout
        pagelayout.addWidget(self.mainWidget.varTable)  # add element variation table to page layout

        self.mainWidget.varTable.setStyleSheet('background-color: white;')
        self.setStyleSheet('background-color: lightgrey;')
        self.setLayout(pagelayout)

    def changeElements(self):
        """
        Purpose: changes the elements based on current layer
        """

        # prevents adding elements to combobox from triggering other signals
        self.mainWidget.change_elements = True  # are we currently changing the elements

        idx = self.mainWidget.layerBox.currentIndex()  # retrieves the current layer
        self.mainWidget.elementBox.clear()  # clears the element comboBox

        # adds the elements found in the current layer
        for j in range(len(self.mainWidget.structTableInfo[idx])):
            ele = self.mainWidget.structTableInfo[idx][j][0]
            self.mainWidget.elementBox.addItem(ele)

        self.mainWidget.change_elements = False  # no longer changing the elements
        self.mainWidget.elementBox.setCurrentIndex(self.mainWidget.element_index) # sets current index to previous index

    def addVarEle(self):
        """
        Purpose: Make appropriate changes to adding new element variation
        :return:
        """

        self.mainWidget.parameterFit = []  # resets fitting parameters
        self.mainWidget.currentVal = []  # resets current fitting parameter values

        current_layer = self.mainWidget.layerBox.currentIndex()  # retrieves the current layer
        current_element = self.mainWidget.elementBox.currentIndex()  # retrieves the current element

        element = self.mainWidget.structTableInfo[current_layer][current_element][0]  # retrieves the element symbol

        # adds an empty spot in the variation data and magnetic data
        row = len(self.mainWidget.varData[element][current_layer][0])
        for lay in range(len(self.mainWidget.varData[element])):
            if type(self.mainWidget.varData[element][lay][0]) == np.ndarray:  # element already initialized for element variation
                self.mainWidget.varData[element][lay][0] = np.append(self.mainWidget.varData[element][lay][0], '')
                self.mainWidget.varData[element][lay][1] = np.append(self.mainWidget.varData[element][lay][1], '')
                self.mainWidget.varData[element][lay][2] = np.append(self.mainWidget.varData[element][lay][2], '')

                self.mainWidget.magData[element][lay][0] = np.append(self.mainWidget.magData[element][lay][0], '')
                self.mainWidget.magData[element][lay][1] = np.append(self.mainWidget.magData[element][lay][1], '')
                self.mainWidget.magData[element][lay][2] = np.append(self.mainWidget.magData[element][lay][2], '')
            else:  # newly added element variation
                self.mainWidget.varData[element][lay][0].append('')  # add another element to name list
                self.mainWidget.varData[element][lay][1].append('')  # add another element to name list
                self.mainWidget.varData[element][lay][2].append('')  # add another element to name list

                self.mainWidget.magData[element][lay][0].append('')  # make appropriate changes to magnetic data
                self.mainWidget.magData[element][lay][1].append('')
                self.mainWidget.magData[element][lay][2].append('')

        # row = self.mainWidget.varTable.rowCount()
        self.mainWidget.varTable.setRowCount(row + 1)  # reset the number of rows in the element variation table

    def deleteVarEle(self):
        """
        Purpose: Remove the last element variation
        """

        # resets the fitting parameters
        self.mainWidget.parameterFit = []
        self.mainWidget.currentVal = []

        current_layer = self.mainWidget.layerBox.currentIndex()  # retrieves the current layer
        current_element = self.mainWidget.elementBox.currentIndex()  # retrieves the current element

        element = self.mainWidget.structTableInfo[current_layer][current_element][0]  # retrieves element symbol

        row = len(self.mainWidget.varData[element][current_layer][0])  # retrieves the current number of rows

        # removes the element variation from the variation data and magnetic data
        if row != 2:  # no longer a polymorphous element
            for lay in range(len(self.mainWidget.varData[element])):
                if type(self.mainWidget.varData[element][lay][0]) == np.ndarray:
                    self.mainWidget.varData[element][lay][0] = self.mainWidget.varData[element][lay][0][:-1]
                    self.mainWidget.varData[element][lay][1] = self.mainWidget.varData[element][lay][1][:-1]
                    self.mainWidget.varData[element][lay][2] = self.mainWidget.varData[element][lay][2][:-1]

                    self.mainWidget.magData[element][lay][0] = self.mainWidget.magData[element][lay][0][:-1]
                    self.mainWidget.magData[element][lay][1] = self.mainWidget.magData[element][lay][1][:-1]
                    self.mainWidget.magData[element][lay][2] = self.mainWidget.magData[element][lay][2][:-1]
                else: # still a polymorphous element
                    self.mainWidget.varData[element][lay][0].pop()  # add another element to name list
                    self.mainWidget.varData[element][lay][1].pop()  # add another element to name list
                    self.mainWidget.varData[element][lay][2].pop()  # add another element to name list

                    self.mainWidget.magData[element][lay][0].pop()  # make changes to magnetic data
                    self.mainWidget.magData[element][lay][1].pop()
                    self.mainWidget.magData[element][lay][2].pop()

            self.mainWidget.varTable.setRowCount(row - 1)  # set the variation table with the correct number of rows


class ReadOnlyDelegate(QStyledItemDelegate):
    # This class is used to set an item delegate to read only
    def createEditor(self, parent, option, index):
        return


class magneticWidget(QDialog):
    """
    Purpose: This widget is used to set the magnetic info based on the user's input
    """
    def __init__(self, mainWidget, sample):
        """
        :param mainWidget: This is the sampleWidget which allows for access of information from this widget
        :param sample: This is the most recent sample model (slab class)
        """
        super().__init__()

        pagelayout = QHBoxLayout()  # page layout

        self.mainWidget = mainWidget  # sampleWidget

        self.sample = sample  # resets sample

        idx = self.mainWidget.layerBox.currentIndex()  # retrieves current index
        self.mainWidget.layerBox.currentIndexChanged.connect(self.mainWidget.setTableMag)  # set table signal

        # Magnetization direction Widget format
        magLabel = QLabel('Magnetization Direction')
        magLayout = QVBoxLayout()

        # magnetization direction (use phi and theta in future versions)
        # - instead of using a combobox consider using two QLineEdit widgets (one for phi and one for theta)
        #
        self.mainWidget.magDirBox.addItem('x-direction')
        self.mainWidget.magDirBox.addItem('y-direction')
        self.mainWidget.magDirBox.addItem('z-direction')

        # magnetic layout
        magLayout.addWidget(magLabel)
        magLayout.addWidget(self.mainWidget.magDirBox)
        magLayout.addStretch(1)

        # magnetic direction signal setup
        self.mainWidget.magDirBox.currentIndexChanged.connect(self.magDirectionChange)

        # pre-sets the magTable and its headers
        self.mainWidget.magTable.setRowCount(3)
        self.mainWidget.magTable.setColumnCount(2)
        self.mainWidget.magTable.setHorizontalHeaderLabels(
            ['Magnetic Density (mol/cm^3)', 'Form Factor'])

        afont = QtGui.QFont()
        afont.setBold(True)
        self.mainWidget.magTable.horizontalHeader().setFont(afont)
        self.mainWidget.magTable.verticalHeader().setFont(afont)

        self.mainWidget.setTableMag()  # set magTable

        # set the page layout
        pagelayout.addWidget(self.mainWidget.magTable)
        pagelayout.addLayout(magLayout)
        self.setLayout(pagelayout)

        # setting the
        self.mainWidget.magTable.setStyleSheet('background-color: white;')
        self.setStyleSheet('background-color: lightgrey;')


    def magDirectionChange(self):
        """
        Purpose: change the magnetization direction for sampleWidget
        :return:
        """
        #lay = self.mainWidget.layerBox.currentIndex()  # retrieves layer
        mag = self.mainWidget.magDirBox.currentIndex()  # retrieves mag-direction index
        m = len(self.mainWidget.structTableInfo)

        # changes the magnetization direction for all layer
        for i in range(m):
            if mag == 0:
                self.mainWidget.magDirection[i] = 'x'
            elif mag == 1:
                self.mainWidget.magDirection[i] = 'y'
            elif mag == 2:
                self.mainWidget.magDirection[i] = 'z'



class sampleWidget(QWidget):
    """
    Purpose: This widget contains all the information about the sample and all it's properties.
    """
    # sample widget that contains all the information about the sample parameters
    def __init__(self, sample):
        super(sampleWidget, self).__init__()

        # ------------------------------- parameter initialization -------------------------------------#
        self.data_dict = {}
        self.sample = sample  # variable used to define sample info
        self.structTableInfo = []  # used to keep track of the table info instead of constantly switching
        self.parameterFit = []  # keeps track of parameters to fit
        self.varData = {ele: [[['', ''], ['', ''], ['', '']] for i in range(len(sample.structure))] for ele in
                        sample.myelements}  # element variation data [identifier, ratio, form factor]
        self.eShift = dict()  # keep track of the energy shift
        self.ffScale = dict()  # keeps track of the form factor scaling information
        self.currentVal = []  # keep track of the fitting parameter values
        self.change_eShift = True  # boolean used to determine if eShift is changed
        self.varTable = QTableWidget()  # element variation table
        self.elementBox = QComboBox()  # used to select which element to change element variation
        self.elementBox.setStyleSheet('background-color: white;')
        self.variationElements = sample.poly_elements  # retrieve the variation elements from the sample

        self.struct_ff = []  # list that contains the structural form factors
        self.mag_ff = []  # list that contains the magnetic form factors

        self.magData = {ele: [[[''], [''], ['']] for i in range(len(sample.structure))] for ele in
                        sample.myelements}  # initialize the magnetic data [identifier, density, form factor]
        self.magGo = True  # boolean currently not in use
        self.magDirection = ['z' for i in range(len(sample.structure))] # initialize the magnetization direction
        self.magDirBox = QComboBox()  # magnetization direction selection (x,y,z direction)
        self.magDirBox.setStyleSheet('background-color: white;')
        self.getData()  # gets the element variation and magnetic information

        self.magTable = QTableWidget()  # magnetization property table

        self.resetX = False  # boolean used to determine if fitting parameters are reset

        self.change_elements = False  # boolean used to determine if elements are currently being changed
        self.element_index = 0  # keeps track of which positional element is being used (A-site, B-site, or O-site)

        self.previousLayer = 0  # what was the previous layer
        self.changeLayer = False  # no longer in use
        self.firstStruct = True # no longer in use

        # ------------------------------- Widget Layout -------------------------------------#
        # setting up step size
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

        # initializing layerList based on number of layers
        layerList = []
        for i in range(len(self.sample.structure)):
            if i == 0:
                layerList.append('Substrate')
            else:
                layerList.append('Layer ' + str(i))

        # change this for an arbitrary sample model
        self.layerBox.addItems(layerList)
        self.layerBox.currentIndexChanged.connect(self.setTable)
        # changes the table on the screen when new layer selected

        # buttons for adding and removing layers
        cblayout.addStretch(1)
        cblayout.addWidget(addlayerButton)
        cblayout.addWidget(copylayerButton)
        cblayout.addWidget(deletelayerButton)
        cblayout.addSpacing(50)
        cblayout.addLayout(step_size_layout)

        # layer combo box
        cblayout.addWidget(self.layerBox)
        cblayout.addStretch(1)

        self.sampleInfoLayout = QStackedLayout()  # stacked layout for the different parameter types

        # setting up structural parameter table
        self.structTable = QTableWidget()
        self.structTable.setRowCount(3)
        self.structTable.setColumnCount(7)
        self.structTable.setHorizontalHeaderLabels(
            ['Element', 'Thickness (Å)', 'Density (mol/cm^3)', 'Roughness (Å)', 'Linked Roughness (Å)', 'Scattering Factor', 'Stoichiometry'])

        afont = QtGui.QFont()
        afont.setBold(True)
        self.structTable.horizontalHeader().setSectionResizeMode(1)
        self.structTable.horizontalHeader().setFont(afont)
        # <sup>2</sup>
        self._setStructFromSample(sample)  # setting the structural table

        # setTable
        self.setTable()

        # initializing energy shift table (include ff scaling later)
        self.energyShiftTable = QTableWidget()
        self.energyShiftTable.setColumnCount(2)
        #self.energyShiftTable.setHorizontalHeaderLabels(['Energy Shift (eV)'])
        self.energyShiftTable.setVerticalHeaderLabels(['Energy Shift (eV)', 'Scale'])

        self.energyShiftTable.verticalHeader().setFont(afont)
        self.energyShiftTable.horizontalHeader().setFont(afont)

        # variationWidget and magneticWidget initialization
        self.elementVariation = variationWidget(self, self.sample)
        self.elementMagnetic = magneticWidget(self, self.sample)

        # adding widgets to the stacked layout
        self.sampleInfoLayout.addWidget(self.structTable)
        self.sampleInfoLayout.addWidget(self.elementVariation)
        self.sampleInfoLayout.addWidget(self.elementMagnetic)
        self.sampleInfoLayout.addWidget(self.energyShiftTable)

        # setting up data fitting implementation
        self.structTable.viewport().installEventFilter(self)
        self.varTable.viewport().installEventFilter(self)
        self.magTable.viewport().installEventFilter(self)
        self.energyShiftTable.viewport().installEventFilter(self)

        selectlayout = QVBoxLayout()
        selectlayout.addStretch(1)

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

        self.shiftButton = QPushButton('Form Factor')  # energy shift button
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.shiftButton.clicked.connect(self._energy_shift)
        selectlayout.addWidget(self.shiftButton)

        # determine which widget is currently in use
        self.structBool = True
        self.polyBool = False
        self.magBool = False

        selectlayout.addSpacing(50)

        # button used to plot the density profile
        dpButton = QPushButton('Density Profile')
        dpButton.clicked.connect(self._densityprofile)
        dpButton.setStyleSheet("background-color : cyan")
        selectlayout.addWidget(dpButton)
        selectlayout.addStretch(1)

        # set the page layout
        pagelayout.addLayout(cblayout)
        pagelayout.addLayout(self.sampleInfoLayout)
        pagelayout.addLayout(selectlayout)

        mylayout = QVBoxLayout()  # layout that includes the plotting layout and the graph
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
        self.structTable.setItemDelegateForColumn(0, delegate)
        self.setLayout(mylayout)

    def check_element_number(self):
        """
        Purpose: Checks to make sure that there is a sample defined and the number of elements in each layer is the same
        :return: booleans not_empty and equal_elements where True states there are no defined layers and there
                 are not an equal number of elements in each layer, respectively.
        """

        not_empty = True
        equal_elements = True
        n = len(self.structTableInfo)
        if n == 0:  # no defined layers
            not_empty = False
        else:
            substrate = len(self.structTableInfo[0])
            for i in range(n):
                if len(self.structTableInfo[0]) != substrate: # unequal number of elements throughout the sample layers
                    equal_elements = False

        return not_empty, equal_elements

    def _energy_shift(self):
        """
        Purpose: Change to energy shift info
        :return:
        """
        self.sampleInfoLayout.setCurrentIndex(3)
        self.structButton.setStyleSheet('background: lightGrey')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.magButton.setStyleSheet('background: lightGrey')
        self.shiftButton.setStyleSheet('background: blue; color: white')
        self.setTableEShift()

    def setTableEShift(self):
        """
        Purpose: Reset energy table
        """

        self.energyShiftTable.blockSignals(True)  # block energy shift signals
        keys = list(self.eShift.keys())  # retrieve energy shift keys

        # set the energy shift form factor names and their values
        self.energyShiftTable.setColumnCount(len(keys))
        self.energyShiftTable.setRowCount(2)
        self.energyShiftTable.setHorizontalHeaderLabels(keys)
        self.energyShiftTable.setVerticalHeaderLabels(['E (eV)','Scale'])
        for column, key in enumerate(keys):  # setting energy shift
            item = QTableWidgetItem(str(self.eShift[key]))
            self.energyShiftTable.setItem(0, column, item)
        for column, key in enumerate(keys):  # setting scaling factor
            item = QTableWidgetItem(str(self.ffScale[key]))
            self.energyShiftTable.setItem(1, column, item)

        # set the parameter fit colors
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
            self.energyShiftTable.item(0, col).setBackground(QtGui.QColor(0, 255, 0))

        # checks to see if user has any form factors in the current working directory
        # - sets form factors in directory to blue
        for col, key in enumerate(keys):
            if key.startswith('ff-'):
                name = key.strip('ff-')
                if name in self.struct_ff:
                    self.energyShiftTable.horizontalHeaderItem(col).setForeground(QtGui.QColor(0, 0, 255))

            elif key.startswith('ffm'):
                name = key.strip('ffm-')
                if name in self.mag_ff:
                    self.energyShiftTable.horizontalHeaderItem(col).setForeground(QtGui.QColor(0, 0, 255))

        self.energyShiftTable.blockSignals(False)
    def eShiftFromSample(self, sample):
        """
        Purpose: Sets the energy shift from the input sample. This function is mostly used for the loading functions.
        :param sample: slab class
        :return:
        """
        # loop through structural form factors
        for key in list(sample.eShift.keys()):
            my_key = 'ff-' + key
            self.eShift[my_key] = sample.eShift[key]

        # loop through magnetic form factors
        for key in list(sample.mag_eShift.keys()):
            my_key = 'ffm-' + key
            self.eShift[my_key] = sample.mag_eShift[key]

    def changeEShiftValues(self):
        """
        Purpose: Change the energy shift values when signaled
        """
        column = self.energyShiftTable.currentColumn()  # retrieves current column
        row = self.energyShiftTable.currentRow()  # retrieves current row
        key = self.energyShiftTable.horizontalHeaderItem(column).text()  # retrieves the key
        value = self.energyShiftTable.item(row, column).text()  # retrieves the value

        if row == 0:  # energy shift case
            self.eShift[key] = float(value)  # change the value in memory
        elif row == 1: # ff scaling
            self.ffScale[key] = float(value)

        if row == 0:
            # make changes to the sample
            if key[:3] == 'ff-':
                sf = key[3:]
                self.sample.eShift[sf] = float(value)
            elif key[:3] == 'ffm':
                sf = key[4:]
                self.sample.mag_eShift[sf] = float(value)

        elif row == 1:  # need to check if form factor name changed!!
            # make changes to the sample
            if key[:3] == 'ff-':
                sf = key[3:]
                self.sample.ff_scale[sf] = float(value)
            elif key[:3] == 'ffm':
                sf = key[4:]
                self.sample.ffm_scale[sf] = float(value)

        # update changes to the fitting parameter values and their boundaries
        # ff scaling is not included in fitting
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
        """
        Purpose: Change the magnetization values when a change signaled
        """

        layer = self.layerBox.currentIndex()  # retrieve current layer

        column = self.magTable.currentColumn()  # retrieve current column
        row = self.magTable.currentRow()  # retrieve current row

        # make the changes
        if self.magTable.item(row, column) is not None and self.magTable.verticalHeaderItem(row) is not None:
            name = self.magTable.verticalHeaderItem(row).text()  # retrieves name
            element = ''  # pre-initializes element
            idx = 0  # pre-initializes index

            # Determines if the current magnetic element already exists and to what element it belongs to
            for i in range(len(self.structTableInfo[layer])):
                ele = self.structTableInfo[layer][i][0]
                if name in self.magData[ele][layer][0]:
                    idx = list(self.magData[ele][layer][0]).index(name)
                    element = copy.copy(ele)

            value = self.magTable.item(row, column).text()  # retrieves current value

            prev_value = self.magData[element][layer][column + 1][idx]  # retrieves previous value
            if column == 1:  # form factor column
                prev_dict_name = 'ffm-' + prev_value
                if prev_value != '':  # removes the previous form factor from the energy shift
                    del self.eShift[prev_dict_name]
                    del self.ffScale[prev_dict_name]

                if value != '':  # sets eShift to 0
                    dict_name = 'ffm-' + value
                    self.eShift[dict_name] = 0
                    self.ffScale[dict_name] = 1

                if prev_value != value:  # case where the name of the form factor has changed from previous value
                    # replaces form factor in every layer which the old form factors appeared
                    for i in range(len(self.magData[element])):
                        inLayer = False
                        for j in range(len(self.structTableInfo[i])):
                            if element == self.structTableInfo[i][j][0]:
                                inLayer = True
                        if inLayer:
                            self.magData[element][i][column + 1][idx] = value
                            if self.magData[element][i][1][idx] == '':
                                self.magData[element][i][1][idx] = '0'
                            if column == 1:
                                if prev_value != '' and value == '':
                                    self.magData[element][i][1][idx] = ''

            self.magData[element][layer][column + 1][idx] = value

            # updates fitting parameters and their boundaries
            copy_of_list = copy.deepcopy(self.parameterFit)
            for fit in copy_of_list:
                if column == 0:  # density
                    if layer == fit[0] and fit[1] == 'MAGNETIC' and fit[-1] == name:
                        idx = self.parameterFit.index(fit)
                        lower = float(value) - 0.01
                        if lower < 0:
                            lower = 0

                        upper = str(float(value) + 0.01)
                        self.currentVal[idx] = [value, [str(lower), upper]]
                elif column == 1 and fit[0] == 'SCATTERING FACTOR' and fit[
                    1] == 'MAGNETIC':  # magnetic scattering factor
                    if value != prev_value and prev_value == fit[2]:
                        self.parameterFit.remove(fit)
                        self.magTable.item(row, column).setBackground(QtGui.QColor(255, 255, 255))

        self.magTable.verticalHeader().setSectionResizeMode(1)

    def changeVarValues(self):
        """
        Purpose: change the variation parameters when signaled
        """

        layer = self.layerBox.currentIndex()  # retrieve current layer
        ele_idx = self.elementBox.currentIndex() # retrieve elementBox index
        element = self.structTableInfo[layer][ele_idx][0] # retrieve current element
        column = self.varTable.currentColumn() # retrieve current column
        row = self.varTable.currentRow() # retrieve current row

        # change the element variation info
        if self.varTable.item(row, column) is not None and not (self.change_elements):

            copy_of_list = copy.deepcopy(self.parameterFit)

            value = self.varTable.item(row, column).text()  # setting varData correctly depending on user input
            prev_value = self.varData[element][layer][column][row]

            if column == 0:  # checking the name
                empty = True  # has the user input an identifier
                my_num = 0  # counter that keeps track if we need to initialize the data
                for lay in range(len(self.varData[element])):
                    for i in range(len(self.varData[element][lay][0])):  # only need to check the names
                        if self.varData[element][lay][0][i] != '':
                            my_num = my_num + 1
                            empty = False

                # no value provided by the user
                if value != '':
                    my_num = my_num + 1
                    empty = False

                # varData needs to be intialized for entered identifier
                if my_num == 1:
                    for lay in range(len(self.varData[element])):
                        c = len(self.varData[element][lay][0])
                        if empty:
                            self.magData[element][lay] = [[element], [''], ['']]
                        else:
                            self.magData[element][lay] = [['' for i in range(c)],['' for i in range(c)],['' for i in range(c)]]


                # previous value and current value not the same
                if prev_value != value or prev_value == '':
                    # checks each layer and changes the previous identifier to the new one
                    for i in range(len(self.varData[element])):
                        inLayer = False

                        for j in range(len(self.structTableInfo[i])):
                            if element == self.structTableInfo[i][j][0]:
                                inLayer = True

                        if inLayer:
                            self.varData[element][i][0][row] = value
                            self.magData[element][i][0][row] = value

                            # preset all ratio values as a zero if not already initialized
                            if self.varData[element][i][1][row] == '':
                                self.varData[element][i][1][row] = 0


                # Now we check if the user has erased the element variation names
                isReset = True
                for e in self.varData[element][layer][0]:
                    if e != '':  # user has not reset all the names
                        isReset = False

                if isReset:  # reset all the values
                    # obtain the old for factors
                    my_ff = copy.copy(self.varData[element][layer][2])
                    my_ffm = copy.copy(self.magData[element][layer][2])
                    for lay in range(len(self.varData[element])):
                        c = len(self.varData[element][lay][0])
                        # make necessary changes to the element variation
                        self.varData[element][lay][1] = ['' for i in range(c)]
                        self.varData[element][lay][2] = ['' for i in range(c)]

                        # reset all magnetic data related to the element variation
                        self.magData[element][lay] = [[element], [''], ['']]


                    # delete the form factor in eShift
                    for ff_key in my_ff:
                        if ff_key != '':
                            ff_name = 'ff-'+ff_key
                            if ff_name in list(self.eShift.keys()):
                                del self.eShift[ff_name]
                                del self.ffScale[ff_name]

                    # delete the form factor in eShift
                    for ffm_key in my_ffm:
                        if ffm_key != '':
                            ffm_name = 'ffm-' + ffm_key
                            if ffm_name in list(self.eShift.keys()):
                                del self.eShift[ffm_name]
                                del self.ffScale[ffm_name]



            if column == 1:  # ratio
                # only sets the value if varData initialized
                if type(self.varData[element][layer][column]) is np.ndarray:
                    self.varData[element][layer][column] = list(self.varData[element][layer][column])
                self.varData[element][layer][column][row] = value

            # changing the scattering factor
            if column == 2:
                if prev_value != value:

                    # takes into account the scattering factor change
                    prev_dict_name = 'ff-' + prev_value
                    if prev_value != '':
                        del self.eShift[prev_dict_name]
                        del self.ffScale[prev_dict_name]

                    # does not include '' as a possible eShift value
                    if value != '':
                        dict_name = 'ff-' + value
                        self.eShift[dict_name] = 0
                        self.ffScale[dict_name] = 1

                    for idx, my_layer in enumerate(self.varData[element]):

                        if prev_value in my_layer[2]:
                            ZQ = list(my_layer[2]).index(prev_value)
                            self.varData[element][idx][2][ZQ] = value

            # updating the fitting parameters and their boundaries
            for fit in copy_of_list:
                if fit[0] == layer and fit[1] == 'POLYMORPHOUS' and column != 2:
                    idx = self.parameterFit.index(fit)
                    if column == 0:  # case where we changed the element variation name
                        if prev_value != value and fit[3] == prev_value:
                            self.parameterFit[idx][3] = value
                    elif column == 1:  # case for change in ratio
                        if fit[2] == element and fit[3] == self.varData[element][layer][0][row]:
                            lower = float(value) - 0.2
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

        self.setTableVar()  # reset the variation table

    def changeStructValues(self):
        """
        Purpose: Change the structural parameters when signaled
        """
        layer = self.layerBox.currentIndex()  # retrieve the current layer
        row = self.structTable.currentRow()  # retrieve the current row
        column = self.structTable.currentColumn()  # retrieve the current column
        if self.structTable.item(row, column) is not None:
            copy_of_list = self.parameterFit
            value = self.structTable.item(row, column).text()  # retrieves the current value
            prev_value = copy.copy(self.structTableInfo[layer][row][column])  # retrieves the previous value

            # change name or scattering factor if it is changed (eShift name change)
            if column == 5:
                name = self.structTable.item(row, 0).text()
                if value != prev_value:
                    # changing scattering factor info
                    prev_dict_name = 'ff-' + prev_value
                    del self.eShift[prev_dict_name]
                    del self.ffScale[prev_dict_name]
                    dict_name = 'ff-' + value
                    self.eShift[dict_name] = 0
                    self.ffScale[dict_name] = 1
                    for i in range(len(self.structTableInfo)):
                        for j in range(len(self.structTableInfo[i])):
                            if self.structTableInfo[i][j][0] == name:
                                self.structTableInfo[i][j][5] = value  # changes the scattering factor for all

            self.structTableInfo[layer][row][column] = self.structTable.item(row, column).text()  # update info

            # update fitting parameters and their info!
            for fit in copy_of_list:
                element = self.structTableInfo[layer][row][0]
                if fit[0] == layer and fit[1] == 'STRUCTURAL':
                    if column == 1:  # thickness
                        if fit[2] == 'COMPOUND' and fit[3] == 'THICKNESS':
                            if [layer, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'THICKNESS'])
                                lower = float(value) - 5
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 5)
                                self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == element and fit[4] == 'THICKNESS':
                            if [layer, 'STRUCTURAL', 'ELEMENT', element, 'THICKNESS'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'THICKNESS'])
                                lower = float(value) - 5
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 5)
                                self.currentVal[idx] = [value, [str(lower), upper]]

                    elif column == 2:  # density
                        if fit[2] == 'COMPOUND' and fit[3] == 'DENSITY':
                            if [layer, 'STRUCTURAL', 'COMPOUND', 'DENSITY'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'DENSITY'])
                                lower = float(value) - 0.01
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 0.01)
                                self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == element and fit[4] == 'DENSITY':
                            if [layer, 'STRUCTURAL', 'ELEMENT', element, 'DENSITY'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'DENSITY'])
                                lower = float(value) - 0.01
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 0.01)
                                self.currentVal[idx] = [value, [str(lower), upper]]
                    elif column == 3:  # roughness
                        if fit[2] == 'COMPOUND' and fit[3] == 'ROUGHNESS':
                            if [layer, 'STRUCTURAL', 'COMPOUND', 'ROUGHNESS'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'ROUGHNESS'])
                                lower = float(value) - 1
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 1)
                                self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == element and fit[3] == 'ROUGHNESS':
                            if [layer, 'STRUCTURAL', 'ELEMENT', element, 'ROUGHNESS'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'ELEMENT', element, 'ROUGHNESS'])
                                lower = float(value) - 1
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 1)
                                self.currentVal[idx] = [value, [str(lower), upper]]
                    elif column == 4:  # linked roughness
                        if fit[2] == 'COMPOUND' and fit[3] == 'LINKED ROUGHNESS':
                            if [layer, 'STRUCTURAL', 'COMPOUND', 'LINKED ROUGHNESS'] in self.parameterFit:
                                idx = self.parameterFit.index([layer, 'STRUCTURAL', 'COMPOUND', 'LINKED ROUGHNESS'])
                                lower = float(value) - 1
                                if lower < 0:
                                    lower = 0
                                upper = str(float(value) + 1)
                                self.currentVal[idx] = [value, [str(lower), upper]]
                        elif fit[2] == 'ELEMENT' and fit[3] == 'LINKED ROUGHNESS':
                            if [layer, 'STRUCTURAL', 'ELEMENT', element, 'LINKED ROUGHNESS'] in self.parameterFit:
                                idx = self.parameterFit.index(
                                    [layer, 'STRUCTURAL', 'ELEMENT', element, 'LINKED ROUGHNESS'])
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
        """
        Purpose: Determine which handler to use when user right clicks on tables
        :param source: the source of the signal
        :param event: the type of event
        """
        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.button() == Qt.RightButton:
                idx = self.sampleInfoLayout.currentIndex()
                if idx == 0:
                    self.structure_handler()  # fit structural parameters
                elif idx == 1:
                    self.var_handler()  # fit element variation parameters
                elif idx == 2:
                    self.mag_handler()  # fit magnetic parameters
                elif idx == 3:
                    self.eShift_handler()  # fit energy shift parameters
        return False

    def eShift_handler(self):
        """
        Purpose: Handler used in setting fitting state of energy shift
        """
        column = self.energyShiftTable.currentColumn()  # which column is to be fit
        row = self.energyShiftTable.currentRow()
        copy_fit_list = copy.deepcopy(self.parameterFit)
        name = self.energyShiftTable.horizontalHeaderItem(column).text()  # name of the header in the selected column
        val = self.energyShiftTable.currentItem().text()  # value of the current energy shift

        top_menu = QMenu()  # initializes the menu

        menu = top_menu.addMenu("Menu")

        _fit = menu.addAction('Fit')
        _remove_fit = menu.addAction('Remove Fit')

        my_items = self.energyShiftTable.selectedIndexes()
        Disable = False
        for i in my_items:
            if i.row() == 1:
                Disable = True


        if row == 1 or Disable:  # disable fitting capability for the form factor scaling (no longer implemented)
            _fit.setDisabled(True)
            _remove_fit.setDisabled(True)

        action = menu.exec_(QtGui.QCursor.pos())  # initialize action to be taken
        for i in my_items:
            column = i.column()
            name = self.energyShiftTable.horizontalHeaderItem(column).text()
            if action == _fit:
                # add energy shift to the fitting parameters
                self.resetX = True
                fit = []
                if name[:3] == 'ff-':
                    fit = ['SCATTERING FACTOR', 'STRUCTURAL', name[3:]]
                elif name[:3] == 'ffm':
                    fit = ['SCATTERING FACTOR', 'MAGNETIC', name[4:]]

                if fit != [] and fit not in self.parameterFit:  # set the boundaries
                    self.parameterFit.append(fit)
                    lower = str(float(val) - 0.5)
                    upper = str(float(val) + 0.5)
                    self.currentVal.append([val, [lower, upper]])

            elif action == _remove_fit:
                # remove the parameter from the fit (if found in parameterFit list)
                self.resetX = True
                fit = []
                if name[:3] == 'ff-':
                    fit = ['SCATTERING FACTOR', 'STRUCTURAL', name[3:]]
                elif name[:3] == 'ffm':
                    fit = ['SCATTERING FACTOR', 'MAGNETIC', name[4:]]

                # remove the fit
                if fit in self.parameterFit:
                    idx = self.parameterFit.index(fit)
                    self.parameterFit.remove(fit)
                    self.currentVal.pop(idx)

        # reset the tables
        self.setTableEShift()
        self.setTable()
        self.setTableMag()
        self.setTableVar()

    def mag_handler(self):
        """
        Purpose: Handler used in setting fitting state of magnetic parameters
        """
        idx = self.sampleInfoLayout.currentIndex()  # keeps track of which parameters are to be fit
        my_layer = self.layerBox.currentIndex()  # layer index
        copy_fit_list = copy.deepcopy(self.parameterFit)
        row = self.magTable.currentRow()  # current layer
        column = self.magTable.currentColumn()  # current column

        # find the element that the variation belongs too
        name = self.magTable.verticalHeaderItem(row).text()
        element = ''
        if name in list(self.magData.keys()):
            element = name
        else:
            for key in list(self.magData.keys()):  # loop through all the keys
                if name in self.magData[key][my_layer][0]:
                    element = key

        top_menu = QMenu()  # initializes the menu

        menu = top_menu.addMenu("Menu")

        _fit = menu.addAction('Fit')
        _remove_fit = menu.addAction('Remove Fit')
        my_items = self.magTable.selectedItems()

        Disable = False
        for i in my_items:
            item = self.magTable.item(i.row(),i.column()).text()
            if item == '':
                Disable = True


        if column == 2 or Disable:  # disable fitting capability for the form factor scaling (no longer implemented)
            _fit.setDisabled(True)
            _remove_fit.setDisabled(True)

        action = menu.exec_(QtGui.QCursor.pos())  # sets action

        for i in my_items:
            row = i.row()
            column = i.column()

            name = self.magTable.verticalHeaderItem(row).text()
            element = ''
            if name in list(self.magData.keys()):
                element = name
            else:
                for key in list(self.magData.keys()):  # loop through all the keys
                    if name in self.magData[key][my_layer][0]:
                        element = key
            if action == _fit:
                self.resetX = True
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
                        scattering_factor = self.magTable.item(row, 1)
                        if fit[0] == 'SCATTERING FACTOR' and fit[1] == 'MAGNETIC' and fit[2] == scattering_factor:
                            alreadySelected = True

                # Check to make sure that parameter not already selected
                if not alreadySelected:
                    if column == 1:
                        scattering_factor = self.magTable.item(row, column).text()
                        self.currentVal.append([0, [-0.5, 0.5]])
                        self.parameterFit.append(['SCATTERING FACTOR', 'MAGNETIC', scattering_factor])

                    elif column == 0:
                        value = self.magTable.item(row, column).text()
                        lower = float(value) - 0.01
                        if lower < 0:
                            lower = 0
                        upper = float(value) + 0.01
                        self.currentVal.append([float(value), [lower, upper]])
                        if name == element:
                            self.parameterFit.append([my_layer, 'MAGNETIC', element])
                        else:
                            self.parameterFit.append([my_layer, 'MAGNETIC', element, name])

            if action == _remove_fit:
                # removes fit from parameterFit
                self.resetX = True
                for fit in copy_fit_list:
                    if column == 0:  # magnetic density
                        if len(fit) == 3 and fit[1] == 'MAGNETIC':  # not polymorphous case
                            if fit[0] == my_layer and fit[2] == name:
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                        elif len(fit) == 4 and fit[1] == 'MAGNETIC':  # polymorphous case
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
        """
        Purpose: Handler used in setting fitting state of element variation
        """
        idx = self.sampleInfoLayout.currentIndex()  # keeps track of which parameters are to be fit
        my_layer = self.layerBox.currentIndex()  # layer index
        my_ele = self.elementBox.currentIndex()  # element index
        copy_fit_list = copy.deepcopy(self.parameterFit)  # fit list

        element = self.structTableInfo[my_layer][my_ele][0]  # retrieves element in layer selected

        top_menu = QMenu()

        menu = top_menu.addMenu("Menu")  # initializes menu

        _fit = menu.addAction('Fit')
        _remove_fit = menu.addAction('Remove Fit')

        row = self.varTable.currentRow()  # retrieves current row
        column = self.varTable.currentColumn()  # retrieves current column

        Disable = False
        my_items = self.varTable.selectedIndexes()
        my_columns = []
        for i in my_items:
            col = i.column()
            if col in my_columns:
                Disable = True
            my_columns.append(row)
            for fit in self.parameterFit:
                if column == 1:
                    if fit[1] == 'POLYMORPHOUS' and fit[0] == my_layer:
                        name = self.varTable.item(row, 0).text()
                        if element == fit[2] and name != fit[3]:
                            Disable = True

        value = self.varData[element][my_layer][1][row]

        if value == 1 and len(self.varData[element][my_layer][1]) > 2:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Fitting Error")
            messageBox.setText('Fitting and element variation with ratio value of 1 will result in a divide by zero. Please select another control variable to use. ')
            messageBox.exec()
            Disable = True

        if column == 0 or Disable:  # disable fitting for identifier (element variation name)
            _fit.setDisabled(True)
            _remove_fit.setDisabled(True)

        if value == 1 and len(self.varData[element][my_layer][1]) > 2:
            Disable = True
        action = menu.exec_(QtGui.QCursor.pos())
        for i in my_items:
            row = i.row()
            column = i.column()
            if action == _fit:
                # add parameter to parameterFit
                self.resetX = True
                alreadySelected = False
                for fit in copy_fit_list:
                    if column == 1:
                        if fit[1] == 'POLYMORPHOUS' and fit[0] == my_layer:
                            name = self.varTable.item(row, 0).text()
                            if element == fit[2] and name == fit[3]:
                                alreadySelected = True

                    elif column == 2:
                        scattering_factor = self.varTable.item(row, column).text()
                        if len(fit) == 3 and scattering_factor == fit[2]:
                            alreadySelected = True

                if not alreadySelected:  # case where parameter has not been selected yet
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

                    elif column == 2:  # fitting the energy shift for the structural scattering factor
                        scattering_factor = self.varTable.item(row, 2).text()
                        if scattering_factor != '':
                            self.parameterFit.append(['SCATTERING FACTOR', 'STRUCTURAL', scattering_factor])
                            self.currentVal.append([0, [-0.5, 0.5]])

            elif action == _remove_fit:
                # remove the fit
                self.resetX = True
                for fit in copy_fit_list:
                    if column == 1:
                        if len(fit) == 4 and fit[1] == 'POLYMORPHOUS':
                            name = self.varTable.item(row, 0).text()
                            if my_layer == fit[0] and element == fit[2] and name == fit[3]:
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                    elif column == 2:
                        if len(fit) == 3:
                            scattering_factor = self.varTable.item(row, 2).text()
                            if scattering_factor == fit[2] and fit[1] == 'STRUCTURAL':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
        self.setTableVar()

    def structure_handler(self):
        """
        Purpose: Handler used in setting fitting state of energy shift
        """

        idx = self.sampleInfoLayout.currentIndex()  # current index
        my_layer = self.layerBox.currentIndex()  # current layer
        copy_fit_list = copy.deepcopy(self.parameterFit)
        top_menu = QMenu()

        menu = top_menu.addMenu("Menu")  # initializing the menu

        # initializing options
        _element_fit = menu.addAction('Element Fit')
        _compound_fit = menu.addAction('Compound Fit')
        _remove_fit = menu.addAction('Remove Fit')

        for i in self.structTable.selectedIndexes():
            row = i.row()  # current row
            column = i.column()  # current column

            # only allows user to fit certain parameters
            if column == 5:  # disable compound fit if scattering factor selected
                _compound_fit.setDisabled(True)
            elif column == 0 or (column == 1 and my_layer == 0): # disable fitting for substrate thickness
                _element_fit.setDisabled(True)
                _compound_fit.setDisabled(True)
                _remove_fit.setDisabled(True)
            elif column == 4:  # disable fit for linked roughness
                LRough = self.structTableInfo[my_layer][row][4]
                if str(LRough).upper() == 'FALSE' or LRough == '':  # there exists a False or '' as input for linked rough
                    _element_fit.setDisabled(True)
                    _compound_fit.setDisabled(True)
                    _remove_fit.setDisabled(True)
                compound = True
                for cool in self.structTableInfo[my_layer]:  # disables coumpound fit if required
                    if str(cool[4]).upper() == 'FALSE' or cool[4] == '':
                        compound = False

                if not compound:
                    _compound_fit.setDisabled(True)
            elif column == 6:  # disables fitting for stoichiometry
                _element_fit.setDisabled(True)
                _compound_fit.setDisabled(True)
                _remove_fit.setDisabled(True)

        action = menu.exec_(QtGui.QCursor.pos())
        my_indices = []
        for i in self.structTable.selectedIndexes():
            row = i.row()  # current row
            column = i.column()  # current column
            # Element Mode
            if action == _element_fit:
                self.resetX = True
                value = self.structTable.item(row,column).text()
                element = self.structTable.item(row, 0).text()
                alreadySelected = False

                # Check to make sure parameter is not already selected
                for fit in copy_fit_list:
                    # check if layer and parameter in compound mode
                    n = len(fit)
                    if n == 3:  # scattering factor
                        if column == 5:  # trying to fit scattering factor
                            item = self.structTable.item(row,column).text()
                            if item == fit[2]:
                                alreadySelected = True

                    elif n == 5:  # element mode
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

                        if layer == my_layer and column == param_num and column not in my_indices:
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)


                    elif n == 5:  # compound mode
                        layer = fit[0]
                        ele = fit[4]
                        param = fit[3]
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


                        if layer == my_layer and column == param_num and my_ele == ele:
                            alreadySelected = True


                if not alreadySelected:  # add parameter to parameterFit and pre-set boundary
                    if column == 1:  # thickness
                        if [my_layer, 'STRUCTURAL', 'ELEMENT', element, 'THICKNESS'] not in self.parameterFit:
                            self.parameterFit.append([my_layer, 'STRUCTURAL', 'ELEMENT', element, 'THICKNESS'])
                            lower = float(value) - 5
                            if lower < 0:
                                lower = 0
                            upper = float(value) + 5
                            self.currentVal.append([float(value), [lower, upper]])
                    elif column == 2:  # density
                        if [my_layer, 'STRUCTURAL', 'ELEMENT', element, 'DENSITY'] not in self.parameterFit:
                            self.parameterFit.append([my_layer, 'STRUCTURAL', 'ELEMENT', element, 'DENSITY'])
                            lower = float(value) - 0.01
                            if lower < 0:
                                lower = 0
                            upper = float(value) + 0.01
                            self.currentVal.append([float(value), [lower, upper]])
                    elif column == 3:  # roughness
                        if [my_layer, 'STRUCTURAL', 'ELEMENT', element, 'ROUGHNESS'] not in self.parameterFit:
                            self.parameterFit.append([my_layer, 'STRUCTURAL', 'ELEMENT', element, 'ROUGHNESS'])
                            lower = float(value) - 1
                            if lower < 0:
                                lower = 0
                            upper = float(value) + 1
                            self.currentVal.append([float(value), [lower, upper]])
                    elif column == 4:  # linked roughness
                        if [my_layer, 'STRUCTURAL', 'ELEMENT', element, 'LINKED ROUGHNESS'] not in self.parameterFit:
                            self.parameterFit.append([my_layer, 'STRUCTURAL', 'ELEMENT', element, 'LINKED ROUGHNESS'])
                            lower = float(value) - 1
                            if lower < 0:
                                lower = 0
                            upper = float(value) + 1
                            self.currentVal.append([float(value), [lower, upper]])
                    elif column == 5:  # scattering factor
                        scattering_factor = self.structTable.item(row, 5).text()
                        if ['SCATTERING FACTOR', 'STRUCTURAL', scattering_factor] not in self.parameterFit:

                            if scattering_factor[0] != '[':
                                self.parameterFit.append(['SCATTERING FACTOR', 'STRUCTURAL', scattering_factor])
                                self.currentVal.append([0, [-0.5, 0.5]])

            elif action == _compound_fit:
                self.resetX = True
                # retrieve minimum value in the row
                my_vals = list()
                for i in range(self.structTable.rowCount()):
                    my_vals.append(float(self.structTable.item(i, column).text()))

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

                        if param_n == column and layer == my_layer and column not in my_indices:
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)

                if not alreadySelected:  # adds parameter to parameterFit and sets boundary
                    if column == 1 and column not in my_indices:  # thickness
                        if my_layer != 0:
                            my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'THICKNESS', my_row]
                            self.parameterFit.append(my_fit)
                            lower = float(value) - 5
                            if lower < 0:
                                lower = 0
                            upper = float(value) + 5
                            self.currentVal.append([float(value), [lower, upper]])
                    elif column == 2 and column not in my_indices:  # density
                        my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'DENSITY', my_row]
                        self.parameterFit.append(my_fit)
                        lower = float(value) - 0.01
                        if lower < 0:
                            lower = 0
                        upper = float(value) + 0.01
                        self.currentVal.append([float(value), [lower, upper]])
                    elif column == 3 and column not in my_indices:  # roughness
                        my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'ROUGHNESS', my_row]
                        self.parameterFit.append(my_fit)
                        lower = float(value) - 1
                        if lower < 0:
                            lower = 0
                        upper = float(value) + 1
                        self.currentVal.append([float(value), [lower, upper]])
                    elif column == 4 and column not in my_indices:  # linked roughness
                        my_fit = [my_layer, 'STRUCTURAL', 'COMPOUND', 'LINKED ROUGHNESS', my_row]
                        self.parameterFit.append(my_fit)
                        lower = float(value) - 1
                        if lower < 0:
                            lower = 0
                        upper = float(value) + 1
                        self.currentVal.append([float(value), [lower, upper]])

            elif action == _remove_fit:
                # removes the parameter from parameterFit
                self.resetX = True
                element = self.structTableInfo[my_layer][row][0]
                scattering_factor = self.structTableInfo[my_layer][row][5]
                for fit in copy_fit_list:
                    n = len(fit)
                    if column == 1:
                        mode = fit[2]
                        if mode == "ELEMENT" and my_layer == fit[0]:
                            ele = fit[3]
                            if ele == element and fit[4] == 'THICKNESS':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)


                        elif mode == 'COMPOUND' and my_layer == fit[0] and column not in my_indices:
                            if fit[3] == 'THICKNESS':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)

                    elif column == 2:
                        mode = fit[2]
                        if mode == "ELEMENT" and my_layer == fit[0]:
                            ele = fit[3]
                            if ele == element and fit[4] == 'DENSITY':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                        elif mode == 'COMPOUND' and my_layer == fit[0] and column not in my_indices:
                            if fit[3] == 'DENSITY':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                    elif column == 3:
                        mode = fit[2]
                        if mode == "ELEMENT" and my_layer == fit[0]:
                            ele = fit[3]
                            if ele == element and fit[4] == 'ROUGHNESS':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                        elif mode == 'COMPOUND' and my_layer == fit[0] and column not in my_indices:
                            if fit[3] == 'ROUGHNESS':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                    elif column == 4:
                        mode = fit[2]
                        if mode == "ELEMENT" and my_layer == fit[0]:
                            ele = fit[3]
                            if ele == element and fit[4] == 'LINKED ROUGHNESS':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                        elif mode == 'COMPOUND' and my_layer == fit[0] and column not in my_indices:
                            if fit[3] == 'LINKED ROUGHNESS':
                                idx = self.parameterFit.index(fit)
                                self.parameterFit.remove(fit)
                                self.currentVal.pop(idx)
                    elif column == 5 and n == 3:
                        if scattering_factor == fit[2] and fit[1] == 'STRUCTURAL':
                            idx = self.parameterFit.index(fit)
                            self.parameterFit.remove(fit)
                            self.currentVal.pop(idx)

            my_indices.append(column)
        self.setTable()  # reset the table

    def changeStepSize(self):
        """
        Purpose: Change step size when signaled
        """
        self._step_size = self.step_size.text()

    def _setVarMagFromSample(self, sample):
        """
        Purpose: Retrieve element variation and magnetic data from slab class
        :param sample: slab class
        """

        # setting up the varData and magData dictionaries
        num_layers = len(sample.structure)
        self.magDirection = ['z' for i in range(num_layers)]
        self.varData = dict()  # resets varData
        self.magData = dict()  # resets magData
        # need to add sample.myelements in the data file
        for ele in sample.myelements:
            self.varData[ele] = [[['', ''], ['', ''], ['', '']] for i in range(num_layers)]
            self.magData[ele] = [[[ele], [''], ['']] for i in range(num_layers)]

        # sets element variation data for all layers
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
                    n = len(element.mag_density)
                    mag_elements = sample.mag_elements[ele]

                    for i in range(n):
                        if self.magData[ele][idx][0][i] in mag_elements:
                            self.magData[ele][idx][1][i] = element.mag_density[i]
                            self.magData[ele][idx][2][i] = element.mag_scattering_factor[i]
                        else:
                            self.magData[ele][idx][1][i] = ''
                            self.magData[ele][idx][2][i] = ''





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
        """
        Purpose: Retrieves structural data from slab class
        :param sample: slab class
        """
        self.sample = sample  # resets internal slab class
        find_sf = self.sample.find_sf  # retrieves form factors
        self.structTableInfo = []  # resets structure info

        # loops through each layer and retrieves info
        for idx in range(len(sample.structure)):
            structInfo = sample.structure
            rows = len(structInfo[idx])  # determines the number of layers
            tempArray = np.empty((rows, 7), dtype=object)
            for row in range(rows):
                for col in range(7):
                    if col == 0:  # name
                        element = list(structInfo[idx].keys())[row]
                        tempArray[row, col] = str(element)
                    elif col == 1:  # thickness
                        element = list(structInfo[idx].keys())[row]
                        thickness = structInfo[idx][element].thickness
                        tempArray[row, col] = str(thickness)
                    elif col == 2:  # density
                        element = list(structInfo[idx].keys())[row]
                        density = structInfo[idx][element].density
                        tempArray[row, col] = str(density)
                    elif col == 3:  # roughness
                        element = list(structInfo[idx].keys())[row]
                        roughness = structInfo[idx][element].roughness
                        tempArray[row, col] = str(roughness)
                    elif col == 4:  # linked roughness
                        element = list(structInfo[idx].keys())[row]
                        linked_roughness = structInfo[idx][element].linked_roughness
                        tempArray[row, col] = str(linked_roughness)
                    elif col == 5:  # scattering factor
                        element = list(structInfo[idx].keys())[row]
                        if element in list(find_sf[0].keys()):
                            scattering_factor = find_sf[0][element]
                        else:
                            scattering_factor = structInfo[idx][element].scattering_factor

                        # keeps track of the scattering factors - implementation will change when loading data
                        if type(scattering_factor) is not list and type(scattering_factor) is not np.ndarray:
                            name = 'ff-' + scattering_factor
                            if name not in list(self.eShift.keys()):
                                self.eShift[name] = 0
                                self.ffScale[name] = 1

                        tempArray[row, col] = str(scattering_factor)
                    elif col == 6:  # stoichiometry
                        element = list(structInfo[idx].keys())[row]
                        stoichiometry = structInfo[idx][element].stoichiometry
                        tempArray[row, col] = str(stoichiometry)

            self.structTableInfo.append(tempArray)

    def clearTable(self):
        """
        Purpose: Clear structure table
        """
        self.structTable.setRowCount(0)
        self.structTable.setColumnCount(7)

    def setTable(self):
        """
        Purpose: set the structure table
        """

        if len(self.structTableInfo) != 0:  # check to make sure there is info
            self.structTable.blockSignals(True)  # disable any signals associated with the table
            idx = self.layerBox.currentIndex()  # retrieve the current index
            tableInfo = self.structTableInfo[idx]  # retrieve table info
            num_rows = len(tableInfo)  # determine number of rows required
            self.structTable.setRowCount(num_rows)  # set correct number of rows
            self.structTable.setColumnCount(7)

            # loops through each parameter
            for col in range(7):
                for row in range(num_rows):
                    item = QTableWidgetItem(str(tableInfo[row][col]))  # converts item to QTableWidgetItem
                    self.structTable.setItem(row, col, item)  # sets item

                    # change color to grey for substrate
                    if col == 0:
                        self.structTable.item(row, col).setBackground(QtGui.QColor('lightGray'))
                    if col == 1 and idx == 0:
                        self.structTable.item(row, col).setBackground(QtGui.QColor('lightGray'))

            # sets green color for element fit and purple for compound fit
            for fit in self.parameterFit:
                layer = fit[0]
                n = len(fit)
                if layer == idx:  # not scattering factor parameters
                    mode = fit[2]
                    if mode == 'COMPOUND':  # compound mode
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
                            self.structTable.item(row, 5).setBackground(QtGui.QColor(150, 255, 150))

            self.structTable.blockSignals(False)

    def clearVarTable(self):
        """
        Purpose: clear element variation table
        """
        self.varTable.setRowCount(0)
        self.varTable.setColumnCount(3)

    def setTableVar(self):
        """
        Purpose: Set element variation table
        """
        if len(self.structTableInfo) != 0:  # checks if there is element variation info
            self.varTable.blockSignals(True)  # blocks varTable signals
            idx = self.sampleInfoLayout.currentIndex()  # gets current index
            layer_idx = self.layerBox.currentIndex()  # gets current layer index
            ele_idx = self.elementBox.currentIndex()  # gets current element index

            # makes sure that when we switch layers we show the same positional element
            if not self.change_elements:
                ele_idx = self.elementBox.currentIndex()
                self.element_index = copy.deepcopy(ele_idx)

            if ele_idx != -1:  # do nothing if element index is -1
                ele = self.structTableInfo[layer_idx][ele_idx][0]  # retrieves element symbol

                info = self.varData[ele][layer_idx]  # gets element variation info

                # sets the element variation info to the table
                if len(info[0]) != 0:
                    self.varTable.setRowCount(len(info[0]))

                    # loops through all the info
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

                        # sets colors for fit parameters
                        if idx == 1:
                            for fit in self.parameterFit:
                                if len(fit) == 3 and fit[2] == info[2][row] and fit[1] == 'STRUCTURAL':
                                    self.varTable.item(row, 2).setBackground(QtGui.QColor(150, 255, 150))
                                elif len(fit) == 4:
                                    if fit[0] == layer_idx and fit[1] == "POLYMORPHOUS" and fit[2] == ele and fit[3] == \
                                            info[0][row]:
                                        self.varTable.item(row, 1).setBackground(QtGui.QColor(150, 255, 150))
                else:  # no element variation info
                    for row in range(self.varTable.rowCount()):
                        item = QTableWidgetItem('')
                        self.varTable.setItem(row, 0, item)

                        item = QTableWidgetItem('')
                        self.varTable.setItem(row, 1, item)

                        item = QTableWidgetItem('')
                        self.varTable.setItem(row, 2, item)

            self.varTable.blockSignals(False)

    def clearMagTable(self):
        """
        Purpose: clear magnetization table
        """
        self.magTable.setRowCount(0)
        self.magTable.setColumnCount(2)

    def setTableMag(self):

        if len(self.structTableInfo) != 0:  # checks if there is magnetization info
            self.magTable.blockSignals(True)
            layer_idx = self.layerBox.currentIndex()

            # set the magnetic direction combobox to the correct magnetization direction
            mag_idx = 0
            dir = self.magDirection[layer_idx]
            if dir == 'x':
                mag_idx = 0
            elif dir == 'y':
                mag_idx = 1
            elif dir == 'z':
                mag_idx = 2

            self.magDirBox.setCurrentIndex(mag_idx)

            labels = []
            density = []
            sf = []

            # Loops through all of the elements
            for ele_idx in range(len(self.structTableInfo[layer_idx])):
                element = self.structTableInfo[layer_idx][ele_idx][0]

                names = self.magData[element][layer_idx][0]  # retrieves the identifiers
                D = self.magData[element][layer_idx][1]  # retrieves the density
                S = self.magData[element][layer_idx][2]  # retrieves the form factors


                num = len(names)
                if num != 0:  # element variation case
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

                else:  # non-element variation
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
                    self.magTable.setItem(row, 0, item)
                else:
                    item = QTableWidgetItem(str(mydensity))
                    self.magTable.setItem(row, 0, item)

                    for fit in self.parameterFit:  # checks to see if element variation in current layer are to be fitted
                        if len(fit) == 3 and fit[0] != 'SCATTERING FACTOR':
                            if layer_idx == fit[0] and fit[1] == 'MAGNETIC' and fit[2] == name:
                                self.magTable.item(row, 0).setBackground(QtGui.QColor(150, 255, 150))
                        elif len(fit) == 4 and fit[1] == 'MAGNETIC':
                            if layer_idx == fit[0] and fit[1] == 'MAGNETIC' and fit[3] == name:
                                self.magTable.item(row, 0).setBackground(QtGui.QColor(150, 255, 150))

                if mysf == '':  # checks to see if scattering factor is filled
                    item = QTableWidgetItem('')
                    self.magTable.setItem(row, 1, item)
                else:
                    item = QTableWidgetItem(mysf)
                    self.magTable.setItem(row, 1, item)
                    for fit in self.parameterFit:
                        if len(fit) == 3:
                            scattering_factor = self.magTable.item(row, 1).text()
                            if fit[0] == 'SCATTERING FACTOR' and fit[1] == 'MAGNETIC' and fit[2] == scattering_factor:
                                self.magTable.item(row, 1).setBackground(QtGui.QColor(150, 255, 150))


            self.magTable.blockSignals(False)

    def _addLayer(self):
        """
        Purpose: initialize compoundInput widget and add layer info
        """
        # start add layer application
        addLayerApp = compoundInput()
        addLayerApp.show()
        addLayerApp.exec_()
        userinput = addLayerApp.val
        addLayerApp.close()

        # check to see if we have a new element

        num_layers = len(self.structTableInfo)

        my_elements = []


        # initializing the element variation and magnetic data
        if len(userinput) != len(self.structTableInfo[0]):  # checks to make sure number of elements is consistent
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Invalid Entry")
            messageBox.setText('Each layer must have the same number of structural elements. If you would like to define a different number of structural elements in a layer then a dummy variable can be used (e.g. A, D, E, G, J, L, M, Q, R, T, X, Z). The corresponding form factors and atomic mass of these variables are 0, respectively.')
            messageBox.exec()
        elif len(userinput) != 0:  # checks if user exit the coumpoundInput widget
            for i in range(len(userinput)):
                # includes new scattering factors into the energy shift
                if userinput[i][5] not in self.sample.eShift.keys():
                    self.sample.eShift[userinput[i][5]] = 0
                    name = 'ff-' + userinput[i][5]
                    self.eShift[name] = 0
                    self.ffScale[name] = 1
                my_elements.append(userinput[i][0])
                if userinput[i][0] not in list(self.varData.keys()):
                    self.varData[userinput[i][0]] = [[['', ''], ['', ''], ['', '']] for j in range(num_layers)]
                    self.magData[userinput[i][0]] = [[[userinput[i][0]], [''], ['']] for j in range(num_layers)]

            self.parameterFit = []  # resets the parameter fit info
            self.currentVal = []  # resets fitting parameter values

            num = self.layerBox.count()  # retrieves current number of layers
            idx = self.layerBox.currentIndex()  # gets current layer
            if num == 0:  # resets everything if there are no layers currently defined
                self.varData = dict()

                self.structTableInfo.insert(0, userinput)
                self.layerBox.blockSignals(True)
                self.layerBox.addItem('Substrate')
                self.layerBox.blockSignals(False)
                for i in range(len(userinput)):
                    self.varData[userinput[i][0]] = [[['', ''], ['', ''], ['', '']]]
                    self.magData[userinput[i][0]] = [[[userinput[i][0]], [''], ['']]]
                self.setTable()
                self.elementVariation.changeElements()
            else:  # there exists other layers
                self.layerBox.blockSignals(True)  # blocks signals from layerBox
                self.layerBox.addItem('Layer ' + str(num))
                self.layerBox.blockSignals(False)
                self.structTableInfo.insert(idx + 1, userinput)  # insert the layer info into the correct position

                # update varData info depending if element already exists in sample
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

                                    data = [info[0], [1, 0], info[2]]
                                    isVar = True
                                else:
                                    data = [['', ''], ['', ''], ['', '']]

                            if not isMag:
                                if info_mag[1][0] != '' or info_mag[2][0] != '':
                                    n = len(info_mag[0])
                                    data_mag = [info_mag[0], [0 for i in range(n)], info_mag[2]]
                                    isMag = True

                                else:
                                    data_mag = [[key], [''], ['']]
                    else:
                        data = [['', ''], ['', ''], ['', '']]
                        data_mag = [key, [''], [''],]

                    self.varData[key].insert(idx + 1, data)
                    self.magData[key].insert(idx + 1, data_mag)
                    self.magDirection.insert(idx + 1, 'z')
                self.layerBox.setCurrentIndex(idx + 1)

    def _removeLayer(self):
        """
        Purpose: Remove selected layer
        """

        num = self.layerBox.count()  # get the number of layers in the material

        idx = self.layerBox.currentIndex()  # determine which layer has been selected to remove
        last_var = dict()
        last_mag = dict()

        self.parameterFit = []
        self.currentVal = []

        if num != 0:
            self.layerBox.removeItem(num - 1)

            last_struct = self.structTableInfo[idx]

            self.structTableInfo.pop(idx)  # removes the information about that layer

            # removes the element variation data
            for key in list(self.varData.keys()):
                last_var[key] = self.varData[key][idx]
                last_mag[key] = self.magData[key][idx]
                self.varData[key].pop(idx)
                self.magData[key].pop(idx)

            if num != 1:
                self.setTable()  # sets the table for layer that replaces the removed layer
            else:
                self.clearTable()
                self.clearVarTable()
                self.clearMagTable()

            ff_to_remove = []
            ffm_to_remove = []

            elements_removed = dict()
            for layer_info in last_struct:
                elements_removed[layer_info[0]] = False

            # find out if the element still exists
            for layer in self.structTableInfo:
                for element_info in layer:
                    my_element = element_info[0]

                    if my_element in list(elements_removed.keys()):
                        elements_removed[my_element] = True

            # determine the scattering factors that we need to remove
            for key in list(elements_removed.keys()):
                n_ff = len(ff_to_remove)
                if not (elements_removed[key]):  # element has been removed
                    del self.magData[key]
                    del self.varData[key]

                    for sf in last_var[key][2]:  # element variation case
                        if sf != '' and sf not in ff_to_remove:
                            ff_to_remove.append(sf)
                    if len(ff_to_remove) == n_ff:
                        for ele_info in last_struct:
                            if ele_info[0] == key:
                                ff_to_remove.append(ele_info[5])

                    for sfm in last_mag[key][2]:
                        if sfm != '' and sfm not in ffm_to_remove:
                            ffm_to_remove.append(sfm)

            for key in list(self.eShift.keys()):
                if key.startswith('ff-'):
                    to_remove = key.strip('ff-')
                    if to_remove in ff_to_remove:
                        del self.eShift[key]
                        del self.ffScale[key]


                elif key.startswith('ffm-'):
                    to_remove = key.strip('ffm-')
                    if to_remove in ffm_to_remove:
                        del self.eShift[key]
                        del self.ffScale[key]

    def _copyLayer(self):
        """
        Purpose: Copy the current selected layer and add above
        """

        # reset fitting parameters
        self.parameterFit = []
        self.currentVal = []

        num = self.layerBox.count()  # current number of layers

        # created proper label
        if num == 0:
            self.layerBox.addItem('Substrate')
        else:
            self.layerBox.addItem('Layer ' + str(num))

        idx = self.layerBox.currentIndex()

        newLayer = copy.deepcopy(self.structTableInfo[idx])  # retrieve new layer info
        self.structTableInfo.insert(idx + 1, newLayer)  # add info to the structural info

        # goes through all the keys and copy info for element variation and magnetization
        for key in list(self.varData.keys()):
            info = copy.deepcopy(self.varData[key][idx])
            info_mag = copy.deepcopy(self.magData[key][idx])

            self.varData[key].insert(idx + 1, info)
            self.magData[key].insert(idx + 1, info_mag)

        new_dir = copy.deepcopy(self.magDirection[idx])
        self.magDirection.insert(idx + 1, new_dir)

    def _structural(self):
        """
        Purpose: Sets workspace to sampleWidget
        """
        self.structButton.setStyleSheet('background: blue; color: white')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.magButton.setStyleSheet('background: lightGrey')
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.sampleInfoLayout.setCurrentIndex(0)

    def _elementVariation(self):
        """
        Purpose: Sets workspace to variationWidget
        """
        self.structButton.setStyleSheet('background: lightGrey')
        self.polyButton.setStyleSheet('background: blue; color: white')
        self.magButton.setStyleSheet('background: lightGrey')
        self.shiftButton.setStyleSheet('background: lightGrey')
        self.sampleInfoLayout.setCurrentIndex(1)

    def _magnetic(self):
        """
        Purpose: Sets workspace to magneticWidget
        """
        self.structButton.setStyleSheet('background: lightGrey')
        self.polyButton.setStyleSheet('background: lightGrey')
        self.magButton.setStyleSheet('background: blue; color: white')
        self.shiftButton.setStyleSheet('background: lightGrey')

        self.sampleInfoLayout.setCurrentIndex(2)

    def _densityprofile(self):

        """
        Purpose: Calculate the density profile from the sample info
        """

        step_size = float(self.step_size.text())  # retrieve step size
        self.sample = self._createSample()  # transform sample information into slab class

        thickness, density, density_magnetic = self.sample.density_profile(step=step_size)  # compute density profile

        self.densityWidget.clear()  # clear graph
        self._plotDensityProfile(thickness, density, density_magnetic)  # plot density profile

    def _plotDensityProfile(self, thickness, density, density_magnetic):
        """
        Purpose: plot the density profile
        :param thickness: thickness numpy array
        :param density: dictionary containing the density of each element as a function of thickness
        :param density_magnetic: dictionary containing the magnetic density of each element as a function of thickness
        """

        # determine total number of density profiles to plot
        num = len(density)
        num = num + len(density_magnetic)

        val = list(density.values())  # density profiles
        mag_val = list(density_magnetic.values())  # magnetic density profile

        # checks to make sure that we are not dealing with surface key
        check = []
        for key in list(density.keys()):
            if key[-1].isdigit():
                check.append(True)
            else:
                check.append(False)

        # plot density profile
        for idx in range(len(val)):
            if check[idx]:
                self.densityWidget.plot(thickness, val[idx], pen=pg.mkPen((idx, num), width=2),
                                        name=list(density.keys())[idx])
            else:
                self.densityWidget.plot(thickness, val[idx], pen=pg.mkPen((idx, num), width=2),
                                        name=list(density.keys())[idx])

        # plot magnetic density profile
        for idx in range(len(mag_val)):
            myname = 'Mag: ' + list(density_magnetic.keys())[idx]
            self.densityWidget.plot(thickness, -mag_val[idx],
                                    pen=pg.mkPen((num - idx, num), width=2, style=Qt.DashLine), name=myname)

        # labels
        self.densityWidget.setLabel('left', "Density (mol/cm^3)")
        self.densityWidget.setLabel('bottom', "Thickness (Å)")

    def _createSample(self):
        """
        Purpose: Takes information from tables and converts into slab class
        :return: Information in format of slab class
        """

        m = len(self.structTableInfo)  # determines how many layers in the sample
        sample = ms.slab(m)  # initializes the slab class

        # loops through each layer and sets the appropriate parameters
        for idx in range(m):
            formula = ''  # used to determine the chemical formula
            thickness = []  # thickness of the layer
            density = []  # density of the layer
            roughness = []  # roughness of the layer
            linked_roughness = []  # linked roughness of the layer
            scat_fact = []

            layer = self.structTableInfo[idx]  # gets the layer information

            for ele in range(len(layer)):
                element = layer[ele]  # element data


                # remove numbers from symbol
                symbol = []
                for c in element[0]:
                    if not c.isdigit():
                        symbol.append(c)

                # recreates the chemical formula string
                symbol = ''.join(symbol)
                formula = formula + symbol

                if element[6] != '1':
                    formula = formula + element[6]

                thickness.append(float(element[1]))  # gets thickness data
                density.append(float(element[2]))  # gets density data
                roughness.append(float(element[3]))  # gets roughness data
                if element[4] != '' and element[4] != 'False' and element[4] != False:
                    linked_roughness.append(float(element[4]))
                else:
                    linked_roughness.append(False)

                if element[5][-1] != '[':
                    scat_fact.append(element[5])

                # scat_fact.append(element[5])
                # still need to take into account sf that are different than element

            sample.addlayer(idx, formula, thickness, density=density, roughness=roughness,
                            linked_roughness=linked_roughness, sf=scat_fact)

        # setting the energy shift value
        for key, value in self.eShift.items():
            if str(key).startswith('ff-'):
                sample.eShift[str(key)[3:]] = value
            elif str(key).startswith('ffm-'):
                sample.mag_eShift[str(key)[4:]] = value

        # setting the form factor scaling values
        for key, value in self.ffScale.items():
            if str(key).startswith('ff-'):
                sample.ff_scale[str(key)[3:]] = value
            elif str(key).startswith('ffm-'):
                sample.ffm_scale[str(key)[4:]] = value

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

                        if type(scattering_factor) is str:
                            sample.polymorphous(idx, ele_name, names, ratio, ast.literal_eval(scattering_factor))
                        else:
                            sample.polymorphous(idx, ele_name, names, ratio, scattering_factor)


        # determines the magnetization direction
        for idx in range(m):
            layer = self.structTableInfo[idx]  # gets the layer information
            for ele in range(len(layer)):
                ele_name = layer[ele][0]

                mag = self.magData[ele_name][idx]
                magDir = self.magDirection[idx]

                gamma = 0
                phi = 0
                if magDir == 'y':
                    gamma = 90
                    phi = 90
                elif magDir == 'x':
                    gamma = 0
                    phi = 90
                elif magDir == 'z':
                    gamma = 0
                    phi = 0

                names = mag[0]
                ratio = mag[1]
                scattering_factor = mag[2]

                isMagnetized = True
                nameBool = True
                notMagnetic = True
                Magnetic = True
                count = 0

                # loops through all the elements and determines if the layer is magnetized
                for i in range(len(scattering_factor)):
                    nameBool = True
                    notMagnetic = True
                    Magnetic = True
                    if names[i] == '':
                        nameBool = False

                    if not(ratio[i] == '' and scattering_factor[i] == ''):
                        notMagnetic = True
                    else:
                        count = count + 1

                    if not(ratio[i] != '' and scattering_factor[i] != ''):
                        Magnetic = False

                    if not(notMagnetic or Magnetic):
                        isMagnetized = False

                    if not(nameBool):
                        isMagnetized = False

                if count == len(scattering_factor):
                    isMagnetized = False


                #if ratio[0] != '' and names[0] != '' and scattering_factor[0] != '':
                if isMagnetized:

                    # handle case where not every polymorph is magnetic
                    my_idx = [j for j in range(len(ratio)) if ratio[j] != '' and (scattering_factor[j] != 0 and scattering_factor[j] != '')]

                    ratio = np.array(ratio)[my_idx]
                    names = np.array(names)[my_idx]
                    scattering_factor = np.array(scattering_factor)[my_idx]
                    ratio = [float(ratio[i]) for i in range(len(ratio))]

                    sample.magnetization(idx, names, ratio, scattering_factor, gamma=gamma,phi=phi)

                # setting phi and gamma for all elements in the layer
                # required for the construction of the density profile!

                sample.structure[idx][ele_name].gamma = gamma
                sample.structure[idx][ele_name].phi = phi

        # changing the form factors
        for idx in range(m):
            layer = self.structTableInfo[idx]  # gets the layer information
            for ele in range(len(layer)):
                ele_name = layer[ele][0]

                poly = self.varData[ele_name][idx]  # retrieves the element variation data for particular layer
                names = poly[0]

                if len(names) > 1:
                    if names[0] == '':
                        sample._set_form_factors(layer[ele][0], layer[ele][5])

        sample.energy_shift()  # this is done to make sure all the form factors are accounted for

        # retrieving the energy shift of the form factors
        for e in self.eShift.keys():
            if e.startswith('ff-'):
                key = e.strip('ff-')
                sample.eShift[key] = self.eShift[e]
                sample.ff_scale[key] = self.ffScale[e]
            elif e.startswith('ffm-'):
                key = e.strip('ffm-')
                sample.mag_eShift[key] = self.eShift[e]
                sample.ffm_scale[key] = self.ffScale[e]

        return sample

    def getData(self):
        """
        Purpose: retrieves the sample information for the slab class
        """
        # loops through each layer
        for j in range(len(self.sample.structure)):
            layer = self.sample.structure[j]  # retrieves the layer
            elekeys = list(layer.keys()) # retrieves the elements in the layer

            # loops through each element
            for ele in elekeys:

                if len(layer[ele].polymorph) != 0:  # checks if element is a polymorph
                    mag_density = ['' for i in range(len(layer[ele].polymorph))]
                    mag_sf = ['' for i in range(len(layer[ele].polymorph))]
                    self.varData[ele][j] = [layer[ele].polymorph, list(layer[ele].poly_ratio),
                                            layer[ele].scattering_factor]

                    if len(layer[ele].mag_density) != 0:
                        mag_density = list(layer[ele].mag_density)
                    if len(layer[ele].mag_scattering_factor) != 0:
                        mag_sf = layer[ele].mag_scattering_factor

                    # make sure that magData is in format that software can understand
                    for i in range(len(mag_sf)):
                        if mag_sf[i] == 0 or mag_sf[i] == '0':
                            mag_sf[i] = ''
                            mag_density[i] = ''
                    self.magData[ele][j] = [layer[ele].polymorph, mag_density, mag_sf]


                else:  # element is not a polymorph
                    mag_density = ['']
                    mag_sf = ['']
                    if len(layer[ele].mag_density) != 0:
                        mag_density = list(layer[ele].mag_density)
                    if len(layer[ele].mag_scattering_factor) != 0:
                        mag_sf = layer[ele].mag_scattering_factor

                        # implementation will change
                        for scat in mag_sf:
                            name = 'ffm-' + scat
                            if name not in list(self.eShift.keys()):
                                self.eShift[name] = 0
                                self.ffScale[name] = 1

                    self.magData[ele][j] = [[ele], mag_density, mag_sf]

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
    """
    Purpose: This widget handles the reflectivity workspace
    """
    def __init__(self, sWidget, data, data_dict, sim_dict):
        """
        :param sWidget: sampleWidget
        :param data: data information
        :param data_dict: data dictionary
        :param sim_dict: simulation dictionary
        """
        super().__init__()

        # --------------------------- Parameter Definition --------------------------------- #
        self.rom = [True, False, False]  # [reflectivity, optics magneto-optics]
        self.romSim = [True, False, False]  # same thing but for simulation
        self.isFit = True

        self.scanType = True  # if true do reflectivity scan

        self.bs = dict()  # background shift
        self.sf = dict()  # scaling factor

        self.sfBsFitParams = []  # variable used to keep track of the fitting parameters
        self.currentVal = []  # value of current fitting parameter

        self.sWidget = sWidget  # allows for reflectivityWidget to have access to sampleWidget information
        self.sample = sWidget.sample  # retrieves slab class form sampleWidget
        self.data = data  # information on the data
        self.data_dict = data_dict  # data info dictionary

        self.fit = []  # scans to be fit
        self.bounds = []  # bounds
        self.weights = []  # weights of the bounds

        self.isChangeTable = False
        self.axis_state = True
        self.scan_state = True
        self.previousIdx = 0

        # --------------------------- Layout Definition --------------------------------- #
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

        self.whichScan.setCurrentIndex(0)



        whichScanLayout = QHBoxLayout()

        whichScanLayout.addWidget(whichScanLabel)
        whichScanLayout.addWidget(self.whichScan)

        self.plot_scans()

        # button used to select which scans to fit
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
        self.reflectButton.setFixedWidth(150)
        self.energyButton = QPushButton('Energy Scan')
        self.energyButton.clicked.connect(self.changeToEnergyScan)
        self.energyButton.setStyleSheet('background: grey')
        self.energyButton.setFixedWidth(150)

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

        # layout and widgets for selected scans
        self.selectedScans = QComboBox()

        self.selectedScans.activated.connect(self.changeColorFit)
        self.whichScan.activated.connect(self.changeColorScan)
        #self.selectedScans.activated.connect(self.readTable)  # creating issues
        self.selectedScans.activated.connect(self.myPlotting)
        self.whichScan.activated.connect(self.myPlotting)
        self.selectedScans.activated.connect(self.setTable)
        #self.selectedScans.activated.connect(self.changeFitColor)

        # adding and removing boundaries
        self.addBoundaryButton = QPushButton('Add Boundary')
        self.addBoundaryButton.setFixedWidth(250)
        self.removeBoundaryButton = QPushButton('Remove Boundary')
        self.removeBoundaryButton.setMaximumWidth(250)
        self.removeScan = QPushButton('Remove Scan')
        self.removeScan.setMaximumWidth(250)
        self.removeScan.clicked.connect(self._removeScanSelection)
        self.removeScan.setStyleSheet("background-color : cyan")

        # radio button to use qz or angle
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

        # step size widget
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

        ERButtonLayout = QHBoxLayout()
        ERButtonLayout.addWidget(self.reflectButton)
        ERButtonLayout.addWidget(self.energyButton)

        # the line definitions are simply used to add line separators in the widget appearance
        line1 = QFrame()
        line1.setFrameShape(QFrame.HLine)
        line1.setLineWidth(3)

        line2 = QFrame()
        line2.setFrameShape(QFrame.HLine)
        line2.setLineWidth(3)

        line3 = QFrame()
        line3.setFrameShape(QFrame.HLine)
        line3.setLineWidth(3)

        line4 = QFrame()
        line4.setFrameShape(QFrame.HLine)
        line4.setLineWidth(3)

        sideline1 = QFrame()
        sideline1.setFrameShape(QFrame.VLine)
        sideline1.setLineWidth(3)

        sideline2 = QFrame()
        sideline2.setFrameShape(QFrame.VLine)
        sideline2.setLineWidth(3)

        scanSelectionLayout = QVBoxLayout()
        scanSelectionLayout.addWidget(line1)
        scanSelectionLayout.addSpacing(20)

        scanSelectionLayout.addLayout(whichScanLayout)
        scanSelectionLayout.addWidget(self.fitButton)
        scanSelectionLayout.addLayout(buttonLayout)

        scanSelectionLayout.addSpacing(20)
        scanSelectionLayout.addWidget(line2)
        scanSelectionLayout.addSpacing(20)

        scanSelectionLayout.addLayout(hbox)
        scanSelectionLayout.addLayout(stepLayout)

        scanSelectionLayout.addSpacing(20)
        scanSelectionLayout.addWidget(line3)
        scanSelectionLayout.addSpacing(20)

        scanSelectionLayout.addWidget(simulationLabel)
        scanSelectionLayout.addLayout(simButtonLayout)
        scanSelectionLayout.addSpacing(20)

        scanSelectionLayout.addLayout(ERButtonLayout)
        scanSelectionLayout.addLayout(simAngleLayout)
        scanSelectionLayout.addLayout(simEnergyLayout)
        scanSelectionLayout.addLayout(energyRangeLayout)
        scanSelectionLayout.addWidget(self.polarBox)
        scanSelectionLayout.addSpacing(20)
        scanSelectionLayout.addWidget(line4)

        scanLayout = QHBoxLayout()
        scanLayout.addWidget(sideline1)
        scanLayout.addLayout(scanSelectionLayout)
        scanLayout.addWidget(sideline2)


        toplayout = QHBoxLayout()
        toplayout.addWidget(self.spectrumWidget)
        toplayout.addLayout(scanLayout)

        # setting up the scan boundary workspace -------------------------------------------------------------
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

        # setting up scaling factor and background shift layout and widgets
        sfbsLayout = QVBoxLayout()

        bsLayout = QHBoxLayout()
        bsLabel = QLabel('Background Shift:')
        bsLabel.setFixedWidth(90)
        self.backgroundShift = QLineEdit()
        self.scalingFactor = QLineEdit()
        self.backgroundShift.editingFinished.connect(self.changeSFandBS)
        self.backgroundShift.setFixedWidth(100)
        self.backgroundShift.setFixedHeight(25)
        bsLayout.addWidget(bsLabel)
        bsLayout.addWidget(self.backgroundShift)

        sfLayout = QHBoxLayout()
        sfLabel = QLabel('Scaling Factor:')
        sfLabel.setFixedWidth(90)
        # self.scalingFactor.textChanged.connect(self.changeSFandBS)
        self.scalingFactor.editingFinished.connect(self.changeSFandBS)
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
        pagelayout.addSpacing(10)

        pagelayout.addLayout(bottomlayout)

        self.backgroundShift.editingFinished.connect(self.bsChange)
        self.scalingFactor.editingFinished.connect(self.sfChange)
        self.setLayout(pagelayout)

    def changeToReflScan(self):
        """
        Purpose: Plot reflectivity scan simulation
        """
        self.energyButton.setStyleSheet('background: grey')
        self.reflectButton.setStyleSheet('background: cyan')
        self.scanType = True
        self.mySimPlotting()

    def changeToEnergyScan(self):
        """
        Purpose: Plot energy scan simulation
        """
        self.energyButton.setStyleSheet('background: cyan')
        self.reflectButton.setStyleSheet('background: grey')
        self.scanType = False
        self.mySimPlotting()

    def rPlotSim(self):
        """
        Purpose: Plot the simulation scan
        """
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
            Theta = np.linspace(0.1, 89.1, num=n)
            qz = np.sin(Theta * np.pi / 180) * (energy * 0.001013546143)
            qz, R = self.sample.reflectivity(energy, qz, s_min=step_size)
            R = R[polName]
            if self.axis_state:  # momentum transfer
                self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 1), width=2), name='Simulation')
                self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
            else:  # angle
                self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0, 1), width=2), name='Simulation')
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
            E = np.linspace(lw, up, num=n)
            E, R = self.sample.energy_scan(angle, E, s_min=step_size)
            R = R[polName]

            self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 1), width=2), name='Simulation')
            self.spectrumWidget.setLabel('bottom', "Energy (eV)")
            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            self.spectrumWidget.setLogMode(False, False)

    def opPlotSim(self):
        """
        Purpose: plot the optical profile
        """
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
            scale = float(self.sWidget.ffScale[name])
            sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E + dE, False)*scale
        # Magnetic Scattering Factor
        for em in self.sample.find_sf[1].keys():
            name = 'ffm-' + self.sample.find_sf[1][em]
            dE = float(self.sWidget.eShift[name])
            scale = float(self.sWidget.ffScale[name])
            sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E + dE, True)*scale

        delta, beta = ms.index_of_refraction(density, sf, E)  # calculates dielectric constant for structural component

        self.spectrumWidget.plot(thickness, delta, pen=pg.mkPen((0, 2), width=2), name='delta')
        self.spectrumWidget.plot(thickness, beta, pen=pg.mkPen((1, 2), width=2), name='beta')
        self.spectrumWidget.setLabel('left', "Reflectivity, R")
        self.spectrumWidget.setLabel('bottom', "Thickness, Å")
        self.spectrumWidget.setLogMode(False, False)

    def opmPlotSim(self):
        """
        Purpose: plot the magnetic optical profile
        """
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
            scale = float(self.sWidget.ffScale[name])
            sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E + dE, False)*scale
        # Magnetic Scattering Factor
        for em in self.sample.find_sf[1].keys():
            name = 'ffm-' + self.sample.find_sf[1][em]
            dE = float(self.sWidget.eShift[name])
            scale = float(self.sWidget.ffScale[name])
            sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E + dE, True)*scale

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
        """
        Purpose: Makes sure angle is within the appropriate range
        """
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
        """
        Purpose: Makes sure the energy is set within the appropriate range
        """
        energy = self.simEnergyEdit.text()

        if energy != '':
            energy = float(energy)
            lower = str(energy - 50)
            upper = str(energy + 50)
            self.simLowEnergy.setText(lower)
            self.simUpEnergy.setText(upper)
            self.mySimPlotting()

    def setLowerEnergy(self):
        """
        Purpose: set the lower energy for the simulation plot
        """
        energy = float(self.simEnergyEdit.text())
        lower = float(self.simLowEnergy.text())

        if energy < lower:
            lower = energy - 50

        self.simLowEnergy.setText(str(lower))
        self.mySimPlotting()

    def setUpperEnergy(self):
        """
        Purpose: set the upper energy for the simulation plot
        """
        energy = float(self.simEnergyEdit.text())
        upper = float(self.simUpEnergy.text())

        if energy > upper:
            upper = energy + 50

        self.simUpEnergy.setText(str(upper))
        self.mySimPlotting()

    def value_changed(self):
        """
        Purpose: Changes scan boundaries
        """

        idx = self.selectedScans.currentIndex()
        row = self.boundWeightTable.currentRow()
        column = self.boundWeightTable.currentColumn()
        item = self.boundWeightTable.currentItem().text()
        if row == 0 or row == 1:  # upper or lower bounds
            self.bounds[idx][column][row] = item
        elif row == 2:  # weights
            self.weights[idx][column] = item

    def bsChange(self):
        """
        Purpose: change the background shift when signaled
        """
        name = self.selectedScans.currentText()  # scan name
        idx = self.selectedScans.currentIndex()  # scan index
        var = self.backgroundShift.text()  # value
        copy_of_fit = copy.deepcopy(self.sfBsFitParams)
        old_var = ''
        if name != '':
            # first change the value
            if self.allScan.checkState() == 0:  # not for all scans
                old_var = self.bs[name]
                self.bs[name] = var
                self.data_dict[name]['Background Shift'] = float(var)
            else:  # for all scans
                old_var = self.bs[name]
                for key in list(self.bs.keys()):
                    self.bs[key] = var
                    self.data_dict[name]['Background Shift'] = float(var)

            # update fitting parameters
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
        """
        Purpose: Change the scaling factor when signaled
        """
        name = self.selectedScans.currentText()  # scan name
        idx = self.selectedScans.currentIndex()  # scan index
        var = self.scalingFactor.text()  # value
        copy_of_fit = copy.deepcopy(self.sfBsFitParams)

        if name != '':
            # first change the value
            if self.allScan.checkState() == 0:  # not all scans
                self.sf[name] = var
                self.data_dict[name]['Scaling Factor'] = float(var)
            else:  # all scans
                for key in list(self.sf.keys()):
                    self.sf[key] = var
                    self.data_dict[name]['Scaling Factor'] = float(var)

            # update fitting parameters
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
        """
        Purpose: change the color of working space if parameter is to be fit
        """

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
        """
        Purpose: erase the fitting parameters previously input by the user
        """
        self.sfBsFitParams = []
        self.changeFitColor()

    def fitBackgroundShift(self):
        """
        Purpose: Add background shift to fitting parameters
        """
        state = self.allScan.checkState()  # fit all scans or current scan
        name = self.selectedScans.currentText()  # name of current scan
        value = self.backgroundShift.text()  # current value of background shift
        if name != '':
            if state == 2:  # all scans
                fit = ['BACKGROUND SHIFT', 'ALL SCANS']
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 5e-8
                    upper = float(value) + 5e-8
                    self.currentVal.append([float(value), [lower, upper]])
            else:  # current scan
                fit = ['BACKGROUND SHIFT', name]
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 5e-8
                    upper = float(value) + 5e-8
                    self.currentVal.append([float(value), [lower, upper]])

        self.changeFitColor()

    def fitScalingFactor(self):
        """
        Purpose: Add scaling factor to fitting parameters
        """
        state = self.allScan.checkState()  # all scans (2) or current scans (0)
        name = self.selectedScans.currentText()  # name of current scan
        value = self.scalingFactor.text()  # name of current value
        if name != '':
            if state == 2:  # all scans
                fit = ['SCALING FACTOR', 'ALL SCANS']
                if fit not in self.sfBsFitParams:
                    self.sfBsFitParams.append(fit)
                    lower = float(value) - 0.2
                    if lower < 0:
                        lower = 0
                    upper = float(value) + 0.2
                    self.currentVal.append([float(value), [lower, upper]])
            else:  # current scan
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
        """
        Purpose: Remove background shift from fitting parameters
        """
        state = self.allScan.checkState()  # all scans or current scan
        name = self.selectedScans.currentText()  # current scan name
        if name != '':
            if state == 2:  # remove fit for all scans
                if ['BACKGROUND SHIFT', 'ALL SCANS'] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['BACKGROUND SHIFT', 'ALL SCANS'])
                    self.sfBsFitParams.remove(['BACKGROUND SHIFT', 'ALL SCANS'])
                    self.currentVal.pop(idx)
            else:  # remove fit for current scan
                if ['BACKGROUND SHIFT', name] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['BACKGROUND SHIFT', name])
                    self.sfBsFitParams.remove(['BACKGROUND SHIFT', name])
                    self.currentVal.pop(idx)
        self.changeFitColor()

    def unfitScalingFactor(self):
        """
        Purpose: Remove scaling factor from fitting parameters
        """
        state = self.allScan.checkState()  # all scans or current scan
        name = self.selectedScans.currentText()  # name of scan
        if name != '':
            if state == 2:  # remove all scans from fit
                if ['SCALING FACTOR', 'ALL SCANS'] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['SCALING FACTOR', 'ALL SCANS'])
                    self.sfBsFitParams.remove(['SCALING FACTOR', 'ALL SCANS'])
                    self.currentVal.pop(idx)
            else:  # remove current scan from fit
                if ['SCALING FACTOR', name] in self.sfBsFitParams:
                    idx = self.sfBsFitParams.index(['SCALING FACTOR', name])
                    self.sfBsFitParams.remove(['SCALING FACTOR', name])
                    self.currentVal.pop(idx)

        self.changeFitColor()

    def changeSFandBS(self):
        """
        Purpose: change scaling factor and background shift when signaled
        """
        idx = self.selectedScans.currentIndex()  # current index
        name = self.selectedScans.currentText()  # current scan name
        bs = self.backgroundShift.text()  # current background shift
        sf = self.scalingFactor.text()  # current scaling factor
        if bs != '' and sf != '':  # checks to make sure values are not empty
            if self.allScan.checkState() == 0:  # case where all scans have different bs and sf
                self.bs[name] = bs  # set background shift
                self.sf[name] = sf  # set scaling factor
                self.data_dict[name]['Background Shift'] = float(bs)  # update data_dict
                self.data_dict[name]['Scaling Factor'] = float(sf)
            else:  # case where all scans have same bs and sf
                for key in list(self.bs.keys()):
                    self.data_dict[key]['Background Shift'] = float(bs)
                    self.data_dict[key]['Scaling Factor'] = float(sf)

    def changeStepSize(self):
        """
        Purpose: Change the step size
        """
        self.sWidget._step_size = self.stepWidget.text()

    def mySimPlotting(self):
        """
        Purpose: determine to plot reflectometry, optical profile, or magnetic optical profile for simulation
        :return:
        """
        self.isFit = False
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        # self.sample = self.sWidget.sample
        idx = self.romSim.index(True)
        if idx == 0:  # reflectometry
            self.rPlotSim()
        elif idx == 1: # optical profile
            self.opPlotSim()
        elif idx == 2:  # magneto-optical profile
            self.opmPlotSim()

    def myPlotting(self):
        """
        Purpose: determine to plot relfectometry, optical profile, magneto-optical profile
        :return:
        """

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')
        self.isFit = True

        # self.sample = self.sWidget.sample
        idx = self.rom.index(True)
        if idx == 0:  # relfectometry
            self.rPlot()
        elif idx == 1:  # optical profile
            self.opPlot()
        elif idx == 2:  # magneto-optical profile
            self.opmPlot()

    def rPlot(self):
        """
        Purpose: plot reflectometry
        """
        self.rom = [True, False, False]

        self.rButton.setStyleSheet('background: cyan')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        if self.scan_state:  # plot from all scans
            self.plot_scans()
        else:  # plot from selected scans
            self.plot_selected_scans()
            self.setTable()

    def opPlot(self):
        """
        Purpose: plot optical profile
        """
        self.rom = [False, True, False]
        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: cyan')
        self.opmButton.setStyleSheet('background: grey')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        name = ''
        self.spectrumWidget.clear()

        if self.scan_state:  # from all scans
            name = self.whichScan.currentText()
        else:  # from selected scans
            name = self.selectedScans.currentText()
        if name != '':
            E = self.data_dict[name]['Energy']  # scan energy

            step_size = float(self.sWidget._step_size)  # get step size
            thickness, density, density_magnetic = self.sample.density_profile(
                step=step_size)  # Computes the density profile

            sf = dict()  # form factors of non-magnetic components
            sfm = dict()  # form factors of magnetic components

            # Non-Magnetic Scattering Factor
            for e in self.sample.find_sf[0].keys():
                name = 'ff-' + self.sample.find_sf[0][e]
                dE = float(self.sWidget.eShift[name])
                scale = float(self.sWidget.ffScale[name])
                sf[e] = ms.find_form_factor(self.sample.find_sf[0][e], E + dE, False)*scale
            # Magnetic Scattering Factor
            for em in self.sample.find_sf[1].keys():
                name = 'ffm-' + self.sample.find_sf[1][em]
                dE = float(self.sWidget.eShift[name])
                scale = float(self.sWidget.ffScale[name])
                sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E + dE, True)*scale

            delta, beta = ms.index_of_refraction(density, sf,
                                                 E)  # calculates dielectric constant for structural component

            self.spectrumWidget.plot(thickness, delta, pen=pg.mkPen((0, 2), width=2), name='delta')
            self.spectrumWidget.plot(thickness, beta, pen=pg.mkPen((1, 2), width=2), name='beta')
            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            self.spectrumWidget.setLabel('bottom', "Thickness, Å")
            self.spectrumWidget.setLogMode(False, False)
            # delta_m, beta_m = ms.magnetic_optical_constant(density_magnetic, sfm, E)  # calculates dielectric constant for magnetic component

    def opmPlot(self):
        """
        Purpose: plot magneto-optical profile
        """
        self.rom = [False, False, True]

        self.rButton.setStyleSheet('background: grey')
        self.opButton.setStyleSheet('background: grey')
        self.opmButton.setStyleSheet('background: cyan')

        self.rButtonSim.setStyleSheet('background: grey')
        self.opButtonSim.setStyleSheet('background: grey')
        self.opmButtonSim.setStyleSheet('background: grey')

        name = ''
        self.spectrumWidget.clear()

        if self.scan_state:  # from all scans
            name = self.whichScan.currentText()
        else:  # from selected scans
            name = self.selectedScans.currentText()

        if name != '':
            E = self.data_dict[name]['Energy']  # energy of scan

            step_size = float(self.sWidget._step_size)
            thickness, density, density_magnetic = self.sample.density_profile(
                step=step_size)  # Computes the density profile

            sf = dict()  # form factors of non-magnetic components
            sfm = dict()  # form factors of magnetic components

            # Non-Magnetic Scattering Factor

            for em in self.sample.find_sf[1].keys():
                sfm[em] = ms.find_form_factor(self.sample.find_sf[1][em], E, True)

            delta_m, beta_m = ms.magnetic_optical_constant(density_magnetic, sfm,
                                                           E)  # calculates dielectric constant for magnetic component
            self.spectrumWidget.plot(thickness, delta_m, pen=pg.mkPen((0, 2), width=2), name='delta_m')
            self.spectrumWidget.plot(thickness, beta_m, pen=pg.mkPen((1, 2), width=2), name='beta_m')
            self.spectrumWidget.setLabel('left', "Reflectivity, R")
            self.spectrumWidget.setLabel('bottom', "Thickness, Å")
            self.spectrumWidget.setLogMode(False, False)

    def updateAxis(self):
        """
        Purpose: Update x-axis for qz and angle
        """
        self.readTable()
        rbtn = self.sender()

        if rbtn.isChecked() == True:
            if rbtn.text() == 'qz (A)':  # axis qz
                self.axis_state = True
            else:  # axis theta
                self.axis_state = False

        # plot reflectometry
        if self.isFit:
            self.myPlotting()
        else:
            self.mySimPlotting()

    def changeColorScan(self):
        """
        Purpose: Changes all scan combobox to red demonstrating that the current scan being shown is from that selection
        """
        self.selectedScans.setStyleSheet('background: white; selection-background-color: grey')
        self.whichScan.setStyleSheet('background: red; selection-background-color: red')
        self.scan_state = True

    def changeColorFit(self):
        """
        Purpose: Change color of comboBox showing the fit combobox scan is currently being shown on the plot
        """
        self.selectedScans.setStyleSheet('background: red; selection-background-color: red')
        self.whichScan.setStyleSheet('background: white; selection-background-color: grey')
        self.scan_state = False

    def setTable(self):
        """
        Purpose: set background shift and scaling factor table
        """

        self.isChangeTable = True
        self.boundWeightTable.blockSignals(True)
        idx = self.selectedScans.currentIndex()  # index of scan

        name = self.selectedScans.currentText()  # name of scan

        if name != '':
            self.scalingFactor.setText(self.sf[name])  # setting the appropriate scaling factor
            self.backgroundShift.setText(self.bs[name])  # setting the appropriate background shift

            E = self.data_dict[name]['Energy']  # energy of scan
            mykeys = list(self.data_dict[name].keys())

            bound = self.bounds[idx]
            weight = self.weights[idx]
            col = len(bound)

            row = 3

            self.boundWeightTable.setRowCount(row)
            self.boundWeightTable.setColumnCount(col)

            self.boundWeightTable.setVerticalHeaderLabels(['Lower Bound', 'Upper Bound', 'Weight'])

            # loop through all columns and rows
            for i in range(row):
                for j in range(col):
                    if i == 0:
                        myitem = ''
                        if 'Angle' not in mykeys:
                            if not self.axis_state:
                                if len(bound[j][0]) != 0:  # show boundary as an angle
                                    myitem = str(np.arcsin(float(bound[j][0]) / (E * 0.001013546143)) * 180 / np.pi)[:7]
                            else:
                                myitem = copy.copy(bound[j][0])
                        else:
                            myitem = copy.copy(bound[j][0])

                        item = QTableWidgetItem(str(myitem))
                        self.boundWeightTable.setItem(i, j, item)

                        # self.boundWeightTable.setItem(i, j, item)
                    elif i == 1:
                        myitem = ''
                        if 'Angle' not in mykeys:
                            if not self.axis_state:
                                if len(bound[j][1]) != 0:  # show boundary as an angle
                                    myitem = str(np.arcsin(float(bound[j][1]) / (E * 0.001013546143)) * 180 / np.pi)[:7]
                            else:
                                myitem = copy.copy(bound[j][1])
                        else:
                            myitem = copy.copy(bound[j][1])

                        item = QTableWidgetItem(str(myitem))
                        self.boundWeightTable.setItem(i, j, item)
                    elif i == 2:

                        item = QTableWidgetItem(str(weight[j]))
                        self.boundWeightTable.setItem(i, j, item)
        self.isChangeTable = False
        self.boundWeightTable.blockSignals(False)

    def _scanSelection(self):
        """
        Purpose: add scan to scan selection list when signaled
        """
        idx = self.whichScan.currentIndex()  # index of current scan
        name = self.whichScan.currentText()  # name of current scan

        if name not in self.fit:  # checks if scan already selected
            # pre-initialize boundary
            if 'Angle' in list(self.data_dict[name].keys()):
                lower = str(self.data_dict[name]['Data'][3][0])[0:7]
                upper = str(self.data_dict[name]['Data'][3][-1])[0:7]
            else:
                lower = str(self.data_dict[name]['Data'][0][0])[0:7]
                upper = str(self.data_dict[name]['Data'][0][-1])[0:7]

            self.fit.append(name)  # add to fitting list
            self.bounds.append([[lower, upper]])  # initialize boundary
            self.weights.append(['1'])  # add pre-set weight
            self.selectedScans.addItem(name)

            # sets the bs and sf values  (background shift and scaling factor)
            if self.allScan.checkState() == 2:  # all scans state checked
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

            m = len(self.fit)

            if m != 0:
                my_idx = m-1
                self.selectedScans.setCurrentIndex(my_idx)
                self.setTable()

    def _removeScanSelection(self):
        """
        Purpose: Remove scan from selection
        """
        idx = self.selectedScans.currentIndex()  # current scan index
        name = self.selectedScans.currentText()  # scan name

        if name != '':
            # takes care of proper indexing
            if idx == self.previousIdx and self.previousIdx != 0:
                self.previousIdx = self.previousIdx - 1

            self.selectedScans.removeItem(idx)  # selected scans case where all scans have same bs and sf
            self.fit.pop(idx)
            self.bounds.pop(idx)
            self.weights.pop(idx)

            # remove background shift and scaling factor
            del self.bs[name]
            del self.sf[name]

            if len(self.fit) != 0:
                self.setTable()  # makes sure that the table is switched
                self.myPlotting()
            else:
                self.spectrumWidget.clear()


    def addBoundWeight(self):
        """
        Purpose: Add boundary weight
        """
        col = self.boundWeightTable.columnCount()
        idx = self.selectedScans.currentIndex()  # gets scan index
        n = len(self.bounds[idx])
        upper = self.bounds[idx][n - 1][1]  # gets the last boundary
        self.bounds[idx][n - 1][1] = ''
        self.bounds[idx].append(['', upper])
        self.weights[idx].append('1')
        self.boundWeightTable.setColumnCount(col + 1)

        self.setTable()

    def removeBoundWeight(self):
        """
        Purpose: Remove boundary weight
        """
        col = self.boundWeightTable.columnCount()  # current column
        idx = self.selectedScans.currentIndex()  # gets the selected scan

        if col != 1:
            n = len(self.bounds[idx])  # get the number of boundaries
            upper = self.bounds[idx][n - 1][1]  # gets the proper upper boundary
            self.bounds[idx][n - 2][1] = upper
            self.bounds[idx].pop()
            self.weights[idx].pop()
            self.boundWeightTable.setColumnCount(col - 1)

        self.setTable()

    def readTable(self):
        """
        Purpose: Read the table for scan boundaries
        """
        idx = self.selectedScans.currentIndex()  # current scan index
        name = self.selectedScans.currentText()  # current scan name

        row = self.boundWeightTable.rowCount()  # current row
        column = self.boundWeightTable.columnCount()  # current column
        if name != '':
            E = self.data_dict[name]['Energy']
            for i in range(row):
                for j in range(column):

                    item = self.boundWeightTable.item(i, j).text()
                    if i == 0:
                        if 'Angle' not in list(self.data_dict[name].keys()):
                            if not (self.axis_state):  # we are in angle state
                                if len(item) != 0:
                                    item = str(np.sin(float(item) * np.pi / 180) * (E * 0.001013546143))
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
                                    item = str(np.sin(float(item) * np.pi / 180) * (E * 0.001013546143))

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
        """
        Purpose: plot scans from whichScan
        """
        self.spectrumWidget.clear()
        if len(self.data) != 0:
            # self.sample = self.sWidget.sample

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

                    qz, Rsim = self.sample.reflectivity(E, qz, s_min=step_size, bShift=background_shift,
                                                        sFactor=scaling_factor)
                    Theta = np.arcsin(qz / (E * 0.001013546143)) * 180 / np.pi
                    Rsim = Rsim[pol]
                    n = len(qz)
                    if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                        if self.axis_state:
                            self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                            self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((2, 3), width=2), name='Simulation')

                        else:
                            self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                            self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((2, 3), width=2), name='Simulation')

                        self.spectrumWidget.setLabel('left', "Reflectivity, R")
                        self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                        self.spectrumWidget.setLogMode(False, True)

                    elif pol == 'AL' or pol == 'AC':
                        rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                        if self.axis_state:

                            self.spectrumWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2), name='Data')
                            self.spectrumWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((2, 3), width=2),
                                                     name='Simulation')
                        else:
                            self.spectrumWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2),
                                                     name='Data')
                            self.spectrumWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((2, 3), width=2),
                                                     name='Simulation')

                        self.spectrumWidget.setLogMode(False, False)
                        self.spectrumWidget.setLabel('left', "Reflectivity, R")
                        self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")

                elif scan_type == 'Energy':
                    E = dat[3]
                    R = dat[2]
                    Theta = self.data_dict[name]['Angle']
                    E, Rsim = self.sample.energy_scan(Theta, E, s_min=step_size, bShift=background_shift,
                                                      sFactor=scaling_factor)
                    Rsim = Rsim[pol]
                    self.spectrumWidget.setLogMode(False, False)
                    self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                    self.spectrumWidget.plot(E, Rsim, pen=pg.mkPen((2, 3), width=2), name='Simulation')
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Energy, E (eV)")
        self.spectrumWidget.enableAutoRange()  # resets the range such that we view everything

    def plot_selected_scans(self):
        """
        Purpose: Plot scan from selectedScans
        """
        step_size = float(self.sWidget._step_size)
        self.sample = self.sWidget._createSample()
        self.spectrumWidget.clear()
        name = self.selectedScans.currentText()
        b_idx = self.selectedScans.currentIndex()

        if name != '':
            bound = self.bounds[b_idx]
            lower = float(bound[0][0])
            upper = float(bound[-1][-1])
            background_shift = self.data_dict[name]['Background Shift']
            scaling_factor = self.data_dict[name]['Scaling Factor']

            idx = 0
            notDone = True
            while notDone and idx == len(self.data) - 1:
                temp_name = self.data[idx][2]
                if temp_name == name:
                    notDone = False
                else:
                    idx = idx + 1

            dat = self.data_dict[name]['Data']
            pol = self.data_dict[name]['Polarization']

            if 'Angle' not in list(self.data_dict[name].keys()):
                qz = dat[0]
                R = dat[2]
                E = self.data_dict[name]['Energy']
                qz, Rsim = self.sample.reflectivity(E, qz, s_min=step_size, bShift=background_shift,
                                                    sFactor=scaling_factor)
                Theta = np.arcsin(qz / (E * 0.001013546143)) * 180 / np.pi

                Rsim = Rsim[pol]
                if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                    if self.axis_state:
                        self.spectrumWidget.plot(qz, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.spectrumWidget.plot(qz, Rsim, pen=pg.mkPen((2, 3), width=2), name='Simulation')
                    else:
                        self.spectrumWidget.plot(Theta, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.spectrumWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                        lower = np.arcsin(lower / (E * 0.001013546143)) * 180 / np.pi
                        upper = np.arcsin(upper / (E * 0.001013546143)) * 180 / np.pi
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.spectrumWidget.setLogMode(False, True)

                    self.spectrumWidget.setXRange(lower, upper)
                elif pol == 'AL' or pol == 'AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    if self.axis_state:
                        self.spectrumWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.spectrumWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((2, 3), width=2),
                                                 name='Simulation')
                    else:
                        self.spectrumWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.spectrumWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((2, 3), width=2),
                                                 name='Simulation')

                    self.spectrumWidget.setLogMode(False, False)
                    self.spectrumWidget.setLabel('left', "Reflectivity, R")
                    self.spectrumWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.spectrumWidget.setXRange(lower, upper)
            else:
                E = dat[3]
                R = dat[2]
                Theta = self.data_dict[name]['Angle']
                E, Rsim = self.sample.energy_scan(Theta, E, s_min=step_size, bShift=background_shift,
                                                  sFactor=scaling_factor)
                Rsim = Rsim[pol]
                self.spectrumWidget.setLogMode(False, False)
                self.spectrumWidget.plot(E, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                self.spectrumWidget.plot(E, Rsim, pen=pg.mkPen((2, 3), width=2), name='Simulation')
                self.spectrumWidget.setLabel('left', "Reflectivity, R")
                self.spectrumWidget.setLabel('bottom', "Energy, E (eV)")


class Worker(QObject):
    """
    Purpose: Worker used to allow for GUI not to freeze while data fitting running
    """
    finished = pyqtSignal()  # signal that process finished
    progress = pyqtSignal(int)  # signal used for progress

    def __init__(self, function):
        super().__init__()
        self.function = function  # globalOptimizationWidget

    def run(self):
        x, fun = self.function._optimizer()  # run the data fitting function
        self.function.x = x  # parameters
        self.function.fun = fun  # cost function
        self.finished.emit()  # let process know that data fitting has terminated


class UpdateWorker(QObject):
    """
    Purpose: Worker used to update cost function after each data fitting callback
    """
    finished = pyqtSignal()  # program finished
    progress = pyqtSignal(int)  # update

    def __init__(self, function):
        super().__init__()
        self.function = function  # globalOptimizationWidget

    def run(self):
        # first check script
        # if script good continue, otherwise apport
        self.function.update_optimization()  # run update cost function
        self.finished.emit()  # process finished

    def stop(self):
        self.function.stop()  # stop process when signaled


class callback():
    """
    Purpose: terminate data fitting when user chooses
    """
    def __init__(self):
        self.Finish = False

    def stop_evolution(self, x, convergence):
        # end differential evolution properly
        x_vars.append(x)
        if stop:
            return True
        else:
            return False

    def stop_simplicial(self, x):
        # end simplicial homology properly
        x_vars.append(x)
        if stop:
            return True
        else:
            return False

    def stop_annealing(self, x, f, connect):
        # end simulated annealing properly
        x_vars.append(x)
        if stop:
            return True
        else:
            return False


class GlobalOptimizationWidget(QWidget):
    """
    Purpose: Widget used to setup a data fitting
    """
    def __init__(self, sWidget, rWidget, nWidget, pWidget, rApp):
        super().__init__()


        # ------------------------- Parameter Definitions ----------------------------------#

        self.sWidget = sWidget  # sampleWidget
        self.rWidget = rWidget  # reflectivityWidget
        self.nWidget = nWidget  # smoothingWidget
        self.pWidget = pWidget  # progressWidget
        self.rApp = rApp  # Application widget

        self.sample = copy.deepcopy(self.sWidget.sample)  # update slab class
        self.temp_sample = copy.deepcopy(self.sample)  # temporary slab class
        self.sampleBounds = []  # sample boundaries
        self.sfBounds = []  # scattering factor boundaries
        self.otherBounds = []  # other boundaries

        self.x = []  # parameters fitting values
        self.fun = 0  # cost function
        self.callback = callback()  # callback function
        self.progressFinished = True  # has progress finished
        self.objective = 'Chi-Square'  # initialized objective function
        self.shape_weight = 0  # initialized total variation weight


        plotLayout = QHBoxLayout()  # plotting layout

        self.goParameters = {
            'differential evolution': ['currenttobest1bin', 2, 15, 1e-6, 0, 0.5, 1, 0.7, True, 'latinhypercube',
                                       'immediate'],
            'simplicial homology': ['None', 1, 'simplicial'],
            'dual annealing': [150, 5230.0, 2e-5, 2.62, 5.0, 10000000.0, True],
            'least squares': ['2-point', 'trf', 1e-8, 1e-8, 1e-8, 1.0, 'linear', 1.0, 'None', 'None']}

        # Uncomment this is the direct algorithm is used
        """
        self.goParameters = {
            'differential evolution': ['currenttobest1bin', 2, 15, 1e-6, 0, 0.5, 1, 0.7, True, 'latinhypercube',
                                       'immediate'],
            'simplicial homology': ['None', 1, 'simplicial'],
            'dual annealing': [150, 5230.0, 2e-5, 2.62, 5.0, 10000000.0, True],
            'least squares': ['2-point', 'trf', 1e-8, 1e-8, 1e-8, 1.0, 'linear', 1.0, 'None', 'None'],
            'direct': [0.0001, 'None', 1000, False, 0.0001,1e-16,1e-6]}
        """

        # ------------------------------- Layout Definition -------------------------------- #

        # determine R transformation
        isLogLayout = QVBoxLayout()
        isLogLabel = QLabel('Optimization Scale:')
        isLogLabel.setFixedWidth(200)
        self.isLogWidget = QComboBox()
        self.isLogWidget.addItems(['log(x)', 'ln(x)', 'x', 'qz^4'])
        self.isLogWidget.currentIndexChanged.connect(self._set_y_scale)
        isLogLayout.addWidget(isLogLabel)
        isLogLayout.addWidget(self.isLogWidget)

        self.paramChange = True
        self.parameters = []
        selectedScansLayout = QVBoxLayout()
        self.selectedScans = QComboBox()  # shows the scans selected for the data fitting process
        self.selectedScans.activated.connect(self.plot_scan)
        selectedScansLabel = QLabel('Fit Scans:')
        selectedScansLabel.setFixedWidth(200)

        selectedScansLayout.addWidget(selectedScansLabel)
        selectedScansLayout.addWidget(self.selectedScans)

        # Adding the plotting Widget
        self.plotWidget = pg.PlotWidget()
        self.plotWidget.setBackground('w')

        self.plotWidget.addLegend()


        # Global optimization parameters and fitting
        buttonLayout = QVBoxLayout()

        mylayout = QHBoxLayout()
        self.checkBoxLabel = QLabel('Run Script: ')
        self.checkBoxLabel.setFixedWidth(60)
        self.checkBox = QCheckBox()
        mylayout.addWidget(self.checkBoxLabel)
        mylayout.addSpacing(5)
        mylayout.addWidget(self.checkBox)


        self.runButton = QPushButton('Run Optimization')
        self.runButton.pressed.connect(self._run_global_optimization)
        self.runButton.setStyleSheet('background: green')

        self.stopButton = QPushButton('Stop Optimization')
        self.stopButton.clicked.connect(self._stop_optimization)

        self.optButton = QPushButton('Update Sample')
        self.optButton.clicked.connect(self._save_optimization)
        self.optButton.setStyleSheet('background: cyan')

        self.clearFitButton = QPushButton('Clear Fit')
        self.clearFitButton.clicked.connect(self._clear_fit)

        buttonLayout.addLayout(selectedScansLayout)
        buttonLayout.addLayout(isLogLayout)
        buttonLayout.addStretch(1)
        buttonLayout.addLayout(mylayout)
        buttonLayout.addWidget(self.runButton)
        buttonLayout.addWidget(self.stopButton)
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(self.optButton)
        buttonLayout.addWidget(self.clearFitButton)

        # Adding Widgets to plotting layout
        plotLayout.addWidget(self.plotWidget, 5)

        plotLayout.addLayout(buttonLayout)

        self.fittingParamTable = QTableWidget()
        self.fittingParamTable.setColumnCount(5)
        self.fittingParamTable.setHorizontalHeaderLabels(
            ['Name', 'Current Value', 'Lower Boundary', 'Upper Boundary', 'New'])
        delegate = ReadOnlyDelegate()
        self.fittingParamTable.setItemDelegateForColumn(0, delegate)
        self.fittingParamTable.setItemDelegateForColumn(1, delegate)
        self.fittingParamTable.setItemDelegateForColumn(4, delegate)
        tableLayout = QVBoxLayout()
        tableLayout.addWidget(self.fittingParamTable)
        self.fittingParamTable.itemChanged.connect(self._changeFitVar)
        self.fittingParamTable.viewport().installEventFilter(self)

        # Include the different input parameters for the global optimization and their algorithms
        self.goParamWidget = QWidget()

        # Adding objective function parameters

        # start with the layout of the "main" window
        algorithmLayout = QVBoxLayout()
        algorithmLabel = QLabel('Algorithm Selection')
        self.algorithmSelect = QComboBox()
        #self.algorithmSelect.addItems(
        #    ['differential evolution', 'simplicial homology', 'dual annealing', 'least squares', 'direct'])
        self.algorithmSelect.addItems(
            ['differential evolution', 'simplicial homology', 'dual annealing', 'least squares'])
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
        costButton = QPushButton('Cost Value')
        costButton.setFixedWidth(80)
        costButton.clicked.connect(self.calculateCost)
        self.costValue = QLineEdit('0')
        self.costValue.setFixedWidth(150)
        costlayout = QHBoxLayout()
        costlayout.addWidget(costButton)
        costlayout.addWidget(self.costValue)
        costlayout.addStretch(1)
        vbox = QVBoxLayout()
        vbox.addWidget(objectiveFunction)
        vbox.addWidget(self.chi)
        vbox.addWidget(self.L1)
        vbox.addWidget(self.L2)
        vbox.addLayout(totLayout)
        vbox.addLayout(costlayout)

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
        self.eStrategy.addItems(['best1bin', 'best1exp', 'rand1exp', 'randtobest1exp', 'best2exp', 'rand2exp',
                                 'randtobest1bin', 'currenttobest1bin', 'best2bin', 'rand2bin', 'rand1bin'])
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
        self.eInit.addItems(['latinhypercube', 'sobol', 'halton', 'random'])
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
        self.shgoSampling.addItems(['simplicial', 'halton', 'sobol'])
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

        # least squares
        lsLayout = QVBoxLayout()
        self.lsWidget = QWidget()

        lsJacLayout = QHBoxLayout()
        lsJacLabel = QLabel('Jac')
        lsJacLabel.setFixedWidth(70)
        self.lsJac = QComboBox()
        self.lsJac.addItems(['2-point', '3-point', 'cs'])
        self.lsJac.currentIndexChanged.connect(self.getGOParameters)
        lsJacLayout.addWidget(lsJacLabel)
        lsJacLayout.addWidget(self.lsJac)
        lsLayout.addLayout(lsJacLayout)

        lsMethodLayout = QHBoxLayout()
        lsMethodLabel = QLabel('Method')
        lsMethodLabel.setFixedWidth(70)
        self.lsMethod = QComboBox()
        self.lsMethod.addItems(['trf', 'dogbox', 'lm'])
        self.lsMethod.currentIndexChanged.connect(self.getGOParameters)
        lsMethodLayout.addWidget(lsMethodLabel)
        lsMethodLayout.addWidget(self.lsMethod)
        lsLayout.addLayout(lsMethodLayout)

        lsFtolLayout = QHBoxLayout()
        lsFtolLabel = QLabel('ftol')
        lsFtolLabel.setFixedWidth(70)
        self.lsFtol = QLineEdit()
        self.lsFtol.textChanged.connect(self.getGOParameters)
        lsFtolLayout.addWidget(lsFtolLabel)
        lsFtolLayout.addWidget(self.lsFtol)
        lsLayout.addLayout(lsFtolLayout)

        lsXtolLayout = QHBoxLayout()
        lsXtolLabel = QLabel('xtol')
        lsXtolLabel.setFixedWidth(70)
        self.lsXtol = QLineEdit()
        self.lsXtol.textChanged.connect(self.getGOParameters)
        lsXtolLayout.addWidget(lsXtolLabel)
        lsXtolLayout.addWidget(self.lsXtol)
        lsLayout.addLayout(lsXtolLayout)

        lsGtolLayout = QHBoxLayout()
        lsGtolLabel = QLabel('gtol')
        lsGtolLabel.setFixedWidth(70)
        self.lsGtol = QLineEdit()
        self.lsGtol.textChanged.connect(self.getGOParameters)
        lsGtolLayout.addWidget(lsGtolLabel)
        lsGtolLayout.addWidget(self.lsGtol)
        lsLayout.addLayout(lsGtolLayout)

        lsXscaleLayout = QHBoxLayout()
        lsXscaleLabel = QLabel('x_scale')
        lsXscaleLabel.setFixedWidth(70)
        self.lsXscale = QLineEdit()
        self.lsXscale.textChanged.connect(self.getGOParameters)
        lsXscaleLayout.addWidget(lsXscaleLabel)
        lsXscaleLayout.addWidget(self.lsXscale)
        lsLayout.addLayout(lsXscaleLayout)

        lsLossLayout = QHBoxLayout()
        lsLossLabel = QLabel('Loss')
        lsLossLabel.setFixedWidth(70)
        self.lsLoss = QComboBox()
        self.lsLoss.addItems(['linear', 'soft_l1', 'huber', 'cauchy', 'arctan'])
        self.lsLoss.currentIndexChanged.connect(self.getGOParameters)
        lsLossLayout.addWidget(lsLossLabel)
        lsLossLayout.addWidget(self.lsLoss)
        lsLayout.addLayout(lsLossLayout)

        lsFscaleLayout = QHBoxLayout()
        lsFscaleLabel = QLabel('f_scale')
        lsFscaleLabel.setFixedWidth(70)
        self.lsFscale = QLineEdit()
        self.lsFscale.textChanged.connect(self.getGOParameters)
        lsFscaleLayout.addWidget(lsFscaleLabel)
        lsFscaleLayout.addWidget(self.lsFscale)
        lsLayout.addLayout(lsFscaleLayout)

        lsDiffLayout = QHBoxLayout()
        lsDiffLabel = QLabel('diff_step')
        lsDiffLabel.setFixedWidth(70)
        self.lsDiff = QLineEdit()
        self.lsDiff.textChanged.connect(self.getGOParameters)
        lsDiffLayout.addWidget(lsDiffLabel)
        lsDiffLayout.addWidget(self.lsDiff)
        lsLayout.addLayout(lsDiffLayout)

        lsMaxLayout = QHBoxLayout()
        lsMaxLabel = QLabel('max_nfev')
        lsMaxLabel.setFixedWidth(70)
        self.lsMax = QLineEdit()
        self.lsMax.textChanged.connect(self.getGOParameters)
        lsMaxLayout.addWidget(lsMaxLabel)
        lsMaxLayout.addWidget(self.lsMax)
        lsLayout.addLayout(lsMaxLayout)
        self.lsWidget.setLayout(lsLayout)

        # direct algorithm (required python 3.8)
        """
        dLayout = QVBoxLayout()
        self.dWidget = QWidget()

        dEpsLayout = QHBoxLayout()
        dEpsLabel = QLabel('eps')
        dEpsLabel.setFixedWidth(70)
        self.dEps = QLineEdit()
        self.dEps.textChanged.connect(self.getGOParameters)
        dEpsLayout.addWidget(dEpsLabel)
        dEpsLayout.addWidget(self.dEps)
        dLayout.addLayout(dEpsLayout)

        dMaxFunLayout = QHBoxLayout()
        dMaxFunLabel = QLabel('maxFun')
        dMaxFunLabel.setFixedWidth(70)
        self.dMaxFun = QLineEdit()
        self.dMaxFun.textChanged.connect(self.getGOParameters)
        dMaxFunLayout.addWidget(dMaxFunLabel)
        dMaxFunLayout.addWidget(self.dMaxFun)
        dLayout.addLayout(dMaxFunLayout)

        dMaxiterLayout = QHBoxLayout()
        dMaxiterLabel = QLabel('maxiter')
        dMaxiterLabel.setFixedWidth(70)
        self.dMaxiter = QLineEdit()
        self.dMaxiter.textChanged.connect(self.getGOParameters)
        dMaxiterLayout.addWidget(dMaxiterLabel)
        dMaxiterLayout.addWidget(self.dMaxiter)
        dLayout.addLayout(dMaxiterLayout)

        dLocalLayout = QHBoxLayout()
        dLocalLabel = QLabel('locally biased')
        dLocalLabel.setFixedWidth(70)
        self.dLocal = QCheckBox()
        self.dLocal.stateChanged.connect(self.getGOParameters)
        dLocalLayout.addWidget(dLocalLabel)
        dLocalLayout.addWidget(self.dLocal)
        dLayout.addLayout(dLocalLayout)

        dFminLayout = QHBoxLayout()
        dFminLabel = QLabel('fmin_rtol')
        dFminLabel.setFixedWidth(70)
        self.dFmin = QLineEdit()
        self.dFmin.textChanged.connect(self.getGOParameters)
        dFminLayout.addWidget(dFminLabel)
        dFminLayout.addWidget(self.dFmin)
        dLayout.addLayout(dFminLayout)

        dVtolLayout = QHBoxLayout()
        dVtolLabel = QLabel('volume tol.')
        dVtolLabel.setFixedWidth(70)
        self.dVtol = QLineEdit()
        self.dVtol.textChanged.connect(self.getGOParameters)
        dVtolLayout.addWidget(dVtolLabel)
        dVtolLayout.addWidget(self.dVtol)
        dLayout.addLayout(dVtolLayout)

        dLtolLayout = QHBoxLayout()
        dLtolLabel = QLabel('length tol.')
        dLtolLabel.setFixedWidth(70)
        self.dLtol = QLineEdit()
        self.dLtol.textChanged.connect(self.getGOParameters)
        dLtolLayout.addWidget(dLtolLabel)
        dLtolLayout.addWidget(self.dLtol)
        dLayout.addLayout(dLtolLayout)
        self.dWidget.setLayout(dLayout)
        """
        # adding the algorithm widgets to stacked layout
        self.goStackLayout.addWidget(self.evolutionWidget)
        self.goStackLayout.addWidget(self.shgoWidget)
        self.goStackLayout.addWidget(self.dualWidget)
        self.goStackLayout.addWidget(self.lsWidget)
        #self.goStackLayout.addWidget(self.dWidget)

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
        #self.checkscript()

    def calculateCost(self):
        """
        Purpose: calculate the cost function of the current data fitting iteration
        :param x_array: current iteration parameters
        """

        sample = self.sample
        bounds = self.rWidget.bounds
        weights = self.rWidget.weights

        scan = self.selectedScans.currentText()
        idx = int(self.selectedScans.currentIndex())
        y_scale = self.isLogWidget.currentText()
        fun = 0



        name = scan

        fun_val = 0
        xbound = bounds[idx]
        weights = weights[idx]

        background_shift = 0
        scaling_factor = 1
        data = self.rWidget.data_dict
        if 'Angle' not in data[name].keys():
            myDataScan = data[name]
            myData = myDataScan['Data']
            E = myDataScan['Energy']
            pol = myDataScan['Polarization']
            Rdat = np.array(myData[2])

            qz = np.array(myData[0])
            qz, Rsim = sample.reflectivity(E, qz, bShift=background_shift, sFactor=scaling_factor)
            Rsim = Rsim[pol]




            if y_scale == 'log(x)':
                Rsim = np.log10(Rsim)
                Rdat = np.log10(Rdat)

            elif y_scale == 'ln(x)':
                Rsim = np.log(Rsim)
                Rdat = np.log(Rdat)

            elif y_scale == 'qz^4':
                Rsim = np.multiply(Rsim, np.power(qz, 4))
                Rdat = np.multiply(Rdat, np.power(qz, 4))

            elif y_scale == 'x':
                pass



            m = 0
            for b in range(len(xbound)):
                lw = float(xbound[b][0])
                up = float(xbound[b][1])
                w = float(weights[b])

                idx = [x for x in range(len(qz)) if
                       qz[x] >= lw and qz[x] < up]  # determines index boundaries

                k = len(idx)
                m = m + k
                if len(idx) != 0:
                    if self.objective == 'Chi-Square':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2 / abs(Rsim[idx])) * w

                    elif self.objective == 'L1-Norm':
                        fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                    elif self.objective == 'L2-Norm':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2) * w

            if m != 0:
                fun = fun + fun_val / m

        else:
            myDataScan = data[name]
            myData = myDataScan['Data']
            Theta = myDataScan['Angle']
            Rdat = np.array(myData[2])
            E = np.array(myData[3])
            pol = myDataScan['Polarization']

            E, Rsim = sample.energy_scan(Theta, E)
            Rsim = Rsim[pol]




            if y_scale == 'log(x)':
                Rsim = np.log10(Rsim)
                Rdat = np.log10(Rdat)

            elif y_scale == 'ln(x)':
                Rsim = np.log(Rsim)
                Rdat = np.log(Rdat)

            elif  y_scale == 'qz^4':
                qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                Rsim = np.multiply(Rsim, np.power(qz, 4))
                Rdat = np.multiply(Rdat, np.power(qz, 4))

            elif y_scale == 'x':
                pass


            m = 0
            for b in range(len(xbound)):
                lw = float(xbound[b][0])
                up = float(xbound[b][1])
                w = float(weights[b])

                idx = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]  # determines index boundaries

                k = len(idx)
                m = m + k
                if len(idx) != 0:
                    if self.objective == 'Chi-Square':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2 / abs(Rsim[idx])) * w
                    elif self.objective == 'L1-Norm':
                        fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                    elif self.objective == 'L2-Norm':
                        fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2) * w

            if m != 0:

                fun = fun + fun_val / m


        self.costValue.setText(str(fun))
    def eventFilter(self, source, event):
        """
        Purpose: Allows user to remove fits directly from globalOptimization widget
        :param source: the source of the signal
        :param event: the type of event
        """

        if event.type() == QtCore.QEvent.MouseButtonPress:
            if event.button() == Qt.RightButton:
                top_menu = QMenu()  # initializes the menu

                menu = top_menu.addMenu("Menu")

                _remove_fit = menu.addAction('Remove Fit')

                action = menu.exec_(QtGui.QCursor.pos())
                my_rows = []
                n = len(self.sWidget.parameterFit)
                my_other_rows = []
                if action == _remove_fit:  # removes the appropriate parameter
                    my_items = self.fittingParamTable.selectedIndexes()
                    for i in my_items:
                        row = i.row()
                        if row < n:
                            if row not in my_rows:
                                my_rows.append(row)
                        else:
                            if row not in my_other_rows:
                                my_other_rows.append(row-n-1)


                my_rows = np.array(my_rows)
                for k in range(len(my_rows)):
                    self.sWidget.parameterFit.pop(my_rows[k])
                    self.sWidget.currentVal.pop(my_rows[k])
                    my_rows = my_rows - 1

                for k in range(len(my_other_rows)):
                    self.rWidget.sfBsFitParams.pop(my_other_rows[k])
                    self.rWidget.currentVal.pop(my_other_rows[k])
                    my_other_rows = my_other_rows - 1

            self.sWidget.setTable()
            self.sWidget.setTableEShift()
            self.sWidget.setTableMag()
            self.sWidget.setTableVar()
            self.setTableFit()

        return False

    def _clear_fit(self):
        """
        Purpose: Allows the user to clear the fitting parameters
        """

        # clears fitting parameters across multiple widgets
        self.sWidget.parameterFit = []
        self.rWidget.sfBsFitParams = []
        self.sWidget.currentVal = []
        self.rWidget.currentVal = []
        self.sWidget.setTable()
        self.sWidget.setTableMag()
        self.sWidget.setTableVar()
        self.setTableFit()

    def _set_y_scale(self):
        """
        Purpose: set R transformation for progress info!
        """

        self.pWidget.y_scale = self.isLogWidget.currentText()

    def _changeFitVar(self):
        """
        Purpose: change the fitting parameter value boundaries
        :return:
        """
        row = self.fittingParamTable.currentRow()  # current row
        col = self.fittingParamTable.currentColumn()  # current column

        ns = len(self.sWidget.currentVal)  # number of sampleWidget fitting parameters
        item = self.fittingParamTable.currentItem().text()  # current item

        if row <= ns - 1:  # checks if current row is from sampleWidget
            if col == 2:
                self.sWidget.currentVal[row][1][0] = float(item)
            elif col == 3:
                self.sWidget.currentVal[row][1][1] = float(item)
        else:  # fitting parameter not from sampleWidget but reflectivityWidget
            if col == 2:
                self.rWidget.currentVal[row][1][0] = float(item)
            elif col == 3:
                self.rWidget.currentVal[row][1][1] = float(item)


    def _changeObjectiveFunction(self):
        # determine which objective function to use
        rbtn = self.sender()

        if rbtn.isChecked() == True:
            self.objective = rbtn.text()

    def _changeShapeWeight(self):
        # change the weight of the total variation parameter
        value = self.totalVarWeight.text()
        if value != '':
            self.shape_weight = float(self.totalVarWeight.text())

    def _stop_optimization(self):
        # stops the optimization or data fitting algorithm from running
        global stop
        stop = True

    def changeFitParameters(self):
        """
        Purpose: Takes all the fitting parameters and save them to the new fi
        """
        # This function simply takes all the fitting parameters and saves the new fit
        for idx, fit in enumerate(self.parameters):
            if type(fit[0]) != str:  # structural, polymorphous, magnetic
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
                                    self.sWidget.structTableInfo[layer][i][1] = float(
                                        self.sWidget.structTableInfo[layer][i][1]) - diff

                        elif char == 'DENSITY':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][2])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                s = float(self.sWidget.structTableInfo[layer][i][6])
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][2] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][2] = float(
                                        self.sWidget.structTableInfo[layer][i][2]) - s * diff

                        elif char == 'ROUGHNESS':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][3])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][3] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][3] = float(
                                        self.sWidget.structTableInfo[layer][i][3]) - diff

                        elif char == 'LINKED ROUGHNESS':
                            p = float(self.sWidget.structTableInfo[layer][ele_idx][4])
                            diff = p - float(self.x[idx])
                            for i in range(len(self.sWidget.structTableInfo[layer])):
                                if i == ele_idx:  # makes sure that we are subtracting the difference value
                                    self.sWidget.structTableInfo[layer][ele_idx][4] = self.x[idx]
                                else:
                                    self.sWidget.structTableInfo[layer][i][4] = float(
                                        self.sWidget.structTableInfo[layer][i][4]) - diff

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
                        self.sWidget.varData[element][layer][1][0] = 1 - float(self.x[idx])
                    elif j == 0:
                        self.sWidget.varData[element][layer][1][1] = 1 - float(self.x[idx])

                elif my_type == 'MAGNETIC':
                    if len(fit) == 3:
                        element = fit[2]
                        self.sWidget.magData[element][layer][1][0] = self.x[idx]
                    elif len(fit) == 4:
                        element = fit[2]
                        poly = fit[3]
                        j = list(self.sWidget.magData[element][layer][0]).index(poly)
                        self.sWidget.magData[element][layer][1][j] = self.x[idx]

            else:  # scattering factor, background shift, scaling factor
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
        """
        Purpose: update optimization to sampleWidget and update boundaries in globalOptimizationWidget
        """

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

        # including scipt implementation
        script, problem, my_error = checkscript(self.sWidget.sample)
        state = self.checkBox.checkState()
        use_script = False
        if not(problem) and state > 0:
            use_script = True

        self.sWidget.sample, self.rWidget.bs, self.rWidget.sf = go.changeSampleParams(self.x, self.parameters,
                                                                                      copy.deepcopy(self.sWidget.sample),
                                                                                      self.rWidget.bs, self.rWidget.sf, script, use_script=use_script)

        # updates all the sample information across all the different Widgets
        self.rWidget.sample = copy.deepcopy(self.sWidget.sample)
        self.sample = copy.deepcopy(self.sWidget.sample)

        # update the background shift and scaling factors
        name = self.rWidget.selectedScans.currentText()
        self.rWidget.backgroundShift.blockSignals(True)
        self.rWidget.scalingFactor.blockSignals(True)
        self.rWidget.backgroundShift.setText(self.rWidget.bs[name])
        self.rWidget.scalingFactor.setText(self.rWidget.sf[name])
        self.rWidget.backgroundShift.blockSignals(False)
        self.rWidget.scalingFactor.blockSignals(False)


        self.setTableFit()
        self.sWidget.sample = copy.deepcopy(self.sample)
        self.rWidget.sample = copy.deepcopy(self.sample)
        self.sWidget._setStructFromSample(self.sample)  # required for when changing to different tab
        self.sWidget._setVarMagFromSample(self.sample)
        self.sWidget.setTable()
        self.sWidget.setTableVar()
        self.sWidget.setTableMag()
        self.sWidget.eShiftFromSample(self.sample)
        self.sWidget.setTableEShift()



    def plot_scan(self):
        """
        Purpose: plot and compare the data, previous simulation, and new fit all in one graph
        """
        script, problem, my_error = checkscript(self.sample)
        use_script = False
        check = self.checkBox.checkState()

        if not(problem) and check>0:
            use_script = True

        self.plotWidget.clear()  # clear current graph
        name = self.selectedScans.currentText()  # name of selected scan

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
            scan_type = 'Reflectivity'
            if 'Angle' in list(self.rWidget.data_dict[name].keys()):
                scan_type = 'Energy'


            step_size = float(self.sWidget._step_size)

            sample1 = copy.deepcopy(self.sample)
            sample2 = self.sample  # just to make python happy
            isGO = False

            backS = copy.deepcopy(self.rWidget.bs)
            scaleF = copy.deepcopy(self.rWidget.sf)
            if len(self.x) != 0:
                sample2, backS, scaleS = go.changeSampleParams(self.x, self.parameters, copy.deepcopy(sample2), backS,
                                                               scaleF, script, use_script=use_script)
                isGO = True

            scaling_factor_old = float(self.rWidget.sf[name])
            background_shift_old = float(self.rWidget.bs[name])

            scaling_factor = float(scaleF[name])
            background_shift = float(backS[name])

            if scan_type == 'Reflectivity':
                qz = dat[0]

                R = dat[2]
                E = self.rWidget.data_dict[name]['Energy']
                qz, Rsim = sample1.reflectivity(E, qz, s_min=step_size, sFactor=scaling_factor_old,
                                                bShift=background_shift_old)

                Theta = np.arcsin(qz / (E * 0.001013546143)) * 180 / np.pi
                Rsim = Rsim[pol]

                if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(qz, Rsim, pen=pg.mkPen((1, 3), width=2), name='Simulation')
                        if isGO:
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size, sFactor=scaling_factor,
                                                           bShift=background_shift)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(qz, Rgo, pen=pg.mkPen((2, 3), width=2), name='Optimized')

                    else:

                        self.plotWidget.plot(Theta, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 3), width=2), name='Simulation')
                        if isGO:
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size, sFactor=scaling_factor,
                                                           bShift=background_shift)
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
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size, bShift=background_shift,
                                                           sFactor=scaling_factor)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(qz[rm_idx], Rgo[rm_idx], pen=pg.mkPen((2, 3), width=2),
                                                 name='Optimized')
                    else:
                        self.plotWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 3), width=2), name='Data')
                        self.plotWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 3), width=2),
                                             name='Simulation')
                        if isGO:
                            qz, Rgo = sample2.reflectivity(E, qz, s_min=step_size, sFactor=scaling_factor,
                                                           bShift=background_shift)
                            Rgo = Rgo[pol]
                            self.plotWidget.plot(Theta[rm_idx], Rgo[rm_idx], pen=pg.mkPen((2, 3), width=2),
                                                 name='Optimized')

                    self.plotWidget.setLogMode(False, False)
                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
            elif scan_type == 'Energy':
                E = dat[3]
                R = dat[2]
                Theta = self.rWidget.data_dict[name]['Angle']

                E, Rsim = sample1.energy_scan(Theta, E, s_min=step_size, sFactor=scaling_factor_old,
                                              bShift=background_shift_old)
                if isGO:
                    qz, Rgo = sample2.energy_scan(Theta, E, s_min=step_size, sFactor=scaling_factor,
                                                  bShift=background_shift)
                    Rgo = Rgo[pol]

                Rsim = Rsim[pol]

                self.plotWidget.plot(E, R, pen=pg.mkPen((0, 3), width=2), name='Data')
                self.plotWidget.plot(E, Rsim, pen=pg.mkPen((1, 3), width=2), name='Simulation')
                if isGO:
                    qz, Rgo = sample2.energy_scan(Theta, E, s_min=step_size, sFactor=scaling_factor,
                                                  bShift=background_shift)
                    Rgo = Rgo[pol]
                    self.plotWidget.plot(E, Rgo, pen=pg.mkPen((2, 3), width=2), name='Optimized')

                self.plotWidget.setLogMode(False, False)
                self.plotWidget.setLabel('left', "Reflectivity, R")
                self.plotWidget.setLabel('bottom', "Energy, E (eV)")

    def run_first(self):
        """
        Purpose: Run this first before a data fitting process begins. Performs the appropriate initialization
        """

        # reset stop parameter
        global stop
        stop = False



        # putting the parameters and their boundaries in the proper format!
        parameters = copy.deepcopy(self.sWidget.parameterFit)
        for fit in self.rWidget.sfBsFitParams:
            parameters.append(fit)

        self.parameters = parameters  # needed for creating new sample
        lw = []
        up = []
        x0 = []

        # sorting lower boundary, upper boundary, and current value
        for b in self.sWidget.currentVal:
            lw.append(float(b[1][0]))
            up.append(float(b[1][1]))
            x0.append(float(b[0]))

        for b in self.rWidget.currentVal:
            lw.append(float(b[1][0]))
            up.append(float(b[1][1]))
            x0.append(float(b[0]))

        bounds = list(zip(lw, up))  # create boundary list

        scans = copy.deepcopy(self.rWidget.fit)  # retrieve scans to fit

        # determines the boundaries of the scans
        sBounds = []
        for bound in self.rWidget.bounds:
            temp = []
            for b in bound:
                temp.append((float(b[0]), float(b[1])))
            sBounds.append(temp)

        # retrieve boundary weights
        sWeights = []
        for weight in self.rWidget.weights:
            temp = []
            for w in weight:
                temp.append(float(w))
            sWeights.append(temp)

        data_dict = self.rWidget.data_dict

        sample = copy.deepcopy(self.sample)

        backS = copy.deepcopy(self.rWidget.bs)
        scaleF = copy.deepcopy(self.rWidget.sf)

        idx = self.algorithmSelect.currentIndex()


        if len(parameters) != 0 and len(scans) != 0:
            if idx == 0:
                # initialize the parameters for optimization saving
                self.pWidget.startSaving(sample, data_dict, scans, backS, scaleF, parameters, sBounds, sWeights,
                                         self.objective, self.shape_weight)
            elif idx == 1:
                self.pWidget.startSaving(sample, data_dict, scans, backS, scaleF, parameters, sBounds, sWeights,
                                         self.objective, self.shape_weight)
            elif idx == 2:
                self.pWidget.startSaving(sample, data_dict, scans, backS, scaleF, parameters, sBounds, sWeights,
                                         self.objective, self.shape_weight)
            elif idx == 3:
                self.pWidget.startSaving(sample, data_dict, scans, backS, scaleF, parameters, sBounds, sWeights,
                                         self.objective, self.shape_weight)

        # update sample slab class
        self.sample = copy.deepcopy(self.sWidget._createSample())
        self.temp_sample = copy.deepcopy(self.sample)

        self.runButton.setStyleSheet('background: red')
        self.runButton.blockSignals(True)

    def optimizationFinished(self):
        """
        Purpose: Peform this after optimization has finished
        """

        # send signal to callback function to stop data fitting
        global stop
        stop = True



        self.update_worker.stop()  # stop update worker
        self.update_thread.quit()  # quit update worker thread

        self.thread.wait()
        self.update_thread.wait()
        # The while loop is used to check and make sure the thread has finished before deleting
        #while not(self.update_thread.isFinished()):
        #    time.sleep(0.5)
        #    pass

        self.update_thread.deleteLater()  # delete update thread
        self.update_worker.deleteLater() # delete update worker

        self.sWidget.resetX = False  # do not reset x
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
        self.update_worker = UpdateWorker(self.pWidget)

        self.plot_scan()
        self.setTableFit()
        self.runButton.setStyleSheet('background: green')
        self.runButton.blockSignals(False)


    def _run_global_optimization(self):
        """
        Purpose: Set up threads to run data fitting and update function in parallel
        """
        # perform the proper checks first so that global optimization runs smoothly!
        empty_fit = False
        empty_scans = False
        boundary_bad = False
        weight_bad = False
        data_bad = False
        asymmetry_check = False
        do_data_fit = True


        if len(self.sWidget.parameterFit) == 0 and len(self.rWidget.sfBsFitParams) == 0:
            empty_fit = True
            do_data_fit = False

        if len(self.rWidget.fit) == 0:
            empty_scans = True
            do_data_fit = False

        for bound in self.rWidget.bounds:
            for b in bound:
                if not(isfloat(b[0])) or not(isfloat(b[1])):
                    boundary_bad = True
                    do_data_fit = False


        for weight in self.rWidget.weights:
            for w in weight:
                if not(isfloat(w)):
                    weight_bad = True
                    do_data_fit = False

        if len(list(self.rWidget.data_dict.keys())) == 0:
            data_bad = True
            do_data_fit = False



        if not(empty_scans):
            noisy_text = self.nWidget.smoothScale.currentText()
            global_text = self.isLogWidget.currentText()
            for name in self.rWidget.fit:
                polarization = self.rWidget.data_dict[name]['Polarization']
                if noisy_text == 'log(x)' or noisy_text == 'ln(x)' or global_text == 'log(x)' or global_text == 'ln(x)':
                    if polarization =='AC' or polarization == 'AL':
                        asymmetry_check = True
                        do_data_fit = False
        if do_data_fit:
            # determines if the script will be run
            state = self.checkBox.checkState()
            my_state = False
            if state > 0:
                my_state = True


            self.pWidget.check_script_state(my_state)

            # runs the optimizer method to perform the global optimization
            self.thread = QThread()  # initialize the thread
            self.update_thread = QThread()

            self.worker = Worker(self)  # data fitting worker
            self.update_worker = UpdateWorker(self.pWidget)  # progress report worker

            # move worker to thread
            self.worker.moveToThread(self.thread)  # starts two threads
            self.update_worker.moveToThread(self.update_thread)

            self.thread.started.connect(self.run_first)  # run run_first when process started
            self.thread.started.connect(self.worker.run) # then run the worker
            self.update_thread.started.connect(self.pWidget.start)
            self.update_thread.started.connect(self.update_worker.run)

            self.worker.finished.connect(self.thread.quit)  # quit
            self.worker.finished.connect(self.thread.deleteLater)
            self.worker.finished.connect(self.worker.deleteLater)
            self.worker.finished.connect(self.optimizationFinished)

            self.thread.start()
            self.update_thread.start()



        else: # these are important checks that makes sure everything is entered properly into the GUI
            if empty_fit:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Fitting Parameters")
                messageBox.setText("A fitting parameter must be selected to perform a data fit.")
                messageBox.exec()
            elif empty_scans:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Data Scan")
                messageBox.setText("A data scan must be selected to perform a data fit.")
                messageBox.exec()
            elif boundary_bad:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Scan Boundary")
                messageBox.setText("There is an error in the scan boundary input. Please check the scan boundaries in the Reflectivity Workspace.")
                messageBox.exec()
            elif weight_bad:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Scan Weight")
                messageBox.setText("There is an error in the scan boundary input. Please check the scan boundaries in the Reflectivity Workspace.")
                messageBox.exec()
            elif data_bad:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Data")
                messageBox.setText("Data must be loaded into the workspace in order to perform a data fit.")
                messageBox.exec()
            elif asymmetry_check:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Optimization and Smoothing Scale")
                messageBox.setText("A logarithmic transformation is being used on an asymmetry scan, which has potential for undefined values (logarithm of negative values). It is suggested to use 'x' as the optimization scale in the Noise Reduction and Optimization workspace.")
                messageBox.exec()

    def _optimizer(self):
        """
        Purpose: Run data fitting algorithms
        """
        script, problem, my_error = checkscript(self.sample)
        use_script = False
        if not(problem):
            use_script = True
        # getting the scans and putting them in their proper format
        # putting the parameters and their boundaries in the proper format!
        parameters = copy.deepcopy(self.sWidget.parameterFit)
        for fit in self.rWidget.sfBsFitParams:
            parameters.append(fit)

        self.parameters = parameters  # needed for creating new sample
        lw = []
        up = []
        x0 = []

        # organizing boundaries and values
        for b in self.sWidget.currentVal:
            lw.append(float(b[1][0]))
            up.append(float(b[1][1]))
            x0.append(float(b[0]))

        for b in self.rWidget.currentVal:
            lw.append(float(b[1][0]))
            up.append(float(b[1][1]))
            x0.append(float(b[0]))

        bounds = list(zip(lw, up))  # boundary list

        scans = copy.deepcopy(self.rWidget.fit)  # scans for data fitting

        # determines the boundaries of the scans
        sBounds = []
        for bound in self.rWidget.bounds:
            temp = []
            for b in bound:
                temp.append((float(b[0]), float(b[1])))
            sBounds.append(temp)

        # determine the weights of the scans
        sWeights = []
        for weight in self.rWidget.weights:
            temp = []
            for w in weight:
                temp.append(float(w))
            sWeights.append(temp)

        x = []
        fun = 0
        data_dict = self.rWidget.data_dict
        smooth_dict = copy.deepcopy(self.nWidget.smoothScans)
        data = self.rWidget.data
        sample = copy.deepcopy(self.sample)

        backS = copy.deepcopy(self.rWidget.bs)
        scaleF = copy.deepcopy(self.rWidget.sf)

        idx = self.algorithmSelect.currentIndex()

        r_scale = self.isLogWidget.currentText()  # determines what scale to use in the global optimization

        # run the selected data fitting algorithm
        if len(parameters) != 0 and len(scans) != 0:
            if idx == 0:
                x, fun = go.differential_evolution(sample, data, data_dict, scans, backS, scaleF, parameters, bounds,
                                                   sBounds, sWeights,
                                                   self.goParameters['differential evolution'], self.callback,
                                                   self.objective, self.shape_weight, r_scale, smooth_dict, script, use_script=use_script)
            elif idx == 1:
                x, fun = go.shgo(sample, data, data_dict, scans, backS, scaleF, parameters, bounds, sBounds, sWeights,
                                 self.goParameters['simplicial homology'], self.callback,
                                 self.objective, self.shape_weight, r_scale, smooth_dict, script, use_script=use_script)
            elif idx == 2:
                x, fun = go.dual_annealing(sample, data, data_dict, scans, backS, scaleF, parameters, bounds, sBounds,
                                           sWeights,
                                           self.goParameters['dual annealing'], self.callback,
                                           self.objective, self.shape_weight, r_scale, smooth_dict, script, use_script=use_script)
            elif idx == 3:
                bounds = (lw, up)
                x, fun = go.least_squares(x0, sample, data, data_dict, scans, backS, scaleF, parameters, bounds,
                                          sBounds, sWeights,
                                          self.goParameters['least squares'], self.callback,
                                          self.objective, self.shape_weight, r_scale, smooth_dict, script, use_script)
            """
            elif idx == 4:
                bounds = (lw,up)
                x,fun = go.direct(sample, data, data_dict, scans, backS, scaleF, parameters, bounds,
                                                   sBounds, sWeights,
                                                   self.goParameters['direct'], self.callback,
                                                 self.objective, self.shape_weight, r_scale, smoothe_dict, script, use_script)
            """
        else:
            print('Try again')

        return x, fun

    def getGOParameters(self):
        """
        Purpose: Retrieve the data fitting algorithm parameters from globalOptimizationWidget
        """

        # block all signals
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
        self.lsJac.blockSignals(True)
        self.lsMethod.blockSignals(True)
        self.lsFtol.blockSignals(True)
        self.lsXtol.blockSignals(True)
        self.lsGtol.blockSignals(True)
        self.lsXscale.blockSignals(True)
        self.lsLoss.blockSignals(True)
        self.lsFscale.blockSignals(True)
        self.lsDiff.blockSignals(True)
        self.lsMax.blockSignals(True)
        """
        self.dEps.blockSignals(True)
        self.dMaxFun.blockSignals(True)
        self.dMaxiter.blockSignals(True)
        self.dLocal.blockSignals(True)
        self.dFmin.blockSignals(True)
        self.dVtol.blockSignals(True)
        self.dLtol.blockSignals(True)
        """

        # retrieve the data fitting parameters depending on fitting algorithm selected by user
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
        elif idx == 3:  # least square
            self.goParameters['least squares'][0] = self.lsJac.currentText()
            self.goParameters['least squares'][1] = self.lsMethod.currentText()
            self.goParameters['least squares'][2] = self.lsFtol.text()
            self.goParameters['least squares'][3] = self.lsXtol.text()
            self.goParameters['least squares'][4] = self.lsGtol.text()
            self.goParameters['least squares'][5] = self.lsXscale.text()
            self.goParameters['least squares'][6] = self.lsLoss.currentText()
            self.goParameters['least squares'][7] = self.lsFscale.text()
            self.goParameters['least squares'][8] = self.lsDiff.text()
            self.goParameters['least squares'][9] = self.lsMax.text()
        """
        elif idx == 4: # direct algorithm
            self.goParameters['direct'][0] = self.dEps.text()
            self.goParameters['direct'][1] = self.dMaxFun.text()
            self.goParameters['direct'][2] = self.dMaxiter.text()
            if self.dLocal.checkState() == 0:
                self.goParameters['direct'][3] = 'False'
            else:
                self.goParameters['direct'][3] = 'True'
            self.goParameters['direct'][4] = self.dFmin.text()
            self.goParameters['direct'][5] = self.dVtol.text()
            self.goParameters['direct'][6] = self.dLtol.text()
        """

        # unblock all signals
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
        self.lsJac.blockSignals(False)
        self.lsMethod.blockSignals(False)
        self.lsFtol.blockSignals(False)
        self.lsXtol.blockSignals(False)
        self.lsGtol.blockSignals(False)
        self.lsXscale.blockSignals(False)
        self.lsLoss.blockSignals(False)
        self.lsFscale.blockSignals(False)
        self.lsDiff.blockSignals(False)
        self.lsMax.blockSignals(False)
        """
        self.dEps.blockSignals(False)
        self.dMaxFun.blockSignals(False)
        self.dMaxiter.blockSignals(False)
        self.dLocal.blockSignals(False)
        self.dFmin.blockSignals(False)
        self.dVtol.blockSignals(False)
        self.dLtol.blockSignals(False)
        """
        self.setGOParameters()

    def setGOParameters(self):
        """
        Purpose: set the algorithm fitting parameters
        """

        # block all signals
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
        self.lsJac.blockSignals(True)
        self.lsMethod.blockSignals(True)
        self.lsFtol.blockSignals(True)
        self.lsXtol.blockSignals(True)
        self.lsGtol.blockSignals(True)
        self.lsXscale.blockSignals(True)
        self.lsLoss.blockSignals(True)
        self.lsFscale.blockSignals(True)
        self.lsDiff.blockSignals(True)
        self.lsMax.blockSignals(True)
        """
        self.dEps.blockSignals(True)
        self.dMaxFun.blockSignals(True)
        self.dMaxiter.blockSignals(True)
        self.dLocal.blockSignals(True)
        self.dFmin.blockSignals(True)
        self.dVtol.blockSignals(True)
        self.dLtol.blockSignals(True)
        """
        idx = self.algorithmSelect.currentIndex()

        # differential evolution
        self.eStrategy.setCurrentText(self.goParameters['differential evolution'][0])
        self.eMaxiter.setText(str(self.goParameters['differential evolution'][1]))
        self.ePopsize.setText(str(self.goParameters['differential evolution'][2]))
        self.eTol.setText(str(self.goParameters['differential evolution'][3]))
        self.eAtol.setText(str(self.goParameters['differential evolution'][4]))
        self.eMinMutation.setText(str(self.goParameters['differential evolution'][5]))
        self.eMaxMutation.setText(str(self.goParameters['differential evolution'][6]))
        self.eRecomb.setText(str(self.goParameters['differential evolution'][7]))
        self.eInit.setCurrentText(str(self.goParameters['differential evolution'][9]))
        self.eUpdating.setCurrentText(str(self.goParameters['differential evolution'][10]))

        # simplicial homology
        self.shgoN.setText(str(self.goParameters['simplicial homology'][0]))
        self.shgoIter.setText(str(self.goParameters['simplicial homology'][1]))
        self.shgoSampling.setCurrentText(str(self.goParameters['simplicial homology'][2]))

        # dual annealing
        self.dualMaxiter.setText(str(self.goParameters['dual annealing'][0]))
        self.dualInitTemp.setText(str(self.goParameters['dual annealing'][1]))
        self.dualRestartTemp.setText(str(self.goParameters['dual annealing'][2]))
        self.dualVisit.setText(str(self.goParameters['dual annealing'][3]))
        self.dualAccept.setText(str(self.goParameters['dual annealing'][4]))
        self.dualMaxfun.setText(str(self.goParameters['dual annealing'][5]))

        # least squares
        self.lsJac.setCurrentText(str(self.goParameters['least squares'][0]))
        self.lsMethod.setCurrentText(str(self.goParameters['least squares'][1]))
        self.lsFtol.setText(str(self.goParameters['least squares'][2]))
        self.lsXtol.setText(str(self.goParameters['least squares'][3]))
        self.lsGtol.setText(str(self.goParameters['least squares'][4]))
        self.lsXscale.setText(str(self.goParameters['least squares'][5]))
        self.lsLoss.setCurrentText(str(self.goParameters['least squares'][6]))
        self.lsFscale.setText(str(self.goParameters['least squares'][7]))
        self.lsDiff.setText(str(self.goParameters['least squares'][8]))
        self.lsMax.setText(str(self.goParameters['least squares'][9]))

        """
        # direct
        self.dEps.setText(str(self.goParameters['direct'][0]))
        self.dMaxFun.setText(str(self.goParameters['direct'][1]))
        self.dMaxiter.setText(str(self.goParameters['direct'][2]))
        self.dFmin.setText(str(self.goParameters['direct'][4]))
        self.dVtol.setText(str(self.goParameters['direct'][5]))
        self.dLtol.setText(str(self.goParameters['direct'][6]))
        """

        # unblock all signals
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
        self.lsJac.blockSignals(False)
        self.lsMethod.blockSignals(False)
        self.lsFtol.blockSignals(False)
        self.lsXtol.blockSignals(False)
        self.lsGtol.blockSignals(False)
        self.lsXscale.blockSignals(False)
        self.lsLoss.blockSignals(False)
        self.lsFscale.blockSignals(False)
        self.lsDiff.blockSignals(False)
        self.lsMax.blockSignals(False)
        """
        self.dEps.blockSignals(False)
        self.dMaxFun.blockSignals(False)
        self.dMaxiter.blockSignals(False)
        self.dLocal.blockSignals(False)
        self.dFmin.blockSignals(False)
        self.dVtol.blockSignals(False)
        self.dLtol.blockSignals(False)
        """
    def change_algorithm(self):
        # set the proper algorithm widget
        idx = self.algorithmSelect.currentIndex()
        self.goStackLayout.setCurrentIndex(idx)

    def updateScreen(self):

        # shows only the selected scans for fitting
        self.selectedScans.clear()
        for text in self.rWidget.fit:
            self.selectedScans.addItem(text)

    def clearTableFit(self):
        # clear the fitting table
        self.fittingParamTable.clear()

    def setTableFit(self):
        """
        Purpose: Set the fitting table
        """

        # set the headers
        self.fittingParamTable.blockSignals(True)
        self.fittingParamTable.setHorizontalHeaderLabels(
            ['Name', 'Current Value', 'Lower Boundary', 'Upper Boundary', 'New'])

        # initialize the boundaries
        rows = len(self.sWidget.parameterFit) + len(self.rWidget.sfBsFitParams)  # total number of rows required
        self.fittingParamTable.setRowCount(rows)

        if self.sWidget.resetX:
            self.x = []

        # create the names and set number of rows
        row = 0
        for idx, param in enumerate(self.sWidget.parameterFit):
            name = self.getName(param)
            value = str(self.sWidget.currentVal[idx][0])
            lower = str(self.sWidget.currentVal[idx][1][0])
            upper = str(self.sWidget.currentVal[idx][1][1])
            item1 = QTableWidgetItem(name)
            item2 = QTableWidgetItem(value)

            item3 = QTableWidgetItem(lower)
            item4 = QTableWidgetItem(upper)

            self.fittingParamTable.setItem(row, 0, item1)
            self.fittingParamTable.setItem(row, 1, item2)
            self.fittingParamTable.setItem(row, 2, item3)
            self.fittingParamTable.setItem(row, 3, item4)

            row = row + 1

        # create name and set row numbers
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

        # set the table
        for idx in range(len(self.x)):
            item = QTableWidgetItem(str(self.x[idx]))
            self.fittingParamTable.setItem(idx, 4, item)

        self.fittingParamTable.blockSignals(False)

    def getName(self, p):
        """
        Purpose: create the name to display on data fitting table
        :param p: the parameter
        :return: a sting that demonstrates the parameter in an efficient way
        """
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
                    name = 'bShift-' + p[1]
            elif p[0] == 'SCALING FACTOR':  # scaling factor
                if p[1] == 'ALL SCANS':
                    name = 'sFactor-All'
                else:
                    name = 'sFactor-' + p[1]

        return name

class dataSmoothingWidget(QWidget):
    """
    Purpose: Widget used to smooth data for total variation penalty
    """
    def __init__(self):
        super(dataSmoothingWidget, self).__init__()

        self.selectedScans = []  # scans to be used in the data fitting
        self.data_dict = dict()  # used to retrieve the data
        self.boundaries = []

        # will contain the parameters and the smoothed data
        self.smoothScans = dict()  # used to store the smoothed out scans

        pagelayout = QHBoxLayout()

        optionLayout = QVBoxLayout()

        self.parameterLayout = QStackedLayout()

        # button to change the smoothing parameters
        self.setSmooth = QPushButton('Set Noise Reduction')
        self.setSmooth.setFixedWidth(200)
        self.setSmooth.clicked.connect(self._selectNoiseReduction)

        self.scanBoxLayout = QVBoxLayout()
        scanBoxLabel = QLabel('Selected Scans: ')
        scanBoxLabel.setFixedWidth(200)
        self.scanBox = QComboBox()
        self.scanBox.setFixedWidth(200)
        self.scanBox.activated.connect(self._setSmoothingVariables)
        self.scanBoxLayout.addWidget(scanBoxLabel)
        self.scanBoxLayout.addWidget(self.scanBox)

        self.smoothScaleLayout = QVBoxLayout()
        smoothScaleLabel = QLabel('Smoothing Scale: ')
        smoothScaleLabel.setFixedWidth(200)
        self.smoothScale = QComboBox()
        self.smoothScale.addItems(['log(x)','ln(x)','x','qz^4'])
        self.smoothScale.setFixedWidth(200)
        self.smoothScaleLayout.addWidget(smoothScaleLabel)
        self.smoothScaleLayout.addWidget(self.smoothScale)

        # create a versatile workspace where other data smoothing options can be chosen
        self.optionBoxLayout = QVBoxLayout()
        optionBoxLabel = QLabel('Methodology: ')
        optionBoxLabel.setFixedWidth(200)

        self.optionBox = QComboBox()
        self.optionBox.addItems(['Spline', 'Fourier Filter'])
        self.optionBox.setFixedWidth(200)
        self.optionBox.currentIndexChanged.connect(self._changeSmooth)
        self.optionBoxLayout.addWidget(optionBoxLabel)
        self.optionBoxLayout.addWidget(self.optionBox)

        # spline options
        self.splineLayout = QVBoxLayout()
        self.splineWidget = QWidget()

        splineSmoothLayout = QVBoxLayout()
        splineSmoothLabel = QLabel('Smoothing: ')
        splineSmoothLabel.setFixedWidth(200)
        splineSmoothLayout.addWidget(splineSmoothLabel)
        self.splineSmooth = QLineEdit()  # quantifies the amount of smoothing to do
        self.splineSmooth.setFixedWidth(200)
        self.splineSmooth.setText('1')
        self.splineSmooth.editingFinished.connect(self._plotGraph)
        splineSmoothLayout.addWidget(self.splineSmooth)

        splineKLayout = QVBoxLayout()
        splineKLabel = QLabel('Order, k (0-5): ')
        splineKLabel.setFixedWidth(200)
        splineKLayout.addWidget(splineKLabel)
        self.splineK = QLineEdit()  # quantifies the polynomial to use in the fitting
        self.splineK.setFixedWidth(200)
        self.splineK.setText('3')
        self.splineK.editingFinished.connect(self._plotGraph)
        splineKLayout.addWidget(self.splineK)

        self.splineLayout.addLayout(splineSmoothLayout)
        self.splineLayout.addLayout(splineKLayout)
        self.splineWidget.setLayout(self.splineLayout)

        # ----------------------- Fourier Filter Layout ---------------------------------------------------------------#
        self.fourierLayout = QVBoxLayout()
        self.fourierWidget = QWidget()

        # spline subtraction
        splineLayout = QVBoxLayout()
        splineLabel = QLabel('Spline Order (1-5): ')
        note = QLabel('*0 for no spline')
        splineLabel.setFixedWidth(200)
        note.setFixedWidth(200)
        self.splineOrder = QLineEdit('3')
        self.splineOrder.setFixedWidth(200)
        splineLayout.addWidget(splineLabel)
        splineLayout.addWidget(note)
        splineLayout.addWidget(self.splineOrder)

        # Button Layout
        buttonLayout = QVBoxLayout()
        self.splineView = QPushButton('View Spline')
        self.splineView.setFixedWidth(200)
        self.splineView.clicked.connect(self._plot_spline)
        self.splineSub = QPushButton("Spline Subtraction")
        self.splineSub.setFixedWidth(200)
        self.splineSub.clicked.connect(self._plot_spline_sub)
        buttonLayout.addWidget(self.splineView)
        buttonLayout.addWidget(self.splineSub)

        # Filter window
        filterLayout = QVBoxLayout()
        filterWindowLabel = QLabel('Window: ')
        filterWindowLabel.setFixedWidth(200)
        self.filterWindow = QLineEdit('100')
        self.filterWindow.setFixedWidth(200)
        filterLayout.addWidget(filterWindowLabel)
        filterLayout.addWidget(self.filterWindow)
        self.viewFilter = QPushButton('View Filter')
        self.viewFilter.setFixedWidth(200)
        self.viewFilter.clicked.connect(self._plot_filter)
        filterLayout.addWidget(self.viewFilter)

        # view filtered signal
        self.viewSignal = QPushButton('View Filtered Signal')
        self.viewSignal.clicked.connect(self._plotGraph)
        self.viewSignal.setFixedWidth(200)

        self.fourierLayout.addLayout(splineLayout)
        self.fourierLayout.addLayout(buttonLayout)
        self.fourierLayout.addStretch(0)
        self.fourierLayout.addLayout(filterLayout)
        self.fourierLayout.addStretch(0)
        self.fourierLayout.addWidget(self.viewSignal)

        self.fourierWidget.setLayout(self.fourierLayout)

        self.parameterLayout.addWidget(self.splineWidget)
        self.parameterLayout.addWidget(self.fourierWidget)

        optionLayout.addLayout(self.scanBoxLayout)
        optionLayout.addLayout(self.smoothScaleLayout)
        optionLayout.addStretch(1)
        optionLayout.addWidget(self.setSmooth)
        optionLayout.addLayout(self.optionBoxLayout)
        optionLayout.addStretch(1)
        optionLayout.addLayout(self.parameterLayout)

        myWidget = QWidget()
        myLayout = QHBoxLayout()
        # smoothing dict = {'Name': {'Data':[...],'Smoothing',[...]}}
        # create the plotting area
        self.graph = pg.PlotWidget()
        self.graph.setBackground('w')
        self.graph.addLegend()
        myLayout.addWidget(self.graph)
        myWidget.setLayout(myLayout)


        self.graph.setFixedSize(1000,600)

        pagelayout.addLayout(optionLayout)
        pagelayout.addWidget(myWidget)


        self.setLayout(pagelayout)

    def fourier_filter(self,x,y, order, window):

        #interpolate the data
        x_int = np.linspace(x[0],x[-1], num=len(x)*2+1)
        N = len(x_int)

        if order != 0:
            spl = UnivariateSpline(x, y, k=order)  # spline
            yspline = y-spl(x)  # spline subtraction
            f = interp1d(x,yspline)  # interpolated spline subtraction

            y_int = f(x_int)  # interpolated spline subtraction
            T = np.average(np.diff(x_int))  # spacing


            # take the fourier transform
            yf = fft(y_int)  # Fourier Transform
            xf = fftfreq(N,T)  # frequency
            xf = fftshift(xf)  # shift the x-axis
            yf = fftshift(yf)  # shift the y-axis

            max_val = xf[-1]


            # create the filter
            filter = np.zeros(N)
            for idx in range(N):
                if xf[idx]>-window and xf[idx] < window:
                    filter[idx] = 1

            yf_filtered = yf*filter # filter the signal
            yf_filtered = ifftshift(yf_filtered)  # invert the shift
            yf_filtered = ifft(yf_filtered)  # inverse fourier transform

            # back interpolation
            fb = interp1d(x_int, yf_filtered)
            ynew = fb(x)  # original dataset

            # running average
            #ynew = np.convolve(ynew,np.ones(5)/5,mode='same')
            # #invert
            ynew = ynew+spl(x)
        else:

            f = interp1d(x, y)  # interpolated spline subtraction

            y_int = f(x_int)  # interpolated spline subtraction
            T = np.average(np.diff(x_int))  # spacing

            # take the fourier transform
            yf = fft(y_int)  # Fourier Transform
            xf = fftfreq(N, T)  # frequency
            xf = fftshift(xf)  # shift the x-axis
            yf = fftshift(yf)  # shift the y-axis

            # create the filter
            filter = np.zeros(N)
            for idx in range(N):
                if xf[idx] > -window and xf[idx] < window:
                    filter[idx] = 1

            yf_filtered = yf * filter  # filter the signal
            yf_filtered = ifftshift(yf_filtered)  # invert the shift
            yf_filtered = ifft(yf_filtered)  # inverse fourier transform

            # back interpolation
            fb = interp1d(x_int, yf_filtered)
            ynew = fb(x)  # original dataset

            # running average
            #ynew = np.convolve(np.real(ynew), np.ones(5)/5, mode='same')

        last1 = np.real(ynew)[-1]
        last2 = np.real(ynew)[-2]
        ynew = np.convolve(np.real(ynew), np.ones(5)/5, mode='same')
        ynew[-1] = last1
        ynew[-2] = last2
        return ynew

    def _plot_spline(self):
        scan = self.scanBox.currentText()
        idx = self.scanBox.currentIndex()

        if scan != '' or scan != None:

            boundary = self.boundaries[idx]

            lw = float(boundary[0][0])
            up = float(boundary[-1][-1])

            type = self.smoothScans[scan]['Type']  # energy or reflectivity scan
            my_scale = self.smoothScale.currentText()

            if type == 'Energy':

                E = self.data_dict[scan]['Data'][3]  # energy array
                R = self.data_dict[scan]['Data'][2]  # relfectivity array
                idx = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]

                E = E[idx]
                R = R[idx]
                Theta = self.data_dict[scan]['Angle']

                if my_scale == 'log(x)':
                    R = np.log10(R)
                elif my_scale == 'ln(x)':
                    R = np.log(R)
                elif my_scale == 'x':
                    pass
                elif my_scale == 'qz^4':
                    qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                    R= np.multiply(R, np.power(qz, 4))

                self.splineOrder.blockSignals(True)
                order = int(self.splineOrder.text())
                self.splineOrder.setText(str(order))
                self.splineOrder.blockSignals(False)
                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Order out of range")
                    messageBox.setText("Order must be found between 1 and 5. If order selected to be 0 then no spline can be displayed.")
                    messageBox.exec()
                else:
                    if order == 0:
                        spl = UnivariateSpline(E, R)
                    else:
                        spl = UnivariateSpline(E, R, k=order)
                    self.graph.clear()
                    self.graph.plot()
                    self.graph.plot(E, R, pen=pg.mkPen((1, 4), width=2), name='Data')
                    self.graph.plot(E, spl(E), pen=pg.mkPen((2, 4), width=2), name='Spline')
                    self.graph.setLabel('left', "f(R)")
                    self.graph.setLabel('bottom', "Energy, E (eV)")

            elif type == 'Reflectivity':
                qz = self.data_dict[scan]['Data'][0]  # energy array
                R = self.data_dict[scan]['Data'][2]  # relfectivity array
                idx = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] < up]

                qz = qz[idx]
                R = R[idx]


                if my_scale == 'log(x)':
                    R = np.log10(R)
                elif my_scale == 'ln(x)':
                    R = np.log(R)
                elif my_scale == 'x':
                    pass
                elif my_scale == 'qz^4':
                    R = np.multiply(R, np.power(qz, 4))

                self.splineOrder.blockSignals(True)
                order = int(self.splineOrder.text())
                self.splineOrder.setText(str(order))
                self.splineOrder.blockSignals(False)
                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Order out of range")
                    messageBox.setText(
                        "Order must be found between 1 and 5. If order selected to be 0 then no spline can be displayed.")
                    messageBox.exec()
                else:
                    if order == 0:
                        spl = UnivariateSpline(qz, R)
                    else:
                        spl = UnivariateSpline(qz, R, k=order)

                    self.graph.clear()
                    self.graph.plot()
                    self.graph.plot(qz, R, pen=pg.mkPen((1, 4), width=2), name='Data')
                    self.graph.plot(qz, spl(qz), pen=pg.mkPen((2, 4), width=2), name='Spline')
                    self.graph.setLabel('left', "f(R)")
                    self.graph.setLabel('bottom', "Energy, E (eV)")

    def _plot_spline_sub(self):
        scan = self.scanBox.currentText()
        idx = self.scanBox.currentIndex()
        boundary = self.boundaries[idx]

        lw = float(boundary[0][0])
        up = float(boundary[-1][-1])
        if scan != '' or scan != None:
            type = self.smoothScans[scan]['Type']  # energy or reflectivity scan
            my_scale = self.smoothScale.currentText()

            if type == 'Energy':
                E = self.data_dict[scan]['Data'][3]  # energy array
                R = self.data_dict[scan]['Data'][2]  # relfectivity array
                idx_new = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]

                E = E[idx_new]
                R = R[idx_new]
                Theta = self.data_dict[scan]['Angle']

                if my_scale == 'log(x)':
                    R = np.log10(R)
                elif my_scale == 'ln(x)':
                    R = np.log(R)
                elif my_scale == 'x':
                    pass
                elif my_scale == 'qz^4':
                    qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                    R = np.multiply(R, np.power(qz, 4))

                self.splineOrder.blockSignals(True)
                order = int(self.splineOrder.text())
                self.splineOrder.setText(str(order))
                self.splineOrder.blockSignals(False)
                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Order out of range")
                    messageBox.setText(
                        "Order must be found between 1 and 5. If order selected to be 0 then no spline can be displayed.")
                    messageBox.exec()
                else:
                    if order == 0:
                        Rplot = R
                    else:
                        spl = UnivariateSpline(E, R, k=order)
                        Rplot = R-spl(E)

                    self.graph.clear()
                    self.graph.plot()
                    self.graph.plot(E, Rplot, pen=pg.mkPen((2, 4), width=2), name='Spline Subtraction')
                    self.graph.setLabel('left', "f(R)")
                    self.graph.setLabel('bottom', "Momentum Transfer, qz (1/angstrom)")

            elif type == 'Reflectivity':
                qz = self.data_dict[scan]['Data'][0]  # energy array
                R = self.data_dict[scan]['Data'][2]  # relfectivity array

                idx_new = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] < up]
                qz = qz[idx_new]
                R = R[idx_new]

                if my_scale == 'log(x)':
                    R = np.log10(R)
                elif my_scale == 'ln(x)':
                    R = np.log(R)
                elif my_scale == 'x':
                    pass
                elif my_scale == 'qz^4':
                    R = np.multiply(R, np.power(qz, 4))

                self.splineOrder.blockSignals(True)
                order = int(self.splineOrder.text())
                self.splineOrder.setText(str(order))
                self.splineOrder.blockSignals(False)
                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Order out of range")
                    messageBox.setText(
                        "Order must be found between 1 and 5. If order selected to be 0 then no spline can be displayed.")
                    messageBox.exec()
                else:
                    if order == 0:
                        Rplot = R
                    else:
                        spl = UnivariateSpline(qz, R, k=order)
                        Rplot = R - spl(qz)
                    self.graph.clear()
                    self.graph.plot()
                    self.graph.plot(qz, Rplot, pen=pg.mkPen((2, 4), width=2), name='Spline Subtraction')
                    self.graph.setLabel('left', "f(R)")
                    self.graph.setLabel('bottom', "Momentum Transfer, qz (1/angstrom)")

    def _plot_filter(self):
        scan = self.scanBox.currentText()
        idx = self.scanBox.currentIndex()
        boundary = self.boundaries[idx]

        lw = float(boundary[0][0])
        up = float(boundary[-1][-1])
        if scan != '' or scan != None:
            type = self.smoothScans[scan]['Type']  # energy or reflectivity scan
            my_scale = self.smoothScale.currentText()

            if type == 'Energy':
                E = self.data_dict[scan]['Data'][3]  # energy array
                R = self.data_dict[scan]['Data'][2]  # relfectivity array
                idx_new = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]

                E = E[idx_new]
                R = R[idx_new]
                Theta = self.data_dict[scan]['Angle']

                if my_scale == 'log(x)':
                    R = np.log10(R)
                elif my_scale == 'ln(x)':
                    R = np.log(R)
                elif my_scale == 'x':
                    pass
                elif my_scale == 'qz^4':
                    qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                    R = np.multiply(R, np.power(qz, 4))

                self.splineOrder.blockSignals(True)
                order = int(self.splineOrder.text())
                self.splineOrder.setText(str(order))
                self.splineOrder.blockSignals(False)
                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Order out of range")
                    messageBox.setText(
                        "Order must be found between 1 and 5. If order selected to be 0 then no spline can be displayed.")
                    messageBox.exec()
                else:
                    spl = UnivariateSpline(E, R, k=order)
                    Rspline = R-spl(E)

                    # for equal spacing
                    Emin = E[0]
                    Emax = E[-1]
                    Eint = np.linspace(Emin, Emax, num=len(R)*2)
                    f = interp1d(E, Rspline)

                    N = len(Eint)
                    T = np.average(np.diff(Eint))

                    Rf = fft(f(Eint))
                    Ef = fftfreq(N,T)
                    Ef = fftshift(Ef)
                    Rplot = fftshift(Rf)

                    # create filter
                    filter = np.zeros(N)

                    self.filterWindow.blockSignals(True)
                    window = abs(float(self.filterWindow.text()))
                    self.filterWindow.blockSignals(False)


                    for idx in range(len(filter)):
                        if Ef[idx] > -window and Ef[idx] < window:
                            filter[idx] = 1
                    max_val = max(1.0/N*np.abs(Rplot))/2
                    self.graph.clear()
                    self.graph.plot(Ef, 1.0/N*np.abs(Rplot), pen=pg.mkPen((1, 4), width=2), name='Fourier Transform')
                    self.graph.plot(Ef, max_val*filter, pen=pg.mkPen((2, 4), width=2), name='Data')
                    self.graph.setLabel('left', "f(R)")
                    self.graph.setLabel('bottom', "Time, t")

            elif type == 'Reflectivity':
                qz = self.data_dict[scan]['Data'][0]  # energy array
                R = self.data_dict[scan]['Data'][2]  # relfectivity array
                idx_new = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] < up]

                qz = qz[idx_new]
                R = R[idx_new]

                if my_scale == 'log(x)':
                    R = np.log10(R)
                elif my_scale == 'ln(x)':
                    R = np.log(R)
                elif my_scale == 'x':
                    pass
                elif my_scale == 'qz^4':
                    R = np.multiply(R, np.power(qz, 4))

                self.splineOrder.blockSignals(True)
                order = int(self.splineOrder.text())
                self.splineOrder.setText(str(order))
                self.splineOrder.blockSignals(False)
                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Order out of range")
                    messageBox.setText(
                        "Order must be found between 1 and 5. If order selected to be 0 then no spline can be displayed.")
                    messageBox.exec()
                else:

                    # for equal spacing
                    qmin = qz[0]
                    qmax = qz[-1]
                    qz_int = np.linspace(qmin, qmax, num=len(R) * 2)
                    if order == 0:
                        f = interp1d(qz, R)
                    else:
                        spl = UnivariateSpline(qz, R, k=order)
                        Rspline = R - spl(qz)
                        f = interp1d(qz, Rspline)

                    N = len(qz_int)
                    T = np.average(np.diff(qz_int))

                    Rf = fft(f(qz_int))
                    qzf = fftfreq(N, T)
                    qzf = fftshift(qzf)
                    Rplot = fftshift(Rf)


                    # create filter
                    filter = np.zeros(N)

                    self.filterWindow.blockSignals(True)
                    window = abs(float(self.filterWindow.text()))
                    self.filterWindow.blockSignals(False)

                    for idx in range(len(filter)):
                        if qzf[idx] > -window and qzf[idx] < window:
                            filter[idx] = 1

                    max_val = max(1.0 / N * np.abs(Rplot)) / 2
                    self.graph.clear()
                    self.graph.plot(qzf, 1.0/N*np.abs(Rplot), pen=pg.mkPen((1, 4), width=2), name='Fourier Transform')
                    self.graph.plot(qzf, max_val*filter, pen=pg.mkPen((2, 4), width=2), name='Spline')
                    self.graph.setLabel('left', "f(R)")
                    self.graph.setLabel('bottom', "Spatial Frequency")



    def _plotGraph(self):
        """
        Purpose: plot the data, previous smoothing iteration, and current smoothing iteration
        """
        self.graph.clear()

        scan = self.scanBox.currentText()  # scan to plot
        smooth = self.optionBox.currentText()  # which methodology to use
        type = self.smoothScans[scan]['Type']
        my_scale = self.smoothScale.currentText()
        idx = self.scanBox.currentIndex()
        boundary = self.boundaries[idx]

        lw = float(boundary[0][0])
        up = float(boundary[-1][-1])

        if type == 'Energy':
            # get data
            E = self.data_dict[scan]['Data'][3]
            idx_new = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]
            E = E[idx_new]
            Theta = self.data_dict[scan]['Angle']
            Rdata = self.data_dict[scan]['Data'][2][idx_new]  # data
            Rprev = self.smoothScans[scan]['Data'][2][idx_new] # previous smooth
            if my_scale == 'log(x)':
                Rdata = np.log10(Rdata)
                Rprev = np.log10(Rprev)
            elif my_scale == 'ln(x)':
                Rdata = np.log(Rdata)
                Rprev = np.log(Rprev)
            elif my_scale == 'x':
                pass
            elif my_scale == 'qz^4':
                qz = np.sin(Theta*np.pi/180)*(E * 0.001013546143)
                Rdata = np.multiply(Rdata,np.power(qz,4))
                Rprev = np.multilpy(Rprev, np.power(qz, 4))

            Rsmooth = copy.copy(Rdata)
            # calculate current smooth
            if smooth == 'Spline':
                s = float(self.splineSmooth.text())
                k = self.splineK.text()
                if '.' in k:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Assumption Made")
                    messageBox.setText("Program assumed integer value as polynomial order must be an integer.")
                    messageBox.exec()
                k = int(k)

                if k < 0 or k>5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Invalid Entry")
                    messageBox.setText("The polynomial order must be found between 0 and 5. Values not saved. Assumed k = 3.")
                    messageBox.exec()
                    k = 3


                Rsmooth = self._noiseRemoval(E, Rdata, s=s,k=k)
            elif smooth == "Fourier Filter":
                self.splineOrder.blockSignals(True)
                self.filterWindow.blockSignals(True)

                order = self.splineOrder.text()
                window = float(self.filterWindow.text())

                self.splineOrder.blockSignals(False)
                self.filterWindow.blockSignals(False)

                if '.' in order:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Assumption Made")
                    messageBox.setText("Program assumed integer value as spline order must be an integer.")
                    messageBox.exec()
                order = int(order)

                if order < 0 or order>5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Invalid Entry")
                    messageBox.setText("The polynomial order must be found between 0 and 5. Values not saved. Assumed k = 3.")
                    messageBox.exec()
                    order = 3

                Rsmooth = self.fourier_filter(E, Rdata, order, window)

            self.graph.plot(E, Rdata, pen=pg.mkPen((0,4), width=2), name='Data')
            self.graph.plot(E, Rprev, pen=pg.mkPen((3, 4), width=2), name='Previous')
            self.graph.plot(E, Rsmooth, pen=pg.mkPen((1, 4), width=2), name='Current')



        elif type == 'Reflectivity':
            # get data
            qz = self.data_dict[scan]['Data'][0]
            idx_new = [x for x in range(len(qz)) if qz[x] >= lw and qz[x] < up]
            qz = qz[idx_new]
            Rdata = self.data_dict[scan]['Data'][2][idx_new]  # data
            Rprev = self.smoothScans[scan]['Data'][2][idx_new]  # previous smooth
            if my_scale == 'log(x)':
                Rdata = np.log10(Rdata)
                Rprev = np.log10(Rprev)
            elif my_scale == 'ln(x)':
                Rdata = np.log(Rdata)
                Rprev = np.log(Rprev)
            elif my_scale == 'x':
                pass
            elif my_scale == 'qz^4':
                Rdata = Rdata*np.power(qz, 4)
                Rprev = Rprev*np.power(qz, 4)

            Rsmooth = copy.copy(Rdata)
            # calculate current smooth
            if smooth == 'Spline':
                s = float(self.splineSmooth.text())
                k = self.splineK.text()
                if '.' in k:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Assumption Made")
                    messageBox.setText("Program assumed integer value as polynomial order must be an integer.")
                    messageBox.exec()
                k = int(k)

                if k < 0 or k > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Invalid Entry")
                    messageBox.setText(
                        "The polynomial order must be found between 0 and 5. Values not saved. Assumed k = 3.")
                    messageBox.exec()
                    k = 3

                Rsmooth = self._noiseRemoval(qz, Rdata, s=s, k=k)

            elif smooth == "Fourier Filter":
                self.splineOrder.blockSignals(True)
                self.filterWindow.blockSignals(True)

                order = self.splineOrder.text()
                window = float(self.filterWindow.text())

                self.splineOrder.blockSignals(False)
                self.filterWindow.blockSignals(False)

                if '.' in order:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Assumption Made")
                    messageBox.setText("Program assumed integer value as spline order must be an integer.")
                    messageBox.exec()
                order = int(order)

                if order < 0 or order > 5:
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Invalid Entry")
                    messageBox.setText(
                        "The polynomial order must be found between 0 and 5. Values not saved. Assumed k = 3.")
                    messageBox.exec()
                    order = 3

                Rsmooth = self.fourier_filter(qz, Rdata, order, window)

            self.graph.plot(qz, Rdata, pen=pg.mkPen((0, 4), width=2), name='Data')
            self.graph.plot(qz, Rprev, pen=pg.mkPen((3, 4), width=2), name='Previous')
            self.graph.plot(qz, Rsmooth, pen=pg.mkPen((1, 4), width=2), name='Current')




    def _setSmoothingVariables(self):
        """
        Purpose: set the smoothing variables depending on the data smoothing methodology used
        """
        smooth = self.optionBox.currentText()  # get the smoothing implementation to use
        scan = self.scanBox.currentText()

        if smooth != '' and scan != '':
            if smooth == 'Spline':  # spline case

                s = self.smoothScans[scan][smooth][0]  # smoothing variable
                k = self.smoothScans[scan][smooth][1]  # kth-order

                # block spline signals
                self.splineSmooth.blockSignals(True)
                self.splineK.blockSignals(True)

                # setting saved smoothing and order variables
                self.splineSmooth.setText(str(s))
                self.splineK.setText(str(k))

                # unblock signals
                self.splineSmooth.blockSignals(False)
                self.splineK.blockSignals(False)

                self._plotGraph()
            elif smooth == 'Fourier Filter':
                order = self.smoothScans[scan][smooth][0]
                window = self.smoothScans[scan][smooth][1]

                self.splineOrder.blockSignals(True)
                self.filterWindow.blockSignals(True)

                self.splineOrder.setText(str(order))
                self.filterWindow.setText(str(window))

                self.splineOrder.blockSignals(True)
                self.filterWindow.blockSignals(True)


    def _selectNoiseReduction(self):
        """
        Purpose: Changes work space to smoothing methodology and checks user input
        """
        # sets the current variables being used for the selected scan
        smooth = self.optionBox.currentText()  # get the smoothing implementation to use
        scan = self.scanBox.currentText()

        if smooth != '' and scan != '':  # makes sure a methodology is chosen
            if smooth == 'Spline':  # spline case
                self.splineSmooth.blockSignals(True)
                self.splineK.blockSignals(True)
                s = float(self.splineSmooth.text())
                k = self.splineK.text()
                self.splineSmooth.blockSignals(False)
                self.splineK.blockSignals(False)

                if '.' in k:  # incorrect input
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Assumption Made")
                    messageBox.setText("Program assumed integer value as polynomial order must be an integer.")
                    messageBox.exec()
                k = int(k)

                if k < 0 or k>5:  # incorrect input
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Invalid Entry")
                    messageBox.setText("The polynomial order must be found between 0 and 5. Values not saved.")
                    messageBox.exec()
                else:
                    self.smoothScans[scan][smooth][0] = s
                    self.smoothScans[scan][smooth][1] = k

                # performs the data smoothing
                type = self.smoothScans[scan]['Type']
                if type == 'Energy':
                    E = self.data_dict[scan]['Data'][3]
                    R = copy.copy(self.data_dict[scan]['Data'][2])
                    Theta = self.data_dict[scan]['Angle']

                    # performs transformation of reflectivity for data smoothing
                    my_scale = self.smoothScale.currentText()
                    if my_scale == 'log(x)':
                        R = np.log10(R)
                    elif my_scale == 'ln(x)':
                        R = np.log(R)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                        R = np.multiply(R, np.power(qz,4))


                    Rsmooth = self._noiseRemoval(E,R, s=s,k=k)

                    # transform back to original R-scale
                    if my_scale == 'log(x)':
                        Rsmooth = np.power(10, Rsmooth)
                    elif my_scale == 'ln(x)':
                        Rsmooth = np.exp(Rsmooth)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                        Rsmooth = np.divide(Rsmooth, np.power(qz,4))


                    self.smoothScans[scan]['Data'][2] = copy.copy(Rsmooth)  # save the smoothed data

                elif type == 'Reflectivity':
                    qz = self.data_dict[scan]['Data'][0]
                    R = copy.copy(self.data_dict[scan]['Data'][2])

                    # transform reflectivity data for smoothing
                    my_scale  = self.smoothScale.currentText()
                    if my_scale == 'log(x)':
                        R = np.log10(R)
                    elif my_scale == 'ln(x)':
                        R = np.log(R)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        R = np.multiply(R, np.power(qz,4))

                    Rsmooth = self._noiseRemoval(qz, R,s=s,k=k)

                    #transform back to orginal R-scale
                    if my_scale == 'log(x)':
                        Rsmooth = np.power(10, Rsmooth)
                    elif my_scale == 'ln(x)':
                        Rsmooth = np.exp(Rsmooth)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        Rsmooth = np.divide(Rsmooth, np.power(qz,4))

                    self.smoothScans[scan]['Data'][2] = copy.copy(Rsmooth)

                self._plotGraph()

            elif smooth == "Fourier Filter":
                self.splineOrder.blockSignals(True)
                self.filterWindow.blockSignals(True)
                self.splineK.blockSignals(True)
                order = self.splineOrder.text()
                window = float(self.filterWindow.text())
                self.splineOrder.blockSignals(True)
                self.filterWindow.blockSignals(True)

                if '.' in order:  # incorrect input
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Assumption Made")
                    messageBox.setText("Program assumed integer value as polynomial order must be an integer.")
                    messageBox.exec()
                order = int(order)

                if order < 0 or order > 5:  # incorrect input
                    messageBox = QMessageBox()
                    messageBox.setWindowTitle("Invalid Entry")
                    messageBox.setText("The polynomial order must be found between 0 and 5. Values not saved.")
                    messageBox.exec()
                else:
                    self.smoothScans[scan][smooth][0] = order
                    self.smoothScans[scan][smooth][1] = window

                # performs the data smoothing
                type = self.smoothScans[scan]['Type']
                if type == 'Energy':
                    E = self.data_dict[scan]['Data'][3]
                    R = copy.copy(self.data_dict[scan]['Data'][2])
                    Theta = self.data_dict[scan]['Angle']

                    # performs transformation of reflectivity for data smoothing
                    my_scale = self.smoothScale.currentText()
                    if my_scale == 'log(x)':
                        R = np.log10(R)
                    elif my_scale == 'ln(x)':
                        R = np.log(R)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                        R = np.multiply(R, np.power(qz, 4))

                    Rsmooth = self.fourier_filter(qz,R,order,window)

                    # transform back to original R-scale
                    if my_scale == 'log(x)':
                        Rsmooth = np.power(10, Rsmooth)
                    elif my_scale == 'ln(x)':
                        Rsmooth = np.exp(Rsmooth)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                        Rsmooth = np.divide(Rsmooth, np.power(qz, 4))

                    self.smoothScans[scan]['Data'][2] = copy.copy(Rsmooth)  # save the smoothed data

                elif type == 'Reflectivity':
                    qz = self.data_dict[scan]['Data'][0]
                    R = copy.copy(self.data_dict[scan]['Data'][2])

                    # transform reflectivity data for smoothing
                    my_scale = self.smoothScale.currentText()
                    if my_scale == 'log(x)':
                        R = np.log10(R)
                    elif my_scale == 'ln(x)':
                        R = np.log(R)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        R = np.multiply(R, np.power(qz, 4))

                    Rsmooth = self.fourier_filter(qz,R,order,window)

                    # transform back to orginal R-scale
                    if my_scale == 'log(x)':
                        Rsmooth = np.power(10, Rsmooth)
                    elif my_scale == 'ln(x)':
                        Rsmooth = np.exp(Rsmooth)
                    elif my_scale == 'x':
                        pass
                    elif my_scale == 'qz^4':
                        Rsmooth = np.divide(Rsmooth, np.power(qz, 4))

                    self.smoothScans[scan]['Data'][2] = copy.copy(Rsmooth)

                self._plotGraph()

    def _changeSmooth(self):
        # changes the smoothing algorithm that we will use
        idx = self.optionBox.currentIndex()
        self.parameterLayout.setCurrentIndex(idx)
        if idx == 0:
            self._plotGraph()  # now we are using a completely different smoothing methodology



    def _resetVariables(self,data_dict, fit, boundaries):
        # This will be used in the loading in of data as well as when this tab is activated
        self.data_dict = data_dict
        self.boundaries = boundaries
        # double checking to make sure that the saved fit is found in the dataset provided
        self.selectedScans = [scan for scan in fit if scan in list(self.data_dict.keys())]

        for name in list(self.data_dict.keys()):
            self.smoothScans[name] = dict()
            self.smoothScans[name]['Data'] = copy.copy(self.data_dict[name]['Data'])
            if 'Angle' in list(self.data_dict[name].keys()):
                self.smoothScans[name]['Type'] = 'Energy'
            else:
                self.smoothScans[name]['Type'] = 'Reflectivity'
            self.smoothScans[name]['Spline'] = [1, 3]  # [s, k]
            self.smoothScans[name]['Fourier Filter'] = [3,100]

        # blocks all signals
        self.splineK.blockSignals(True)
        self.splineSmooth.blockSignals(True)
        self.setSmooth.blockSignals(True)
        self.optionBox.blockSignals(True)
        self.scanBox.blockSignals(True)

        # resets selected scans
        self.scanBox.clear()
        self.scanBox.addItems(self.selectedScans)

        # unblocks all signals
        self.splineK.blockSignals(False)
        self.splineSmooth.blockSignals(False)
        self.setSmooth.blockSignals(False)
        self.optionBox.blockSignals(False)
        self.scanBox.blockSignals(False)

        if len(self.selectedScans) != 0:
            if len(self.selectedScans) != 1 and self.selectedScans[0] != '':
                self.scanBox.setCurrentIndex(0)

    def _noiseRemoval(self,x, R, s=1, k=3):
        """
        Purpose: smoothes data using spline interpolation
        :param x: qz or E numpy array
        :param R: reflectivity numpy array
        :param s: smoothing variable
        :param k: order (0-5)
        :return:
        """
        tck = interpolate.splrep(x, R, s=s, k=k)
        return interpolate.splev(x, tck, der=0)

class ReflectometryApp(QMainWindow):
    """
    Purpose: Main Window
    """
    def __init__(self):
        super().__init__()
        cwd = os.getcwd()

        self.version = '0.2'
        self.fname = cwd + '\demo.h5'  # initial sample
        self.data = []  # data info
        self.data_dict = dict()  # dictionary that contains data
        self.sim_dict = dict()  # dictionary that contains simulation

        self.sample = ms.slab(1)  # app is initialized and no project is selected
        self.sample.addlayer(0, 'SrTiO3', 50)
        self.sample.energy_shift()

        # set the title
        my_name = 'GO-RXR (version '+self.version +')'
        self.setWindowTitle(my_name)

        # set the geometry of the window
        self.setGeometry(180, 60, 1400, 800)

        pagelayout = QVBoxLayout()  # page layout
        buttonlayout = QHBoxLayout()
        self.stackedlayout = QStackedLayout()

        pagelayout.addLayout(buttonlayout)
        pagelayout.addLayout(self.stackedlayout)

        afont = QtGui.QFont('Ariel', 11)
        afont.setBold(True)


        # initializing workspace buttons
        self.sampleButton = QPushButton('Sample')
        self.sampleButton.setFont(afont)
        self.reflButton = QPushButton('Reflectivity')
        self.reflButton.setFont(afont)
        self.smoothButton = QPushButton('Smooth Data')
        self.smoothButton.setFont(afont)
        self.goButton = QPushButton('Optimization')
        self.goButton.setFont(afont)
        self.progressButton = QPushButton('Progress')
        self.progressButton.setFont(afont)
        self.scanProgress = QComboBox()

        # initializing workspace widgets
        self._sampleWidget = sampleWidget(self.sample)  # initialize the sample widget
        self._reflectivityWidget = reflectivityWidget(self._sampleWidget, self.data, self.data_dict, self.sim_dict)
        self._noiseWidget = dataSmoothingWidget()
        self._progressWidget = progressWidget(self._sampleWidget, self._reflectivityWidget, self._noiseWidget)
        self._goWidget = GlobalOptimizationWidget(self._sampleWidget, self._reflectivityWidget, self._noiseWidget, self._progressWidget,
                                                  self)

        # initializing workspace button signals and layout
        self.sampleButton.setStyleSheet("background-color : magenta")
        self.sampleButton.clicked.connect(self.activate_tab_1)
        buttonlayout.addWidget(self.sampleButton)
        self.stackedlayout.addWidget(self._sampleWidget)

        self.reflButton.setStyleSheet("background-color : pink")
        self.reflButton.clicked.connect(self.activate_tab_2)
        buttonlayout.addWidget(self.reflButton)
        self.stackedlayout.addWidget(self._reflectivityWidget)

        self.smoothButton.setStyleSheet("background-color : pink")
        self.smoothButton.clicked.connect(self.activate_tab_3)
        buttonlayout.addWidget(self.smoothButton)
        self.stackedlayout.addWidget(self._noiseWidget)

        self.goButton.setStyleSheet("background-color : pink")
        self.goButton.clicked.connect(self.activate_tab_4)
        buttonlayout.addWidget(self.goButton)
        self.stackedlayout.addWidget(self._goWidget)

        self.progressButton.setStyleSheet('background: pink')
        self.progressButton.clicked.connect(self.activate_tab_5)
        buttonlayout.addWidget(self.progressButton)
        self.stackedlayout.addWidget(self._progressWidget)

        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)

        # starting my menu bar
        menuBar = self.menuBar()

        # saving and loading area
        fileMenu = QMenu("&File", self)

        # create a new file
        self.newFile = QAction("&New Workspace", self)
        self.newFile.triggered.connect(self._newFile)
        fileMenu.addAction(self.newFile)
        fileMenu.addSeparator()

        # load a file
        self.loadFile = QAction("&Load Workspace", self)
        self.loadFile.triggered.connect(self._loadFile)
        fileMenu.addAction(self.loadFile)

        # load a file
        self.loadSample= QAction("&Load Sample", self)
        self.loadSample.triggered.connect(self._loadSample)
        fileMenu.addAction(self.loadSample)
        fileMenu.addSeparator()

        # save current working file
        self.saveFile = QAction("&Save Workspace", self)
        self.saveFile.triggered.connect(self._saveFile)
        fileMenu.addAction(self.saveFile)

        # save file as a new project
        self.saveAsFile = QAction("&Save Workspace As", self)
        self.saveAsFile.triggered.connect(self._saveAsFile)
        fileMenu.addAction(self.saveAsFile)

        # save only the sample information
        self.saveSampleFile = QAction("&Save Sample", self)
        self.saveSampleFile.triggered.connect(self._saveSample)
        fileMenu.addAction(self.saveSampleFile)

        # save the simulation
        self.saveSimulationFile = QAction("&Save Simulation", self)
        self.saveSimulationFile.triggered.connect(self._saveSimulation)
        fileMenu.addAction(self.saveSimulationFile)
        fileMenu.addSeparator()
        self.importDataset = QAction("&Import Dataset", self)
        self.importDataset.triggered.connect(self._importDataSet)
        fileMenu.addAction(self.importDataset)

        # load a ReMagX file (only the data)
        self.loadReMagX = QAction("&Load ReMagX", self)
        self.loadReMagX.triggered.connect(self._loadReMagX)
        fileMenu.addAction(self.loadReMagX)

        # save summary of work as a textfile

        self.saveSummary = QAction("&Save Summary", self)
        self.saveSummary.triggered.connect(self._summary)
        fileMenu.addAction(self.saveSummary)
        fileMenu.addSeparator()
        # exit the application
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

        self.about = QAction('About', self)
        self.about.triggered.connect(self._help)
        helpMenu.addAction(self.about)



    def _newFile(self):
        """
        Purpose: Create a new file
        """

        # create a new file or project workspace
        filename, _ = QFileDialog.getSaveFileName()

        fname = filename.split('/')[-1]
        self.fname = filename
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


            # random sample input
            sample = ms.slab(2)
            sample.addlayer(0, 'SrTiO3', 50)
            sample.addlayer(1, 'LaMnO3', 10)
            sample.energy_shift()

            ds.newFileHDF5(self.fname, sample, self.version)

            self.sample = sample
            self._sampleWidget.sample = sample

            # reset data and simulation
            self.data = list()
            self.data_dict = dict()
            self.sim_dict = dict()

            # loading in the background shift and scaling factor
            self._reflectivityWidget.bs = dict()
            self._reflectivityWidget.sf = dict()

            # save sample information to application
            self._sampleWidget._setStructFromSample(sample)
            self._sampleWidget._setVarMagFromSample(sample)

            self._sampleWidget.eShift = dict()
            self._sampleWidget.ffScale = dict()

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
                self._sampleWidget.ffScale[name] = sample.ff_scale[key]

            for key in list(sample.mag_eShift.keys()):
                name = 'ffm-' + key
                self._sampleWidget.eShift[name] = sample.mag_eShift[key]
                self._sampleWidget.ffScale[name] = sample.ffm_scale[key]

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

            self.goParameters = {
                'differential evolution': ['currenttobest1bin', 2, 15, 1e-6, 0, 0.5, 1, 0.7, True, 'latinhypercube',
                                           'immediate'],
                'simplicial homology': ['None', 1, 'simplicial'],
                'dual annealing': [150, 5230.0, 2e-5, 2.62, 5.0, 10000000.0, True],
                'least squares': ['2-point', 'trf', 1e-8, 1e-8, 1e-8, 1.0, 'linear', 1.0, 'None', 'None']}

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
            for param in self._reflectivityWidget.sfBsFitParams:
                self._goWidget.parameters.append(param)

            # reset all of the tables!!!
            # - otherwise we have everything for the load function
        else:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Invalid file name")
            messageBox.setText("Selected file name or path is not valid. Please select a valid file name.")
            messageBox.exec()

        self.activate_tab_1()

    def _loadFile(self):
        """
        Purpose: Load a new file or project workspace
        """

        self.fname, _ = QFileDialog.getOpenFileName(self, 'Open File')  # retrieve file name

        fname = self.fname.split('/')[-1]  # used to check file type

        if self.fname.endswith('.h5') or self.fname.endswith('.all'):  # check for proper file extension
            if self.fname.endswith('.h5') and self.fname != 'demo.h5':

                n = len(fname)
                f = self.fname[:-n]
                struct_names, mag_names = mm._use_given_ff(f)  # look for form factors in directory

                # set form factors found in file directory
                self._sampleWidget.struct_ff = struct_names
                self._sampleWidget.mag_ff = mag_names

                # read in project information
                self.sample = ds.ReadSampleHDF5(self.fname)

                self.sample.energy_shift()
                fitParams = ds.ReadFitHDF5(self.fname)
                self._sampleWidget.sample = self.sample
                self._reflectivityWidget.sample = self.sample
                self._goWidget.sample = self.sample

                self.data, self.data_dict, self.sim_dict = ds.ReadDataHDF5(self.fname)  # reset data and sim info

                # loading in the background shift and scaling factor
                self._reflectivityWidget.bs = dict()
                self._reflectivityWidget.sf = dict()


                for scan_name in fitParams[4]:
                    if scan_name in list(self.data_dict.keys()):
                        self._reflectivityWidget.sf[scan_name] = str(self.data_dict[scan_name]['Scaling Factor'])
                        self._reflectivityWidget.bs[scan_name] = str(self.data_dict[scan_name]['Background Shift'])

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

                self._sampleWidget.eShift = dict()  # make sure we clear the eshift
                self._sampleWidget.ffScale = dict()
                for key in list(self.sample.eShift.keys()):
                    name = 'ff-' + key
                    self._sampleWidget.eShift[name] = self.sample.eShift[key]
                    self._sampleWidget.ffScale[name] = self.sample.ff_scale[key]

                for key in list(self.sample.mag_eShift.keys()):
                    name = 'ffm-' + key
                    self._sampleWidget.eShift[name] = self.sample.mag_eShift[key]
                    self._sampleWidget.ffScale[name] = self.sample.ffm_scale[key]

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

                # for now let's clear all the fitting parameters
                #self._reflectivityWidget.sfBsFitParams = fitParams[0]
                self._reflectivityWidget.sfBsFitParams = []
                #self._reflectivityWidget.currentVal = list(fitParams[1])
                self._reflectivityWidget.currentVal = []
                self._reflectivityWidget.rom = [True, False, False]

                # checks to make sure that saved scans are still in the dataset
                #selectedScans = [scan for scan in fitParams[4] if scan in list(self.data_dict.keys())]
                selectedScans = []
                boundaries = []
                self._reflectivityWidget.fit = selectedScans
                #self._reflectivityWidget.fit = []

                self._reflectivityWidget.selectedScans.blockSignals(True)
                self._reflectivityWidget.selectedScans.addItems(selectedScans)
                self._reflectivityWidget.selectedScans.blockSignals(False)

                #self._reflectivityWidget.bounds = fitParams[5]
                self._reflectivityWidget.bounds = []
                #self._reflectivityWidget.bounds = [[[0.01125,0.55253]],[[0.01125,0.7]]]
                #self._reflectivityWidget.weights = fitParams[6]
                self._reflectivityWidget.weights = []

                # self._reflectivityWidget.setTable()
                self._sampleWidget.parameterFit = fitParams[2]
                self._sampleWidget.currentVal = fitParams[3]

                self._goWidget.x = fitParams[7]
                self._goWidget.fun = fitParams[8]

                self._goWidget.parameters = []
                for param in self._sampleWidget.parameterFit:
                    self._goWidget.parameters.append(param)
                for param in self._reflectivityWidget.sfBsFitParams:
                    self._goWidget.parameters.append(param)

                self._sampleWidget.setTableEShift()  # reset the energy shift table
                self._noiseWidget._resetVariables(self.data_dict, selectedScans, boundaries)
            elif fname.endswith('.all'):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("ReMagX")
                messageBox.setText("Application cannot properly load in ReMagX workspace. To load data from a ReMagX file please select the 'Load ReMagX' in the File tab.")
                messageBox.exec()
            else:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Improper File Type")
                messageBox.setText("File type not supported by the application. Workspace file must be an HDF5 file type. The HDF5 architecture can be found in the user manual. ")
                messageBox.exec()
        self.activate_tab_1()
    def _loadSample(self):
        """
            Purpose: Load a new file or project workspace
        """

        self.fname, _ = QFileDialog.getOpenFileName(self, 'Open File')  # retrieve file name

        fname = self.fname.split('/')[-1]  # used to check file type


        if self.fname.endswith('.h5') or self.fname.endswith('.all'):  # check for proper file extension
            if self.fname.endswith('.h5') and self.fname != 'demo.h5':
                n= len(fname)
                f = self.fname[:-n]

                struct_names, mag_names = mm._use_given_ff(f)  # look for form factors in directory

                # set form factors found in file directory
                self._sampleWidget.struct_ff = struct_names
                self._sampleWidget.mag_ff = mag_names

                # read in project information
                self.sample = ds.ReadSampleHDF5(self.fname)

                self.sample.energy_shift()


                self._sampleWidget.sample = self.sample
                self._reflectivityWidget.sample = self.sample
                self._goWidget.sample = self.sample



                # loading in the background shift and scaling factor
                self._reflectivityWidget.bs = dict()
                self._reflectivityWidget.sf = dict()

                self._sampleWidget._setStructFromSample(self.sample)
                self._sampleWidget._setVarMagFromSample(self.sample)

                ## now it's time to load the other information
                #self._reflectivityWidget.selectedScans.clear()
                #self._reflectivityWidget.whichScan.clear()


                self._sampleWidget.eShift = dict()  # make sure we clear the eshift
                self._sampleWidget.ffScale = dict()
                for key in list(self.sample.eShift.keys()):
                    name = 'ff-' + key
                    self._sampleWidget.eShift[name] = self.sample.eShift[key]
                    self._sampleWidget.ffScale[name] = self.sample.ff_scale[key]

                for key in list(self.sample.mag_eShift.keys()):
                    name = 'ffm-' + key
                    self._sampleWidget.eShift[name] = self.sample.mag_eShift[key]
                    self._sampleWidget.ffScale[name] = self.sample.ffm_scale[key]

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


                # for now let's clear all the fitting parameters
                # self._reflectivityWidget.sfBsFitParams = fitParams[0]
                self._reflectivityWidget.sfBsFitParams = []
                # self._reflectivityWidget.currentVal = list(fitParams[1])
                self._reflectivityWidget.currentVal = []
                self._reflectivityWidget.rom = [True, False, False]

                # checks to make sure that saved scans are still in the dataset
                # selectedScans = [scan for scan in fitParams[4] if scan in list(self.data_dict.keys())]

                #self._reflectivityWidget.fit = []


                # self._reflectivityWidget.bounds = fitParams[5]
                self._reflectivityWidget.bounds = []
                self._reflectivityWidget.weights = []

                # self._reflectivityWidget.setTable()
                self._sampleWidget.parameterFit = []
                self._sampleWidget.currentVal = []

                self._goWidget.x = []
                self._goWidget.fun = 0

                self._goWidget.parameters = []
                self._sampleWidget.parameterFit = []

                self._sampleWidget.setTableEShift()  # reset the energy shift table
                #self._noiseWidget._resetVariables(self.data_dict, [], [])

            elif fname.endswith('.all'):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("ReMagX")
                messageBox.setText(
                    "Application cannot load sample information from ReMagX file.")
                messageBox.exec()
            else:
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Improper File Type")
                messageBox.setText(
                    "File type not supported by the application. Workspace file must be an HDF5 file type. The HDF5 architecture can be found in the user manual. ")
                messageBox.exec()
        self.activate_tab_1()
    def _saveFile(self):
        """
        Purpose: Save the current project space
        """

        filename = self.fname  # retrieve the current file name

        if not(filename.endswith('demo.h5')):  # makes sure the current workspace is not the demo
            self.sample = self._sampleWidget._createSample()  # retrieve sample info
            self._sampleWidget.sample = self.sample
            self._reflectivityWidget.sample = self.sample

            # save the sample information to the file
            # ds.WriteSampleHDF5(self.fname, self.sample, self.version)

            data_dict = self.data_dict

            # retrieve fitting parameters
            fitParams = [self._reflectivityWidget.sfBsFitParams, self._reflectivityWidget.currentVal,
                         self._sampleWidget.parameterFit, self._sampleWidget.currentVal,
                         self._reflectivityWidget.fit, self._reflectivityWidget.bounds,
                         self._reflectivityWidget.weights, self._goWidget.x, self._goWidget.fun]

            optParams = self._goWidget.goParameters  # retrieve data fitting algorithm parameters

            ds.saveFileHDF5(filename, self.sample, data_dict, fitParams, optParams, self.version)  # save the information
        else:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Create New File")
            messageBox.setText("User cannot save work to current file. Please save workspace with a new file name.")
            messageBox.exec()
        self.activate_tab_1()
    def _saveAsFile(self):
        """
        Purpose: Save project worspace to a specified name
        """
        # create a new file with the inputted
        filename, _ = QFileDialog.getSaveFileName()  # retrieves file name from user
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
            # fitting parameter information
            fitParams = [self._reflectivityWidget.sfBsFitParams, self._reflectivityWidget.currentVal,
                         self._sampleWidget.parameterFit, self._sampleWidget.currentVal,
                         self._reflectivityWidget.fit, self._reflectivityWidget.bounds,
                         self._reflectivityWidget.weights, self._goWidget.x, self._goWidget.fun]

            optParams = self._goWidget.goParameters  # data fitting algorithm information

            self.sample = self._sampleWidget._createSample()
            self._sampleWidget.sample = self.sample
            self._reflectivityWidget.sample = self.sample

            ds.saveAsFileHDF5(self.fname, self.sample, data_dict, sim_dict, fitParams, optParams, self.version)  # saving
        self.activate_tab_1()
    def _saveSimulation(self):
        """
        Purpose: Calculate and save the simulation to the current file
        """
        sim_dict = copy.deepcopy(self.data_dict)  # get simulation dictionary

        fname = self.fname  # retrieve filge name

        if len(sim_dict) != 0:
            # initializing the loading screen
            loadingApp = LoadingScreen(self.sample, sim_dict)
            loadingApp.show()
            loadingApp.exec_()
            sim_dict = loadingApp.sim_dict
            loadingApp.close()

            # takes into account user may exit screen
            if type(sim_dict) is not list:
                ds.saveSimulationHDF5(self.fname, sim_dict, self.version)

    def _saveSample(self):
        """
        Purpose: Save the sample information from the current project space

        """
        self.sample = self._sampleWidget._createSample()
        self._sampleWidget.sample = self.sample
        self._reflectivityWidget.sample = self.sample

        # save the sample information to the file
        ds.WriteSampleHDF5(self.fname, self.sample, self.version)
        self.activate_tab_1()

    def _importDataSet(self):
        """
        Purpose: import data from h5 filetype
        :return:
        """
        # Import the data set
        filename, _ = QFileDialog.getOpenFileName(self, 'Open File')
        fname = filename.split('/')[-1]

        # when loading files I need to be able to scan the entire
        if filename.endswith('.h5') or filename.endswith('.all'):
            if fname.endswith('.h5') and fname != 'demo.h5':
                self.data, self.data_dict, self.sim_dict = ds.LoadDataHDF5(filename)
                self._reflectivityWidget.data = self.data
                self._reflectivityWidget.data_dict = self.data_dict
                for scan in self.data:
                    self._reflectivityWidget.whichScan.addItem(scan[2])
        self.activate_tab_1()
    def _loadReMagX(self):
        """
        Purpose: import data from a ReMagX file type
        """
        filename, _ = QFileDialog.getOpenFileName(self, 'Open File')  # retrieves file

        if filename.endswith('.all'):  # checks to make sure it is a ReMagX file type
            data, data_dict = ds.Read_ReMagX(filename)  # read the file
            self.data = copy.deepcopy(data)
            self.data_dict = copy.deepcopy(data_dict)
            self.sim_dict = copy.deepcopy(data_dict)

            self._reflectivityWidget.data = self.data
            self._reflectivityWidget.data_dict = self.data_dict
            self._reflectivityWidget.bounds = []
            self._reflectivityWidget.weights = []
            self._reflectivityWidget.fit = []


            # loading in the background shift and scaling factor
            self._reflectivityWidget.bs = dict()
            self._reflectivityWidget.sf = dict()

            # make sure we are not plotting when we do not need to
            self._reflectivityWidget.whichScan.blockSignals(True)
            self._reflectivityWidget.whichScan.clear()
            self._reflectivityWidget.selectedScans.clear()
            for scan in self.data:
                self._reflectivityWidget.whichScan.addItem(scan[2])
            self._reflectivityWidget.whichScan.blockSignals(False)

            self._noiseWidget._resetVariables(self.data_dict, [''], [])
        else:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("ReMagX File")
            messageBox.setText("Please select a file with .all extension!")
            messageBox.exec()
        self.activate_tab_1()
    def _summary(self):
        """
        Purpose: save summary of project worspace as a textfile
        """
        sample = self._sampleWidget._createSample()

        filename, _ = QFileDialog.getSaveFileName()  # retrieve file name from user
        fname = filename.split('/')[-1]
        cont = True

        if fname == '':
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
                file.write("magEnergyShift = %s \n" % sample.mag_eShift)
                file.write("ffScale = %s \n" % sample.ff_scale)
                file.write("ffmScale = %s \n\n" % sample.ffm_scale)

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
        """
        Purpose: initialize the script widget
        """

        script = scriptWidget(self._sampleWidget)
        script.show()
        script.exec_()
        script.close()
        checkscript(self.sample)

    def _showFormFactor(self):
        """
        Purpose: initialize the form factor widget
        """
        self.sample = copy.deepcopy(self._sampleWidget._createSample())

        formFactor = showFormFactors(self.sample)
        formFactor.show()
        formFactor.exec_()
        formFactor.close()

    def _license(self):
        """
        Purpose: demonstrate license if ever obtained
        """
        license = licenseWidget()
        license.show()
        license.exec_()
        license.close()
        # no current implementation
        pass

    def _help(self):
        """
        Purpose: provide user document of how to use application
        """
        from PyQt5.QtGui import QDesktopServices
        from PyQt5.QtCore import QUrl

        file_path = 'Documentation_Pythonreflectivity.pdf'

        QDesktopServices.openUrl(QUrl.fromLocalFile(file_path))
        # not currently implemented


    def activate_tab_1(self):
        # sample setup
        self.sample = copy.deepcopy(self._sampleWidget._createSample())
        self._reflectivityWidget.sample = self.sample
        self._goWidget.sample = self.sample
        self._goWidget.clearTableFit()
        self._sampleWidget.step_size.setText(self._sampleWidget._step_size)
        self.sampleButton.setStyleSheet("background-color : magenta")
        self.reflButton.setStyleSheet("background-color : pink")
        self.smoothButton.setStyleSheet("background-color : pink")
        self.goButton.setStyleSheet("background-color : pink")
        self.progressButton.setStyleSheet("background-color: pink")
        self.stackedlayout.setCurrentIndex(0)

    def activate_tab_2(self):
        # reflectivity

        not_empty, equal_elements = self._sampleWidget.check_element_number()
        if not_empty and equal_elements:
            self.sample = copy.deepcopy(self._sampleWidget._createSample())
            self._reflectivityWidget.sample = copy.deepcopy(self.sample)
            self._goWidget.sample = copy.deepcopy(self.sample)
            self._goWidget.clearTableFit()
            self._reflectivityWidget.myPlotting()
            self._reflectivityWidget.stepWidget.setText(self._sampleWidget._step_size)
            self.sampleButton.setStyleSheet("background-color : pink")
            self.reflButton.setStyleSheet("background-color : magenta")
            self.smoothButton.setStyleSheet("background-color : pink")
            self.goButton.setStyleSheet("background-color : pink")
            self.progressButton.setStyleSheet("background-color: pink")
            self.stackedlayout.setCurrentIndex(1)
        else:
            if not(not_empty):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Empty Sample Definition")
                messageBox.setText("A sample definition is required to use any other workspace.")
                messageBox.exec()
            elif not(equal_elements):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Invalid Sample Definition")
                messageBox.setText("Each layer must have the same number of elements defined in each layer. A dummy variable can be used to meet these requirements (A,D,E,G,J,L,M,Q,R,T,X,Z). ")
                messageBox.exec()

    def activate_tab_3(self):
        # Best fit curve
        not_empty, equal_elements = self._sampleWidget.check_element_number()
        if not_empty and equal_elements:
            self.sample = copy.deepcopy(self._sampleWidget._createSample())
            self._reflectivityWidget.sample = self.sample
            self._goWidget.sample = self.sample

            self.sampleButton.setStyleSheet("background-color : pink")
            self.reflButton.setStyleSheet("background-color : pink")
            self.smoothButton.setStyleSheet("background-color : magenta")
            self.goButton.setStyleSheet("background-color : pink")
            self.progressButton.setStyleSheet("background-color: pink")

            #self._noiseWidget.selectedScans = self._reflectivityWidget.fit
            self._noiseWidget._resetVariables(self.data_dict, self._reflectivityWidget.fit, self._reflectivityWidget.bounds)

            self.stackedlayout.setCurrentIndex(2)
        else:
            if not(not_empty):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Empty Sample Definition")
                messageBox.setText("A sample definition is required to use any other workspace.")
                messageBox.exec()
            elif not(equal_elements):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Invalid Sample Definition")
                messageBox.setText("Each layer must have the same number of elements defined in each layer. A dummy variable can be used to meet these requirements (A,D,E,G,J,L,M,Q,R,T,X,Z). ")
                messageBox.exec()
    def activate_tab_4(self):
        # global optimization
        not_empty, equal_elements = self._sampleWidget.check_element_number()
        if not_empty and equal_elements:
            self.sample = copy.deepcopy(self._sampleWidget._createSample())
            self._reflectivityWidget.sample = self.sample
            self._goWidget.sample = self.sample

            self.sampleButton.setStyleSheet("background-color : pink")
            self.reflButton.setStyleSheet("background-color : pink")
            self.smoothButton.setStyleSheet("background-color : pink")
            self.goButton.setStyleSheet("background-color : magenta")
            self.progressButton.setStyleSheet("background-color: pink")
            self._goWidget.setTableFit()
            self._goWidget.updateScreen()

            self.stackedlayout.setCurrentIndex(3)
        else:
            if not(not_empty):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Empty Sample Definition")
                messageBox.setText("A sample definition is required to use any other workspace.")
                messageBox.exec()
            elif not(equal_elements):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Invalid Sample Definition")
                messageBox.setText("Each layer must have the same number of elements defined in each layer. A dummy variable can be used to meet these requirements (A,D,E,G,J,L,M,Q,R,T,X,Z). ")
                messageBox.exec()
    def activate_tab_5(self):
        # optimization
        not_empty, equal_elements = self._sampleWidget.check_element_number()
        if not_empty and equal_elements:
            self.sample = copy.deepcopy(self._sampleWidget._createSample())
            self._reflectivityWidget.sample = self.sample
            self._goWidget.sample = self.sample
            self._progressWidget.reset_fit_scans()  # resets the scans in the fit
            state = self._goWidget.checkBox.checkState()
            my_state = False
            if state > 0:
                my_state = True

            self._progressWidget.check_script_state(my_state)
            self.sampleButton.setStyleSheet("background-color : pink")
            self.reflButton.setStyleSheet("background-color : pink")
            self.smoothButton.setStyleSheet("background-color : pink")
            self.goButton.setStyleSheet("background-color : pink")
            self.progressButton.setStyleSheet("background-color: magenta")

            self.stackedlayout.setCurrentIndex(4)

        else:
            if not(not_empty):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Empty Sample Definition")
                messageBox.setText("A sample definition is required to use any other workspace.")
                messageBox.exec()
            elif not(equal_elements):
                messageBox = QMessageBox()
                messageBox.setWindowTitle("Invalid Sample Definition")
                messageBox.setText("Each layer must have the same number of elements defined in each layer. A dummy variable can be used to meet these requirements (A,D,E,G,J,L,M,Q,R,T,X,Z). ")
                messageBox.exec()


class progressWidget(QWidget):
    """
    Purpose: Used to update the user of the data fitting progress
    """
    def __init__(self, sWidget, rWidget, nWidget):
        super().__init__()
        # which plot
        self.sWidget = sWidget
        self.rWidget = rWidget
        self.nWidget = nWidget
        self.y_scale = 'log(x)'
        self.whichPlot = [True, False, False, False, False]
        # parameters required for calculations
        self.run_script = False
        self.sample = None
        self.scans = None
        self.data = None
        self.backS = None
        self.scaleF = None
        self.parameters = None
        self.sBounds = None
        self.sWeights = None
        self.objective = None
        self.shape_weight = None
        self.keep_going = True
        self.script_state = False

        self.objFun = dict()
        self.costFun = dict()
        self.varFun = dict()
        self.par = []  # keep track of the names

        pagelayout = QHBoxLayout()

        buttonLayout = QVBoxLayout()

        # plot objective function
        self.objButton = QPushButton('Total Cost')
        self.objButton.clicked.connect(self._setObj)
        self.objButton.setFixedWidth(200)
        self.objButton.setStyleSheet('background: blue; color: white')
        buttonLayout.addWidget(self.objButton)

        # plot total cost function
        self.costButton = QPushButton('Norm')
        self.costButton.clicked.connect(self._setCost)
        self.costButton.setFixedWidth(200)
        self.costButton.setStyleSheet('background: lightGrey')
        buttonLayout.addWidget(self.costButton)

        # plot shape parameterization
        self.varButton = QPushButton('Variation')
        self.varButton.clicked.connect(self._setVar)
        self.varButton.setFixedWidth(200)
        self.varButton.setStyleSheet('background: lightGrey')
        buttonLayout.addWidget(self.varButton)
        buttonLayout.addSpacing(20)
        # show variation in parameters
        self.parButton = QPushButton('Parameters')
        self.parButton.clicked.connect(self._setPar)
        self.parButton.setFixedWidth(200)
        self.parButton.setStyleSheet('background: lightGrey')
        buttonLayout.addWidget(self.parButton)

        # show density profile
        self.denseButton = QPushButton('Density Profile')
        self.denseButton.clicked.connect(self._setDensityProfile)
        self.denseButton.setFixedWidth(200)
        self.denseButton.setStyleSheet('background: lightGrey')
        buttonLayout.addWidget(self.denseButton)
        buttonLayout.addSpacing(20)
        my_scans = list(self.rWidget.data_dict.keys())  # all scans


        # plot current scan progress
        scanBoxLabel = QLabel('Selected Scans: ')
        self.scanBox = QComboBox()
        self.scanBox.addItems(self.rWidget.fit)
        self.scanBox.activated.connect(self.plot_scan)
        self.scanBox.setFixedWidth(200)
        buttonLayout.addWidget(scanBoxLabel)
        buttonLayout.addWidget(self.scanBox)
        buttonLayout.addSpacing(10)

        allBoxLabel = QLabel('All Scans: ')
        self.allScans = QComboBox()
        self.allScans.addItems(my_scans)
        self.allScans.activated.connect(self.plot_scans_all)
        self.allScans.setFixedWidth(200)
        buttonLayout.addWidget(allBoxLabel)
        buttonLayout.addWidget(self.allScans)

        buttonLayout.addStretch(1)

        # plot that will show current progress
        self.plotWidget = pg.PlotWidget()
        self.plotWidget.setBackground('w')
        self.plotWidget.addLegend()

        pagelayout.addLayout(buttonLayout)
        pagelayout.addWidget(self.plotWidget)

        self.setLayout(pagelayout)


    def check_script_state(self, state):
        self.script_state = state

    def plot_scan(self):
        """
        Purpose: plot the data and current iteration of the data fitting
        """

        # script checker
        script, problem, my_error = checkscript(self.sample)

        use_script=False
        if not(problem) and self.script_state:
            use_script=True

        #script, problem = self.
        self.scanBox.setStyleSheet('background: red; selection-background-color: grey')
        self.allScans.setStyleSheet('background: white; selection-background-color: red')

        # retrieve current iteration of the data fitting
        global x_vars
        x = copy.deepcopy(x_vars)
        if len(x) == 0:
            x = go.return_x()

        self.plotWidget.clear()

        step_size = float(self.sWidget._step_size)
        name = self.scanBox.currentText()
        b_idx = self.scanBox.currentIndex()

        sample, backS, scaleF = go.changeSampleParams(x[-1], self.parameters, self.sample,
                                                      self.backS, self.scaleF, script, use_script=use_script)
        background_shift = float(backS[name])
        scaling_factor = float(scaleF[name])

        if name != '':
            bound = self.rWidget.bounds[b_idx]
            lower = float(bound[0][0])
            upper = float(bound[-1][-1])
            background_shift = self.rWidget.data_dict[name]['Background Shift']
            scaling_factor = self.rWidget.data_dict[name]['Scaling Factor']
            idx = 0
            notDone = True
            while notDone and idx == len(self.data) - 1:
                temp_name = self.data[idx][2]
                if temp_name == name:
                    notDone = False
                else:
                    idx = idx + 1

            dat = self.rWidget.data_dict[name]['Data']
            pol = self.rWidget.data_dict[name]['Polarization']

            if 'Angle' not in list(self.rWidget.data_dict[name].keys()):
                qz = dat[0]
                R = dat[2]
                E = self.rWidget.data_dict[name]['Energy']
                qz, Rsim = sample.reflectivity(E, qz, s_min=step_size, bShift=background_shift,
                                               sFactor=scaling_factor)
                Theta = np.arcsin(qz / (E * 0.001013546143)) * 180 / np.pi

                Rsim = Rsim[pol]
                if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                    else:
                        self.plotWidget.plot(Theta, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                        lower = np.arcsin(lower / (E * 0.001013546143)) * 180 / np.pi
                        upper = np.arcsin(upper / (E * 0.001013546143)) * 180 / np.pi
                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.plotWidget.setLogMode(False, True)

                    self.plotWidget.setXRange(lower, upper)
                elif pol == 'AL' or pol == 'AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2),
                                             name='Simulation')
                    else:
                        self.plotWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2),
                                             name='Simulation')

                    self.plotWidget.setLogMode(False, False)
                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.plotWidget.setXRange(lower, upper)
            else:
                E = dat[3]
                R = dat[2]
                Theta = self.rWidget.data_dict[name]['Angle']
                E, Rsim = sample.energy_scan(Theta, E, s_min=step_size, bShift=background_shift,
                                             sFactor=scaling_factor)
                Rsim = Rsim[pol]
                self.plotWidget.setLogMode(False, False)
                self.plotWidget.plot(E, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                self.plotWidget.plot(E, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                self.plotWidget.setLabel('left', "Reflectivity, R")
                self.plotWidget.setLabel('bottom', "Energy, E (eV)")

    def plot_scans_all(self):
        """
        Purpose: plot the data and current iteration of the data fitting
        """
        script, problem, my_error = checkscript(self.sample)
        use_script=False

        if not(problem) and self.script_state:
            use_script=True

        self.scanBox.setStyleSheet('background: white; selection-background-color: grey')
        self.allScans.setStyleSheet('background: red; selection-background-color: red')
        # retrieve current iteration of the data fitting
        global x_vars
        x = copy.deepcopy(x_vars)
        if len(x) == 0:
            x = go.return_x()

        self.plotWidget.clear()

        step_size = float(self.sWidget._step_size)
        name = self.allScans.currentText()

        sample, backS, scaleF = go.changeSampleParams(x[-1], self.parameters, self.sample,
                                                      self.backS, self.scaleF,script,use_script=use_script)
        background_shift = 0
        scaling_factor = 1

        if name != '':

            background_shift = self.rWidget.data_dict[name]['Background Shift']
            scaling_factor = self.rWidget.data_dict[name]['Scaling Factor']
            idx = 0
            notDone = True
            while notDone and idx == len(self.data) - 1:
                temp_name = self.data[idx][2]
                if temp_name == name:
                    notDone = False
                else:
                    idx = idx + 1

            dat = self.rWidget.data_dict[name]['Data']
            pol = self.rWidget.data_dict[name]['Polarization']

            if 'Angle' not in list(self.rWidget.data_dict[name].keys()):
                qz = dat[0]
                R = dat[2]
                E = self.rWidget.data_dict[name]['Energy']
                qz, Rsim = sample.reflectivity(E, qz, s_min=step_size, bShift=background_shift,
                                               sFactor=scaling_factor)
                Theta = np.arcsin(qz / (E * 0.001013546143)) * 180 / np.pi

                Rsim = Rsim[pol]
                if pol == 'S' or pol == 'P' or pol == 'LC' or pol == 'RC':

                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(qz, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                    else:
                        self.plotWidget.plot(Theta, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(Theta, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')

                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")
                    self.plotWidget.setLogMode(False, True)

                elif pol == 'AL' or pol == 'AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    if self.rWidget.axis_state:
                        self.plotWidget.plot(qz[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(qz[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2),
                                             name='Simulation')
                    else:
                        self.plotWidget.plot(Theta[rm_idx], R[rm_idx], pen=pg.mkPen((0, 2), width=2), name='Data')
                        self.plotWidget.plot(Theta[rm_idx], Rsim[rm_idx], pen=pg.mkPen((1, 2), width=2),
                                             name='Simulation')

                    self.plotWidget.setLogMode(False, False)
                    self.plotWidget.setLabel('left', "Reflectivity, R")
                    self.plotWidget.setLabel('bottom', "Momentum Transfer, qz (Å^{-1})")

            else:
                E = dat[3]
                R = dat[2]
                Theta = self.rWidget.data_dict[name]['Angle']
                E, Rsim = sample.energy_scan(Theta, E, s_min=step_size, bShift=background_shift,
                                             sFactor=scaling_factor)
                Rsim = Rsim[pol]
                self.plotWidget.setLogMode(False, False)
                self.plotWidget.plot(E, R, pen=pg.mkPen((0, 2), width=2), name='Data')
                self.plotWidget.plot(E, Rsim, pen=pg.mkPen((1, 2), width=2), name='Simulation')
                self.plotWidget.setLabel('left', "Reflectivity, R")
                self.plotWidget.setLabel('bottom', "Energy, E (eV)")

    def reset_fit_scans(self):
        """
        Purpose: reset the scans in scanBox
        """

        my_scans = list(self.rWidget.data_dict.keys())
        self.scanBox.blockSignals(True)
        self.allScans.blockSignals(True)
        self.scanBox.clear()
        self.allScans.clear()
        self.scanBox.setFixedWidth(200)
        self.scanBox.addItems(self.rWidget.fit)
        self.allScans.setFixedWidth(200)
        self.allScans.addItems(my_scans)
        self.scanBox.blockSignals(False)
        self.allScans.blockSignals(False)
        self.boundaries = self.rWidget.bounds

    def computeScan(self, x_array, script, use_script=False):
        """
        Purpose: calculate the cost function of the current data fitting iteration
        :param x_array: current iteration parameters
        """

        smooth_dict = self.nWidget.smoothScans  # retrieve the smoothed data

        if len(x_array) != 0:
            # compute the scans for all the new x values
            n = len(self.objFun['total'])

            for x in x_array[n:]:
                sample, backS, scaleF = go.changeSampleParams(x, self.parameters, self.sample,
                                                              self.backS, self.scaleF,script,use_script=use_script)
                gamma = 0
                fun = 0

                for i, scan in enumerate(self.scans):
                    name = scan
                    Rsmooth = smooth_dict[name]['Data'][2]
                    fun_val = 0
                    xbound = self.sBounds[i]
                    weights = self.sWeights[i]

                    background_shift = float(backS[name])
                    scaling_factor = float(scaleF[name])

                    if 'Angle' not in self.data.keys():
                        myDataScan = self.data[name]
                        myData = myDataScan['Data']
                        E = myDataScan['Energy']
                        pol = myDataScan['Polarization']
                        Rdat = np.array(myData[2])

                        qz = np.array(myData[0])
                        qz, Rsim = sample.reflectivity(E, qz, bShift=background_shift, sFactor=scaling_factor)
                        Rsim = Rsim[pol]

                        j = [x for x in range(len(qz)) if qz[x] > xbound[0][0] and qz[x] < xbound[-1][-1]]
                        qz = qz[j]
                        if len(Rsim) != len(j):
                            Rsim = Rsim[j]
                        if len(Rdat) != len(j):
                            Rdat = Rdat[j]
                        if len(Rsmooth) != len(j):
                            Rsmooth = Rsmooth[j]

                        if self.y_scale == 'log(x)':
                            Rsim = np.log10(Rsim)
                            Rdat = np.log10(Rdat)
                            Rsmooth = np.log10(Rsmooth)
                        elif self.y_scale == 'ln(x)':
                            Rsim = np.log(Rsim)
                            Rdat = np.log(Rdat)
                            Rsmooth = np.log(Rsmooth)
                        elif self.y_scale == 'qz^4':
                            Rsim = np.multiply(Rsim, np.power(qz, 4))
                            Rdat = np.multiply(Rdat, np.power(qz, 4))
                            Rsmooth = np.multiply(Rsmooth, np.power(qz,4))
                        elif self.y_scale == 'x':
                            pass

                        #window = 5
                        #R = go.rolling_average(Rdat, window)
                        # total variation
                        #var_idx = [x for x in range(len(qz)) if qz[x] >= xbound[0][0] and qz[x] < xbound[-1][1]]
                        if len(Rsmooth) != 0:
                            val = go.total_variation(Rsmooth, Rsim) / len(Rsmooth)
                        else:
                            val = 0
                        self.varFun[name].append(val * self.shape_weight)
                        gamma = gamma + val

                        m = 0
                        for b in range(len(xbound)):
                            lw = xbound[b][0]
                            up = xbound[b][1]
                            w = weights[b]

                            idx = [x for x in range(len(qz)) if
                                   qz[x] >= lw and qz[x] < up]  # determines index boundaries

                            k = len(idx)
                            m = m + k
                            if len(idx) != 0:
                                if self.objective == 'Chi-Square':
                                    fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2 / abs(Rsim[idx])) * w

                                elif self.objective == 'L1-Norm':
                                    fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                                elif self.objective == 'L2-Norm':
                                    fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2) * w

                        if m != 0:
                            self.costFun[name].append(fun_val / m)
                            self.objFun[name].append(fun_val / m + val * self.shape_weight)
                            fun = fun + fun_val / m

                    else:
                        myDataScan = self.data[name]
                        myData = myDataScan['Data']
                        Theta = myDataScan['Angle']
                        Rdat = np.array(myData[2])
                        E = np.array(myData[3])
                        pol = myDataScan['Polarization']

                        E, Rsim = sample.energy_scan(Theta, E)
                        Rsim = Rsim[pol]

                        j = [x for x in range(len(E)) if E[x] > xbound[0][0] and E[x] < xbound[-1][-1]]
                        E = E[j]
                        if len(Rsim) != len(j):
                            Rsim = Rsim[j]
                        if len(Rdat) != len(j):
                            Rdat = Rdat[j]
                        if len(Rsmooth) != len(j):
                            Rsmooth = Rsmooth[j]

                        if self.y_scale == 'log(x)':
                            Rsim = np.log10(Rsim)
                            Rdat = np.log10(Rdat)
                            Rsmooth = np.log10(Rsmooth)
                        elif self.y_scale == 'ln(x)':
                            Rsim = np.log(Rsim)
                            Rdat = np.log(Rdat)
                            Rsmooth = np.log(Rsmooth)
                        elif self.y_scale == 'qz^4':
                            qz = np.sin(Theta * np.pi / 180) * (E * 0.001013546143)
                            Rsim = np.multiply(Rsim, np.power(qz, 4))
                            Rdat = np.multiply(Rdat, np.power(qz, 4))
                            Rsmooth = np.multiply(Rsmooth, np.power(qz, 4))
                        elif self.y_scale == 'x':
                            pass

                        #window = 5
                        #R = go.rolling_average(Rdat, window)

                        #var_idx = [x for x in range(len(E)) if E[x] >= xbound[0][0] and E[x] < xbound[-1][1]]
                        # total variation
                        if len(Rsmooth) != 0:
                            val = go.total_variation(Rsmooth, Rsim) / len(Rsmooth)
                        else:
                            val = 0
                        self.varFun[name].append(val * self.shape_weight)
                        gamma = gamma + val

                        m = 0
                        for b in range(len(xbound)):
                            lw = xbound[b][0]
                            up = xbound[b][1]
                            w = weights[b]

                            idx = [x for x in range(len(E)) if E[x] >= lw and E[x] < up]  # determines index boundaries

                            k = len(idx)
                            m = m + k
                            if len(idx) != 0:
                                if self.objective == 'Chi-Square':
                                    fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2 / abs(Rsim[idx])) * w
                                elif self.objective == 'L1-Norm':
                                    fun_val = fun_val + sum(np.abs(Rdat[idx] - Rsim[idx])) * w
                                elif self.objective == 'L2-Norm':
                                    fun_val = fun_val + sum((Rdat[idx] - Rsim[idx]) ** 2) * w

                        if m != 0:
                            self.costFun[name].append(fun_val / m)
                            self.objFun[name].append(fun_val / m + val * self.shape_weight)
                            fun = fun + fun_val / m

                self.costFun['total'].append(fun)
                self.varFun['total'].append(gamma * self.shape_weight)
                fun = fun + gamma * self.shape_weight
                self.objFun['total'].append(fun)

    def plotProgress(self):
        """
        Purpose: plot the progress of the scans
        """
        script, problem, my_error = checkscript(self.sWidget.sample)
        use_script=False
        if not(problem) and self.script_state:
            use_script=True

        # retrieve the parameters of the current data fitting iteration
        global x_vars
        x = copy.deepcopy(x_vars)
        if len(x) == 0:
            x = go.return_x()
        self.plotWidget.clear()

        self.computeScan(x, script,use_script=use_script)

        idx = self.whichPlot.index(True)
        n = len(x)
        iterations = np.arange(1, n + 1)

        self.plotWidget.setLogMode(False, False)
        if len(x) != 0:
            if idx == 0:  # total objective function
                m = len(list(self.objFun.keys()))
                for i, key in enumerate(list(self.objFun.keys())):
                    val = self.objFun[key]
                    self.plotWidget.plot(iterations, val, pen=pg.mkPen((i, m), width=2), name=key)
                    self.plotWidget.setLabel('left', "Function")
                    self.plotWidget.setLabel('bottom', "Iteration")

            elif idx == 1:  # cost function
                m = len(list(self.costFun.keys()))
                for i, key in enumerate(list(self.costFun.keys())):
                    val = self.costFun[key]
                    self.plotWidget.plot(iterations, val, pen=pg.mkPen((i, m), width=2), name=key)
                    self.plotWidget.setLabel('left', "Function")
                    self.plotWidget.setLabel('bottom', "Iteration")

            elif idx == 2:  # shape function
                m = len(list(self.varFun.keys()))
                for i, key in enumerate(list(self.varFun.keys())):
                    val = self.varFun[key]
                    self.plotWidget.plot(iterations, val, pen=pg.mkPen((i, m), width=2), name=key)
                    self.plotWidget.setLabel('left', "Function")
                    self.plotWidget.setLabel('bottom', "Iteration")

            elif idx == 3:  # varying parameters
                m = len(self.par)
                for i in range(len(self.par)):
                    x_values = [x_val[i] for x_val in x]
                    self.plotWidget.plot(iterations, x_values, pen=pg.mkPen((i, m), width=2), name=self.par[i])
            elif idx == 4:  # plot the density profile
                sample = copy.deepcopy(self.sample)
                sample, backS, scaleF = go.changeSampleParams(x[-1], self.parameters, sample,
                                                              self.backS, self.scaleF, script, use_script=True)

                thickness, density, density_magnetic = sample.density_profile()

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
                        self.plotWidget.plot(thickness, val[idx], pen=pg.mkPen((idx, num), width=2),
                                             name=list(density.keys())[idx])
                    else:
                        self.plotWidget.plot(thickness, val[idx], pen=pg.mkPen((idx, num), width=2),
                                             name=list(density.keys())[idx])

                for idx in range(len(mag_val)):

                    myname = 'Mag: ' + list(density_magnetic.keys())[idx]
                    self.plotWidget.plot(thickness, -mag_val[idx],
                                         pen=pg.mkPen((num - idx, num), width=2, style=Qt.DashLine), name=myname)
                self.plotWidget.setLabel('left', "Density (mol/cm^3)")
                self.plotWidget.setLabel('bottom', "Thickness (Å)")

    def getName(self, p):
        """
        Purpose: creates efficient names for parameters
        :param p: list containing parameter info
        :return: string
        """
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
                    name = 'bShift-' + p[1]
            elif p[0] == 'SCALING FACTOR':  # scaling factor
                if p[1] == 'ALL SCANS':
                    name = 'sFactor-All'
                else:
                    name = 'sFactor-' + p[1]

        return name

    def startSaving(self, sample, data_dict, scans, backS, scaleF, parameters, sBounds, sWeights,
                    objective, shape_weight):
        """
        Purpose: Used to intialize data fitting updating
        :param sample: slab class
        :param data_dict: data dictionary
        :param scans: scans to fit
        :param backS: background shift
        :param scaleF: scaling factor
        :param parameters: parameters
        :param sBounds: scan boundaries
        :param sWeights: scan weights
        :param objective: objective function to use
        :param shape_weight: total variation weight
        """
        self.sample = sample
        self.scans = scans
        self.data = data_dict
        self.backS = backS
        self.scaleF = scaleF
        self.parameters = parameters
        self.sBounds = sBounds
        self.sWeights = sWeights
        self.objective = objective
        self.shape_weight = shape_weight

        # clear x and start saving as optimization started
        global x_vars
        x_vars = []

        # make call to function that will plot the optimization progress
        self.objFun = dict()
        self.objFun['total'] = []

        self.costFun = dict()
        self.costFun['total'] = []

        self.varFun = dict()
        self.varFun['total'] = []

        # initializing objective, cost, shape lists
        for scan in scans:
            name = scan
            self.objFun[name] = []
            self.costFun[name] = []
            self.varFun[name] = []

        self.par = []
        # initialize the names still
        for p in self.parameters:
            p_name = self.getName(p)
            self.par.append(p_name)

    def _setObj(self):
        """
        Purpose: plot the total cost function
        """
        self.objButton.setStyleSheet('background: blue; color: white')
        self.costButton.setStyleSheet('background: lightGrey')
        self.varButton.setStyleSheet('background: lightGrey')
        self.parButton.setStyleSheet('background: lightGrey')
        self.denseButton.setStyleSheet('background: lightGrey')

        self.whichPlot = [True, False, False, False, False]

        self.plotProgress()
        # change plot

    def _setCost(self):
        """
        Purpose: plot the norm
        """
        self.objButton.setStyleSheet('background: lightGrey')
        self.costButton.setStyleSheet('background: blue; color: white')
        self.varButton.setStyleSheet('background: lightGrey')
        self.parButton.setStyleSheet('background: lightGrey')
        self.denseButton.setStyleSheet('background: lightGrey')

        self.whichPlot = [False, True, False, False, False]

        self.plotProgress()
        # change plot

    def _setVar(self):
        """
        Purpose: plot the total variation penalty
        """
        self.objButton.setStyleSheet('background: lightGrey')
        self.costButton.setStyleSheet('background: lightGrey')
        self.varButton.setStyleSheet('background: blue; color: white')
        self.parButton.setStyleSheet('background: lightGrey')
        self.denseButton.setStyleSheet('background: lightGrey')

        self.whichPlot = [False, False, True, False, False]

        self.plotProgress()

    def _setDensityProfile(self):
        """
        Purpose: plot the density profile of the current data fitting iteration
        """
        self.objButton.setStyleSheet('background: lightGrey')
        self.costButton.setStyleSheet('background: lightGrey')
        self.varButton.setStyleSheet('background: lightGrey')
        self.parButton.setStyleSheet('background: lightGrey')
        self.denseButton.setStyleSheet('background: blue; color: white')

        self.whichPlot = [False, False, False, False, True]

        self.plotProgress()

    def stop(self):
        """
        Purpose: Stop the data fitting update process
        """
        self.keep_going = False

    def start(self):
        """
        Purpose: start the data fitting update process
        """
        self.keep_going = True

    def update_optimization(self):
        """
        Purpose: Data fitting update process
        """
        script, problem, my_error = checkscript(self.sWidget.sample)
        use_script = False
        if not(problem):
            use_script=True

        self.keep_going = True
        idx = 0
        while self.keep_going:
            idx = idx + 1
            time.sleep(0.01)

            if 51 / idx == 1:
                global x_vars
                x = copy.deepcopy(x_vars)
                if len(x) == 0:
                    x = go.return_x()

                self.computeScan(x, script, use_script=use_script)
                idx = 0
        return

    def _setPar(self):
        """
        Purpose: plot the parameter evolution as a function of iteration
        """
        self.objButton.setStyleSheet('background: lightGrey')
        self.costButton.setStyleSheet('background: lightGrey')
        self.varButton.setStyleSheet('background: lightGrey')
        self.parButton.setStyleSheet('background: blue; color: white')
        self.denseButton.setStyleSheet('background: lightGrey')

        self.whichPlot = [False, False, False, True, False]

        self.plotProgress()



class showFormFactors(QDialog):
    """
    Purpose: plot selected form factors
    """
    def __init__(self, sample):
        super().__init__()
        pageLayout = QVBoxLayout()  # page layout
        self.setWindowTitle('Form Factors')
        self.setGeometry(180, 60, 700, 400)  # window geometry

        self.selectedff = []
        self.selectedffm = []

        self.ff = mm.ff
        self.ffm = mm.ffm

        # buttons to choose structural or magnetic form factors
        buttonLayout = QHBoxLayout()
        self.structButton = QPushButton('Structural')
        self.structButton.setStyleSheet('background: magenta')
        self.structButton.clicked.connect(self._structural)
        self.magButton = QPushButton('Magnetic')
        self.magButton.setStyleSheet('background: pink')
        self.magButton.clicked.connect(self._magnetic)
        buttonLayout.addWidget(self.structButton)
        buttonLayout.addWidget(self.magButton)

        self.stackLayout = QStackedLayout()

        # structural ---------------------------------------------------------------------------------------------------
        # Show form factor Widget

        # Adding the plotting Widget
        self.structPlot = pg.PlotWidget()
        self.structPlot.setBackground('w')
        self.structPlot.addLegend()

        self.structWidget = QWidget()
        self.structLayout = QHBoxLayout()

        self.realLayout = QHBoxLayout()
        self.realLabel = QLabel('Real')
        self.real = QCheckBox()
        self.real.stateChanged.connect(self._plot_ff)
        self.realLayout.addWidget(self.realLabel)
        self.realLayout.addWidget(self.real)

        self.imLayout = QHBoxLayout()
        self.imLabel = QLabel('Imaginary')
        self.im = QCheckBox()
        self.im.stateChanged.connect(self._plot_ff)
        self.im.setCheckState(2)
        self.imLayout.addWidget(self.imLabel)
        self.imLayout.addWidget(self.im)

        self.selectFF = QPushButton('Select ff')
        self.selectFF.clicked.connect(self._selectff)
        self.removeFF = QPushButton('Remove ff')
        self.removeFF.clicked.connect(self._removeff)

        fLayout = QVBoxLayout()
        fLabel = QLabel('Form Factor:')
        self.structElements = QComboBox()
        self.structElements.addItems(list(sample.find_sf[0].values()))
        fLayout.addWidget(fLabel)
        fLayout.addWidget(self.structElements)

        fsLayout = QVBoxLayout()
        fsLabel = QLabel('Selected Form Factor:')
        self.structElementsSelect = QComboBox()
        fsLayout.addWidget(fsLabel)
        fsLayout.addWidget(self.structElementsSelect)

        structLeftLayout = QVBoxLayout()

        structLeftLayout.addLayout(fLayout)
        structLeftLayout.addWidget(self.selectFF)
        structLeftLayout.addStretch(1)
        structLeftLayout.addLayout(self.realLayout)
        structLeftLayout.addLayout(self.imLayout)
        structLeftLayout.addStretch(1)
        structLeftLayout.addLayout(fsLayout)
        structLeftLayout.addWidget(self.removeFF)

        self.structLayout.addLayout(structLeftLayout)
        self.structLayout.addWidget(self.structPlot)

        self.structWidget.setLayout(self.structLayout)

        # magnetic -----------------------------------------------------------------------------------------------------
        # Show magnetic form factor Widget
        # Adding the plotting Widget
        self.magPlot = pg.PlotWidget()
        self.magPlot.setBackground('w')
        self.magPlot.addLegend()

        magLeftLayout = QVBoxLayout()

        self.realLayoutMag = QHBoxLayout()
        self.realLabelMag = QLabel('Real')
        self.realMag = QCheckBox()
        self.realMag.stateChanged.connect(self._plot_ffm)
        self.realLayoutMag.addWidget(self.realLabelMag)
        self.realLayoutMag.addWidget(self.realMag)

        self.imLayoutMag = QHBoxLayout()
        self.imLabelMag = QLabel('Imaginary')
        self.imMag = QCheckBox()
        self.imMag.stateChanged.connect(self._plot_ffm)
        self.imMag.setCheckState(2)
        self.imLayoutMag.addWidget(self.imLabelMag)
        self.imLayoutMag.addWidget(self.imMag)

        self.selectFFM = QPushButton('Select ffm')
        self.selectFFM.clicked.connect(self._selectffm)
        self.removeFFM = QPushButton('Remove ffm')
        self.removeFFM.clicked.connect(self._removeffm)

        self.magWidget = QWidget()
        self.magLayout = QHBoxLayout()
        fmLayout = QVBoxLayout()
        fmLabel = QLabel('Form Factors:')
        self.magElements = QComboBox()
        self.magElements.addItems(list(sample.find_sf[1].values()))
        fmLayout.addWidget(fmLabel)
        fmLayout.addWidget(self.magElements)

        fsmLayout = QVBoxLayout()
        fsmLabel = QLabel('Selected Form Factors: ')
        self.magElementsSelect = QComboBox()
        fsmLayout.addWidget(fsmLabel)
        fsmLayout.addWidget(self.magElementsSelect)

        magLeftLayout.addLayout(fmLayout)
        magLeftLayout.addWidget(self.selectFFM)
        magLeftLayout.addStretch(1)
        magLeftLayout.addLayout(self.realLayoutMag)
        magLeftLayout.addLayout(self.imLayoutMag)
        magLeftLayout.addStretch(1)
        magLeftLayout.addLayout(fsmLayout)
        magLeftLayout.addWidget(self.removeFFM)

        self.magLayout.addLayout(magLeftLayout)
        self.magLayout.addWidget(self.magPlot)

        self.magWidget.setLayout(self.magLayout)

        self.stackLayout.addWidget(self.structWidget)
        self.stackLayout.addWidget(self.magWidget)

        pageLayout.addLayout(buttonLayout)
        pageLayout.addLayout(self.stackLayout)

        self._structural()

        self.setLayout(pageLayout)

    def _plot_ff(self):
        """
        Purpose: plot the structural form factors
        """
        self.structPlot.clear()  # clear the plot
        n = 0  # number of data to plot
        m = len(self.selectedff)  # number of form factors
        if self.real.checkState() != 0 and self.im.checkState() != 0:
            n = m * 2  # plot imaginary and real component
        if self.real.checkState() != 0 or self.im.checkState() != 0:
            n = m  # plot just one of real or imaginary component

        # plot the form factors
        for idx, key in enumerate(self.selectedff):
            re = self.ff[key][:, 1]
            E = self.ff[key][:, 0]
            im = self.ff[key][:, 2]

            if self.real.checkState() != 0 and self.im.checkState() != 0:

                my_name = 'r: ' + key
                self.structPlot.plot(E, re, pen=pg.mkPen((idx*2, n*2+1), width=2), name=my_name)
                my_name = 'i: ' + key
                self.structPlot.plot(E, im, pen=pg.mkPen((idx*2+1, n*2+1), width=2), name=my_name)
            elif self.real.checkState() != 0:
                my_name = 'r: ' + key
                self.structPlot.plot(E, re, pen=pg.mkPen((idx, n+1), width=2), name=my_name)
            elif self.im.checkState() != 0:
                my_name = 'i: ' + key
                self.structPlot.plot(E, im, pen=pg.mkPen((idx, n+1), width=2), name=my_name)

        self.structPlot.setLabel('left', "Form Factor")
        self.structPlot.setLabel('bottom', "Energy, (eV))")

    def _plot_ffm(self):
        """
        Purpose: plot the magnetic form factors
        """
        self.magPlot.clear()
        n = 0  # number of data to plot
        m = len(self.selectedffm)

        if self.realMag.checkState() != 0 and self.imMag.checkState() != 0:
            n = m * 2  # plot real and imaginary component
        if self.realMag.checkState() != 0 or self.imMag.checkState() != 0:
            n = m  # plot either real or imaginary component

        # plot the form factors
        for idx, key in enumerate(self.selectedffm):
            re = self.ffm[key][:, 1]
            E = self.ffm[key][:, 0]
            im = self.ffm[key][:, 2]

            if self.real.checkState() != 0 and self.im.checkState() != 0:

                my_name = 'r: ' + key
                self.magPlot.plot(E, re, pen=pg.mkPen((idx*2, n*2+1), width=2), name=my_name)
                my_name = 'i: ' + key
                self.magPlot.plot(E, im, pen=pg.mkPen((idx*2+1, n*2+1), width=2), name=my_name)
            elif self.real.checkState() != 0:
                my_name = 'r: ' + key
                self.magPlot.plot(E, re, pen=pg.mkPen((idx, n), width=2), name=my_name)

            elif self.im.checkState() != 0:
                my_name = 'i: ' + key
                self.magPlot.plot(E, im, pen=pg.mkPen((idx+1, n), width=2), name=my_name)

        self.magPlot.setLabel('left', "Form Factor")
        self.magPlot.setLabel('bottom', "Energy, (eV))")

    def _selectff(self):
        """
        Purpose: user has selected a form factor to plot
        """
        item = self.structElements.currentText()
        if item not in self.selectedff:
            self.structElementsSelect.addItem(item)
            self.selectedff.append(item)
        self._plot_ff()

    def _selectffm(self):
        """
        Purpose: user has selected a magnetic form factor to plot
        """
        item = self.magElements.currentText()
        if item not in self.selectedffm:
            self.magElementsSelect.addItem(item)
            self.selectedffm.append(item)
        self._plot_ffm()

    def _removeff(self):
        """
        Purpose: remove structural form fator from being plotted
        """
        item = self.structElementsSelect.currentText()
        idx = self.structElementsSelect.currentIndex()
        if item != '':
            self.structElementsSelect.removeItem(idx)
            self.selectedff.pop(idx)
            self._plot_ff()

    def _removeffm(self):
        """
        Purpose: remove magnetic form factor from being plotted
        """
        item = self.magElementsSelect.currentText()
        idx = self.magElementsSelect.currentIndex()
        if item != '':
            self.magElementsSelect.removeItem(idx)
            self.selectedff.pop(idx)
            self._plot_ffm()

    def _structural(self):
        """
        Purpose: activate structural form factor workspace
        """
        self.magButton.setStyleSheet('background: pink')
        self.structButton.setStyleSheet('background: magenta')
        self.stackLayout.setCurrentIndex(0)

    def _magnetic(self):
        """
        Purpose: activate magnetic form factor workspace
        """
        self.magButton.setStyleSheet('background: magenta')
        self.structButton.setStyleSheet('background: pink')
        self.stackLayout.setCurrentIndex(1)


class scriptWidget(QDialog):
    """
    Purpose: Initialize script window to allow user to perform special operations
    """
    def __init__(self, sWidget):
        super().__init__()

        cwd = os.getcwd()

        self.fname = cwd + '/default_script.txt'  # obtain current script
        self.setWindowTitle('Script Window')
        self.setGeometry(180, 60, 700, 400)
        self.sWidget = sWidget

        # open and save the script buttons
        pagelayout = QHBoxLayout()
        openButton = QPushButton('Open Script')
        openButton.setFixedWidth(100)
        openButton.clicked.connect(self.open_new_file)
        saveButton = QPushButton('Save Script')
        saveButton.setFixedWidth(100)
        saveButton.clicked.connect(self.save_file)
        checkButton = QPushButton('Check Saved Script')
        checkButton.setFixedWidth(100)
        checkButton.clicked.connect(self.check_script)
        runButton = QPushButton('Run Saved Script')
        runButton.setFixedWidth(100)
        runButton.clicked.connect(self.run_script)

        buttonLayout = QVBoxLayout()
        buttonLayout.addStretch(1)
        buttonLayout.addWidget(openButton)
        buttonLayout.addWidget(saveButton)
        buttonLayout.addWidget(checkButton)
        buttonLayout.addWidget(runButton)
        buttonLayout.addStretch(1)

        # creating the script workspace
        vbox = QVBoxLayout()
        text = 'default_script.txt'
        self.title = QLabel(text)
        self.title.setWordWrap(True)
        self.title.setAlignment(Qt.AlignCenter)
        vbox.addWidget(self.title)

        self.scrollable_text_area = QTextEdit()  # allowing scrollable capabilities

        with open(self.fname) as f:
            file_contents = f.read()
            self.scrollable_text_area.setText(file_contents)

        f.close()

        vbox.addWidget(self.scrollable_text_area)

        pagelayout.addLayout(vbox)
        pagelayout.addLayout(buttonLayout)

        self.setLayout(pagelayout)

    def open_new_file(self):
        """
        Purpose: Open a new script file
        """
        fname, filter_type = QFileDialog.getOpenFileName(self, "Open new file", "", "All files (*)")
        if self.fname.endswith('.txt'):
            with open(fname, "r") as f:
                file_contents = f.read()
                self.title.setText(fname)
                self.scrollable_text_area.setText(file_contents)
        else:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Invalid file")
            messageBox.setText("Selected filename or path is not valid. Please select a valid file.")
            messageBox.exec()

        self.fname = fname  # change the script file that will be altered

    def save_file(self):
        """
        Purpose: save current script
        """
        text = self.scrollable_text_area.toPlainText()
        with open(self.fname, 'w') as f:
            f.write(text)
        f.close()

    def run_script(self):
        """
        Purpose: Checks the script and runs it if all checks passed
        """
        sample = self.sWidget._createSample()
        my_script, problem, my_error = checkscript(sample)

        if problem:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Script unable to execute")
            messageBox.setText("Error: " + my_error['Type'] + '\n' + 'Line: ' + str(my_error['Line']))
            messageBox.exec()
        else:
            sample = go.readScript(sample, my_script)
            self.sWidget.sample = copy.deepcopy(sample)
            self.sWidget._setStructFromSample(sample)
            self.sWidget._setVarMagFromSample(sample)
            self.sWidget.setTable()
            self.sWidget.setTableMag()
            self.sWidget.setTableVar()
            self.sWidget.eShiftFromSample(sample)
            self.sWidget.setTableEShift()



    def check_script(self):
        """
        Purpose: Checks the script
        """

        my_script, problem, my_error = checkscript(self.sWidget._createSample())
        if problem:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("Script unable to execute")
            messageBox.setText("Error: " + my_error['Type'] + '\n' + 'Line: ' + str(my_error['Line']))
            messageBox.exec()


class LoadingScreen(QDialog):
    """
    Purpose: Provides progress of saving simulation
    """
    def __init__(self, sample, sim_dict):
        super().__init__()
        self.setWindowTitle('Saving Simulation')
        self.sample = sample
        self.sim_dict = []
        self.temp_sim = sim_dict
        self.n = len(list(self.temp_sim.keys()))
        layout = QHBoxLayout()
        button = QPushButton('Start Save')
        button.clicked.connect(self.run)
        self.progress = QProgressBar(self)
        self.progress.setRange(0, self.n - 1)
        self.progress.setVisible(True)
        layout.addWidget(button)
        layout.addWidget(self.progress)
        self.setLayout(layout)

    def run(self):
        """
        Purpose: calculate the simulations from current sample model
        """

        my_keys = list(self.temp_sim.keys())
        n = self.n
        if n != 0:
            for idx in range(len(my_keys)):

                if idx % 2 == 0:
                    self.progress.setValue(idx)

                key = my_keys[idx]
                pol = self.temp_sim[key]['Polarization']

                if 'Angle' in self.temp_sim[key].keys():  # energy scan
                    E = self.temp_sim[key]['Data'][3]  # get energy axis
                    Theta = self.temp_sim[key]['Angle']  # get angle
                    E, R = self.sample.energy_scan(Theta, E)  # calculate energy scan
                    R = R[pol]  # polarization
                    self.temp_sim[key]['Data'][2] = list(R)
                else:  # reflectivity scan
                    qz = self.temp_sim[key]['Data'][0]
                    energy = self.temp_sim[key]['Energy']
                    qz, R = self.sample.reflectivity(energy, qz)
                    R = R[pol]
                    self.temp_sim[key]['Data'][2] = list(R)

            self.sim_dict = copy.deepcopy(self.temp_sim)

        self.accept()

def checkbracket(myStr):
    # checks to make sure that the brackets are maintained
    open_list = ["[", "{", "("]
    close_list = ["]", "}", ")"]

    stack = []
    for i in myStr:
        if i in open_list:
            stack.append(i)
        elif i in close_list:
            pos = close_list.index(i)
            if ((len(stack) > 0) and
                    (open_list[pos] == stack[len(stack) - 1])):
                stack.pop()
            else:
                return True

    if len(stack) == 0:
        return False
    else:
        return True

def checkscript(sample):
    script = os.getcwd() + '/default_script.txt'
    my_script = list()
    my_error = dict()
    with open(script, "r") as f:
        my_line = 1
        for line in f.readlines():
            test = line.strip(' ')
            if test != "\n" and not(test.startswith('#')):

                #my_error.append({'Executable': line.strip("\n").split('='), 'Line': my_line})
                my_script.append(line.strip("\n").split('='))
                n = len(my_script)

            my_line = my_line + 1
    f.close()

    my_function_1 = ['setroughness',  'setdensity',  'setthickness', 'setcombinedthickness', 'setratio', 'seteshift', 'setmageshift', 'setmagdensity', 'setvariationconstant']

    my_function_2 = ['getroughness', 'getdensity',  'getthickness', 'gettotalthickness','geteshift', 'getmageshift', 'getmagdensity']


    problem = False

    # checks to make sure that all functions agree with what we expect
    for line in my_script:


        if len(line) == 1:
            function = line[0].split('(')[0].strip(' ')

            if not(problem):
                problem = checkbracket(line[0])
                my_error['Type'] = 'Unbalanced number of brackets'
                my_error['Line'] = my_line

            if function.lower() not in my_function_1:
                problem = True
                my_error['Type'] = 'function defined improperly'
                my_error['Line'] = my_line



            if not(problem):
                params = line[0].strip(' ').strip(function)

                params = params.strip('(')
                params = params.strip(')')
                params = params.strip(' ')
                params = params.split(',')

                if function.lower() == 'setcombinedthickness':
                    if len(params) != 4:  # not expected number of arguments
                        problem = True
                        my_error['Type'] = 'Unexpected number of arguments'
                        my_error['Line'] = my_line

                    if not (params[0].isdigit()) and not (problem):  # first two arguments are not digits
                        problem = True
                        my_error['Type'] = 'variable type'
                        my_error['Line'] = my_line

                    else:
                        if '.' in params[0] or '.' in params[1]:  # first two arguments cannot be floats
                            problem = True
                            my_error['Type'] = 'variable type'
                            my_error['Line'] = my_line



                        if not(problem):
                            start = int(params[0])
                            end = int(params[1])
                            if start > end:
                                problem = True
                                my_error['Type'] = 'layer improperly defined'
                                my_error['Line'] = my_line

                            if start < 0:
                                problem = True
                                my_error['Type'] = 'start layer must be greater than 0'
                                my_error['Line'] = my_line

                            m = len(sample.structure)


                            if m - 1 < end:
                                problem = True
                                my_error['Type'] = 'end layer larger than number of layers in sample'
                                my_error['Line'] = my_line

                            if not (problem):

                                key = params[2].strip(' ')

                                if key.lower() != 'all':
                                    for i in range(start, end + 1, 1):
                                        if key not in list(sample.structure[i].keys()):
                                            problem = True
                                            my_error['Type'] = key + ' not in layer ' + str(i)
                                            my_error['Line'] = my_line
                elif function.lower() in ['seteshift', 'setmageshift']:
                    ffName = params[0]
                    if function.lower() == 'seteshift':
                        if ffName not in list(sample.eShift.keys()):
                            problem = True
                            my_error['Type'] = ffName + ' form factor not found'
                            my_error['Line'] = my_line
                    else:
                        if ffName not in list(sample.mag_eShift.keys()):
                            problem = True
                            my_error['Type'] = ffName + ' form factor not found'
                            my_error['Line'] = my_line

                elif function.lower() == 'setmagdensity':
                    m = len(sample.structure)
                    layer = params[0]

                    if not (layer.isdigit()):
                        problem = True
                        my_error['Type'] = 'layer not entered as an integer'
                        my_error['Line'] = my_line
                    else:
                        layer = int(layer)

                    symbol = params[1].strip(' ')
                    var = params[2].strip(' ')
                    if not problem:
                        if m - 1 < layer or layer < 0:
                            problem = True
                            my_error['Type'] = 'layer improperly defined'
                            my_error['Line'] = my_line
                        if not problem:
                            if symbol not in list(sample.structure[layer].keys()):
                                problem = True
                                my_error['Type'] = str(symbol) + ' not found in layer'
                                my_error['Line'] = my_line

                            if symbol != var:
                                if not (problem) and var not in sample.structure[layer][symbol].polymorph:
                                    problem = True
                                    my_error['Type'] = str(var) + ' not found in layer'
                                    my_error['Line'] = my_line
                elif function.lower() == 'setvariationconstant':
                    if len(params) != 4:
                        problem = True
                        my_error['Type'] = 'unexpected number of parameters'
                        my_error['Line'] = my_line

                    m = len(sample.structure)

                    if int(params[0]) > m - 1 and not (problem):
                        problem = True
                        my_error['Type'] = 'layer larger than number of defined layers'
                        my_error['Line'] = my_line

                    if params[1].isdigit() or params[2].isdigit():
                        problem = True
                        my_error['Type'] = 'values must be entered as integers'
                        my_error['Line'] = my_line

                    if not (problem):
                        params[1] = params[1].strip(' ')
                        params[2] = params[2].strip(' ')
                        if params[1] not in list(sample.structure[int(params[0])].keys()):
                            problem = True
                            my_error['Type'] = str(params[1]) + ' not in layer'
                            my_error['Line'] = my_line


                        if len(sample.structure[int(params[0])][params[1]].polymorph) != 3 and not (problem):
                            problem = True
                            my_error['Type'] = 'only three element variation can exist for this function'
                            my_error['Line'] = my_line

                        if not (problem) and params[2] not in sample.structure[int(params[0])][params[1]].polymorph:
                            problem = True
                            my_error['Type'] = str(params[2]) + ' not found in layer'
                            my_error['Line'] = my_line



                elif function.lower() == 'setratio':
                    if len(params) != 5:
                        problem = True
                        my_error['Type'] = 'Unexpected number of parameters'
                        my_error['Line'] = my_line


                    if not (params[0].isdigit()) and not(problem):
                        problem = True
                        my_error['Type'] = 'variable type must be an integer'
                        my_error['Line'] = my_line

                    m = len(sample.structure)

                    if int(params[0]) > m-1 and not(problem):
                        problem = True
                        my_error['Type'] = 'layer is greater than total number of layers'
                        my_error['Line'] = my_line

                    if params[1].isdigit() or params[2].isdigit() or params[3].isdigit():
                        problem = True
                        my_error['Type'] = 'parameter types must be integers'
                        my_error['Line'] = my_line

                    if not(problem):

                        if params[1] not in list(sample.structure[int(params[0])].keys()):
                            problem = True
                            my_error['Type'] = str(params[1]) + ' not in layer'
                            my_error['Line'] = my_line

                        if len(sample.structure[int(params[0])][params[1]].polymorph) != 3 and not(problem):
                            problem=True
                            my_error['Type'] = 'this function only works with three element variations'
                            my_error['Line'] = my_line


                        if not(problem) and params[2] not in sample.structure[int(params[0])][params[1]].polymorph:
                            problem = True
                            my_error['Type'] = str(params[2]) + ' not found in layer'
                            my_error['Line'] = my_line

                        if not(problem) and params[3] not in sample.structure[int(params[0])][params[1]].polymorph:
                            problem = True
                            my_error['Type'] = str(params[3]) + ' not found in layer'
                            my_error['Line'] = my_line


                else:
                    if len(params) != 3:
                        problem = True
                        my_error['Type'] = 'unexpected number of parameters'
                        my_error['Line'] = my_line



                    if not(problem) and params[0].isdigit():
                        if '.' in params[0]:
                            problem = True
                            my_error['Type'] = 'improper variable type'
                            my_error['Line'] = my_line
                        layer = int(params[0])

                        if layer > len(sample.structure) - 1:
                            problem = True
                            my_error['Type'] = 'layer greater than total number of layers'
                            my_error['Line'] = my_line

                        if not (problem):
                            key = params[1].strip(' ')
                            if key not in list(sample.structure[layer].keys()) and key != 'all':
                                problem = True
                                my_error['Type'] = key + ' not found'
                                my_error['Line'] = my_line

                    else:
                        problem = True
                        my_error['Type'] = 'improper variable type'
                        my_error['Line'] = my_line


        elif len(line) == 2:
            if not(problem):
                problem = checkbracket(line[1])
                my_error['Type'] = 'unbalanced brackets'
                my_error['Line'] = my_line

            function = line[1].split('(')[0].strip(' ')
            if function.lower() not in my_function_2:
                problem = True
                my_error['Type'] = 'function does not exist or spelled incorrectly'
                my_error['Line'] = my_line

            if not(problem):
                params = line[1].strip(' ').strip(function)

                params = params.strip('(')
                params = params.strip(')')
                params = params.strip(' ')
                params = params.split(',')

                if function.lower() == 'gettotalthickness':
                    if len(params) != 3:  # not expected number of arguments
                        problem = True
                        my_error['Type'] = 'unexpected number of parameters'
                        my_error['Line'] = my_line
                    if not(params[0].isdigit()) or not(params[1].isdigit()) and not(problem):  # first two arguments are not digits
                        problem = True
                        my_error['Type'] = 'first two arguments must be digits'
                        my_error['Line'] = my_line
                    else:
                        if '.' in params[0] or '.' in params[1]:  # first two arguments cannot be floats
                            problem=True
                            my_error['Type'] = 'first two parameters cannot be floats'
                            my_error['Line'] = my_line


                        if not(problem):
                            start = int(params[0])
                            end = int(params[1])
                            if start > end:
                                problem = True
                                my_error['Type'] = 'start layer is larger than end layer'
                                my_error['Line'] = my_line
                            if start < 0:
                                problem = True
                                my_error['Type'] = 'start layer must be greater than 0'
                                my_error['Line'] = my_line

                            m = len(sample.structure)

                            if m-1 < end:
                                problem = True
                                my_error['Type'] = 'end layer is too small'
                                my_error['Line'] = my_line

                            if not(problem):
                                key = params[2].strip(' ')
                                if key.lower() != 'all':
                                    for i in range(start,end+1,1):
                                        if key not in list(sample.structure[i].keys()):
                                            problem = True
                                            my_error['Type'] = str(key) + ' not found in layer'
                                            my_error['Line'] = my_line

                elif function.lower() in ['geteshift', 'getmageshift']:
                    ffName = params[0]
                    if function.lower() == 'geteshift':
                        if ffName not in list(sample.eShift.keys()):
                            problem = True
                            my_error['Type'] = ffName + ' form factor not found'
                            my_error['Line'] = my_line
                    else:
                        if ffName not in list(sample.mag_eShift.keys()):
                            problem = True
                            my_error['Type'] = ffName + ' form factor not found'
                            my_error['Line'] = my_line

                elif function.lower() == 'getmagdensity':
                    m = len(sample.structure)
                    layer = params[0]

                    if not(layer.isdigit()):
                        problem = True
                        my_error['Type'] = 'layer must be an integer'
                        my_error['Line'] = my_line
                    else:
                        layer = int(layer)

                    symbol = params[1].strip(' ')
                    var = params[2].strip(' ')
                    if not problem:
                        if m-1 < layer or layer < 0:
                            problem = True
                            my_error['Type'] = 'start and end layer improperly defined'
                            my_error['Line'] = my_line

                        if not problem:
                            if symbol not in list(sample.structure[layer].keys()):
                                problem = True
                                my_error['Type'] = str(symbol) + 'not found in layer'
                                my_error['Line'] = my_line

                            if symbol != var:
                                if not(problem) and var not in sample.structure[layer][symbol].polymorph:
                                    problem = True
                                    my_error['Type'] = str(var) + 'not found in layer'
                                    my_error['Line'] = my_line
                else:
                    if len(params) != 2:
                        problem = True
                        my_error['Type'] = 'unexpected number of parameters'
                        my_error['Line'] = my_line
                    if not(problem) and params[0].isdigit():
                        if '.' in params[0]:
                            problem = True
                            my_error['Type'] = 'layer must be an integer'
                            my_error['Line'] = my_line
                        layer = int(params[0])

                        if layer > len(sample.structure)-1:
                            problem = True
                            my_error['Type'] = 'layer larger than total number of layers'
                            my_error['Line'] = my_line

                        if not(problem):
                            key = params[1].strip(' ')
                            if key.lower() != 'all':
                                if key not in list(sample.structure[layer].keys()):
                                        problem = True
                                        my_error['Type'] = str(key) + ' not found in layer'
                                        my_error['Line'] = my_line
                    else:
                        problem = True
                        my_error['Type'] = 'layer must be integer type'
                        my_error['Line'] = my_line
        else:
            problem = True
            my_error['Type'] = 'get function improperly defined'
            my_error['Line'] = my_line


    return my_script, problem, my_error


def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

class licenseWidget(QDialog):
    """
    Purpose: Initialize script window to allow user to perform special operations
    """
    def __init__(self):
        super().__init__()
        self.setWindowTitle('Script Window')
        self.setGeometry(180, 60, 700, 400)
        pagelayout = QVBoxLayout()


        # creating the script workspace
        vbox = QVBoxLayout()
        text = 'License'
        self.title = QLabel(text)
        self.title.setWordWrap(True)
        self.title.setAlignment(Qt.AlignCenter)
        vbox.addWidget(self.title)

        self.scrollable_text_area = QTextEdit()  # allowing scrollable capabilities
        self.scrollable_text_area.setReadOnly(True)

        with open('license.txt', 'r',encoding='ISO-8859-1') as f:
            file_contents = f.read()
            self.scrollable_text_area.setText(file_contents)

        vbox.addWidget(self.scrollable_text_area)

        pagelayout.addLayout(vbox)

        self.setLayout(pagelayout)

if __name__ == '__main__':
    fname = 'Pim10uc.h5'

    app = QApplication(sys.argv)
    demo = ReflectometryApp()

    demo.show()

    sys.exit(app.exec_())
