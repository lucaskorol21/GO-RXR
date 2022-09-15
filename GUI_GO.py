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
                    print(ele)
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



class reflectivityWidget(QWidget):
    def __init__(self):
        super().__init__()
        pagelayout = QHBoxLayout()
        self.layerBox = QComboBox()
        hello = QLabel('Hello')
        pagelayout.addWidget(self.layerBox)
        pagelayout.addWidget(hello)
        self.setLayout(pagelayout)

class sampleWidget(QWidget):
    def __init__(self, sample):
        super().__init__()
        self.sample = sample  # variable used to define sample info
        self.structTableInfo = []  # used to keep track of the table info instead of constantly switching

        self.previousLayer = 0
        self.firstStruct = True
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
        elementVar = QLabel('hello')
        # Magnetic Stuff
        selectlayout = QVBoxLayout()
        magVar = QLabel('bye')

        self.sampleInfoLayout.addWidget(self.structTable)  # index 1
        self.sampleInfoLayout.addWidget(elementVar)  # index 2
        self.sampleInfoLayout.addWidget(magVar)  # index 3

        # buttons for choosing which parameters to choose
        structButton = QPushButton('Structure')
        structButton.clicked.connect(self._structural)
        selectlayout.addWidget(structButton)
        polyButton = QPushButton('Element Variation')
        polyButton.clicked.connect(self._elementVariation)
        selectlayout.addWidget(polyButton)
        magButton = QPushButton('Magnetic')
        magButton.clicked.connect(self._magnetic)
        selectlayout.addWidget(magButton)
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




        self.setLayout(mylayout)


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

        for row in range(3):
            for col in range(6):
                item = QTableWidgetItem(str(tableInfo[row][col]))
                self.structTable.setItem(row,col, item)

    def changetable(self):
        idx = self.layerBox.currentIndex()

        for row in range(len(self.structTableInfo[self.previousLayer])):
            for col in range(len(self.structTableInfo[self.previousLayer][0])-1):
                item = self.structTable.item(row,col).text()
                self.structTableInfo[self.previousLayer][row][col] = str(item)
        self.previousLayer = idx

        self.setTable(idx)

    def _addLayer(self):

        addLayerApp = compoundInput()
        addLayerApp.show()
        addLayerApp.exec_()
        userinput = addLayerApp.val
        addLayerApp.close()
        print(userinput)
        if len(userinput) != 0:
            num = self.layerBox.count()
            idx = self.layerBox.currentIndex()
            if num == 0:
                self.layerBox.addItem('Substrate')
                self.structTableInfo.insert(0, userinput)

            else:
                self.layerBox.addItem('Layer ' + str(num))
                self.structTableInfo.insert(idx+1, userinput)



    def _removeLayer(self):
        num = self.layerBox.count()

        if num != 0:
            self.layerBox.removeItem(num-1)

        idx = self.layerBox.currentIndex()

        self.structTableInfo.pop(idx)
        self.setTable(idx)

    def _copyLayer(self):
        num = self.layerBox.count()
        if num == 0:
            self.layerBox.addItem('Substrate')
        else:
            self.layerBox.addItem('Layer ' + str(num))

        idx = self.layerBox.currentIndex()

        newLayer = self.structTableInfo[idx]
        self.structTableInfo.insert(idx+1,newLayer)

    def _structural(self):
        print('structural')
        self.sampleInfoLayout.setCurrentIndex(0)

    def _elementVariation(self):
        print('element variation')
        self.sampleInfoLayout.setCurrentIndex(1)

    def _magnetic(self):
        print('magnetic')
        self.sampleInfoLayout.setCurrentIndex(2)

    def _densityprofile(self):

        for row in range(len(self.structTableInfo[self.previousLayer])):
            for col in range(len(self.structTableInfo[self.previousLayer][0])-1):
                item = self.structTable.item(row, col).text()
                self.structTableInfo[self.previousLayer][row][col] = str(item)

        self._createsample()

        thickness, density, density_magnetic = self.sample.density_profile()
        self.densityWidget.clear()

        self._plotDensityProfile(thickness, density, density_magnetic)


    def _plotDensityProfile(self,thickness, density, density_magnetic):


        num = len(density)

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


        """
        for idx in range(len(mag_val)):
            plt.plot(thickness, -mag_val[idx], '--')

        center = np.zeros(len(thickness))
        plt.plot(thickness, center, 'k-.', linewidth=2)
        my_legend = list(density.keys())

        for key in list(density_magnetic.keys()):
            my_legend.append('Mag: ' + key)

        plt.legend(my_legend)
        # plt.legend(my_legend, loc='center left', bbox_to_anchor=(1.02, 0.5))
        plt.xlabel('Thickness (Angstrom)')
        plt.ylabel('Density (mol/cm^3)')

        """
    def _createsample(self):

        m = len(self.structTableInfo)
        self.sample = ms.slab(m)
        for idx in range(m):
            formula = ''
            thickness = []
            density = []
            roughness = []
            linked_roughness = []

            layer = self.structTableInfo[idx]

            for ele in range(len(layer)):
                element = layer[ele]
                formula = formula + element[0]
                if element[6] != '1':
                    formula = formula + element[6]

                thickness.append(float(element[1]))
                density.append(float(element[2]))
                roughness.append(float(element[3]))
                if element[4].isdigit():
                    linked_roughness.append(float(element[4]))
                else:
                    linked_roughness.append(False)


            self.sample.addlayer(idx,formula,thickness,density=density, roughness=roughness, linked_roughness=linked_roughness)

class ReflectometryApp(QMainWindow):
    def __init__(self, sample):
        super().__init__()

        # set the title
        self.setWindowTitle('Reflectometry of Quantum Materials')

        # set the geometry of the window
        self.setGeometry(200,80,1000,600)

        pagelayout = QVBoxLayout()
        buttonlayout = QHBoxLayout()
        self.stackedlayout = QStackedLayout()

        pagelayout.addLayout(buttonlayout)
        pagelayout.addLayout(self.stackedlayout)

        label1 = QLabel('Label 1')
        label2 = QLabel('Label 2')
        label3 = QLabel('Label 3')

        _sampleWidget = sampleWidget(sample)  # initialize the sample widget
        _reflectivityWidget = reflectivityWidget()

        sampleButton = QPushButton('Sample')
        sampleButton.setStyleSheet("background-color : pink")
        sampleButton.clicked.connect(self.activate_tab_1)
        buttonlayout.addWidget(sampleButton)
        self.stackedlayout.addWidget(_sampleWidget)

        reflButton = QPushButton('Reflectivity')
        reflButton.setStyleSheet("background-color : pink")
        reflButton.clicked.connect(self.activate_tab_2)
        buttonlayout.addWidget(reflButton)
        self.stackedlayout.addWidget(label2)

        goButton = QPushButton('Global Optimization')
        goButton.setStyleSheet("background-color : pink")
        goButton.clicked.connect(self.activate_tab_3)
        buttonlayout.addWidget(goButton)
        self.stackedlayout.addWidget(label3)

        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)

    def activate_tab_1(self):
        self.stackedlayout.setCurrentIndex(0)
    def activate_tab_2(self):
        self.stackedlayout.setCurrentIndex(1)
    def activate_tab_3(self):
        self.stackedlayout.setCurrentIndex(2)





if __name__ == '__main__':
    sample = ms.slab(3)
    sample.addlayer(0,'SrTiO3', 50)
    sample.addlayer(1,'LaMnO3', 20)
    sample.addlayer(2, 'LaAlO3', 5)

    app = QApplication(sys.argv)
    demo = ReflectometryApp(sample)

    demo.show()

    sys.exit(app.exec_())