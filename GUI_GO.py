from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt

import sys


class sampleWidget(QWidget):
    def __init__(self):
        super().__init__()
        pagelayout = QHBoxLayout()

        cblayout = QVBoxLayout()

        addlayerButton = QPushButton('Add Layer')
        copylayerButton = QPushButton('Copy Layer')
        deletelayerButton = QPushButton('Delete Layer')


        self.layerBox = QComboBox(self)
        self.layerBox.addItems(['Substrate', 'Layer 1', 'Layer 2', 'Layer 3'])

        # buttons
        cblayout.addWidget(addlayerButton)
        cblayout.addWidget(copylayerButton)
        cblayout.addWidget(deletelayerButton)

        # layer combo box
        cblayout.addWidget(self.layerBox)

        # table widget
        # create table
        self.elementTable = QTableWidget(self)
        #self.elementTable.resize(660, 125)
        self.elementTable.setRowCount(3)
        self.elementTable.setColumnCount(6)

        self.elementTable.setHorizontalHeaderLabels(
            ['Element', 'Thickness', 'Density', 'Roughness', 'Linked Roughness', 'Scattering Factor'])

        selectlayout = QVBoxLayout()

        structButton = QPushButton('Structure')
        selectlayout.addWidget(structButton)
        polyButton = QPushButton('Element Variation')
        selectlayout.addWidget(polyButton)
        magButton = QPushButton('Magnetic')
        selectlayout.addWidget(magButton)

        pagelayout.addLayout(cblayout)
        pagelayout.addWidget(self.elementTable)
        pagelayout.addLayout(selectlayout)

        self.setLayout(pagelayout)

class ReflectometryApp(QMainWindow):
    def __init__(self):
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

        _sampleWidget = sampleWidget()  # initialize the sample widget
        sampleButton = QPushButton('Sample')
        sampleButton.clicked.connect(self.activate_tab_1)
        buttonlayout.addWidget(sampleButton)
        self.stackedlayout.addWidget(_sampleWidget)

        reflButton = QPushButton('Reflectivity')
        reflButton.clicked.connect(self.activate_tab_2)
        buttonlayout.addWidget(reflButton)
        self.stackedlayout.addWidget(label2)

        goButton = QPushButton('Global Optimization')
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


# Element table widget

### Garbage code for layer --------------------------------------------------------------------------------------------
class elementTable(QTableWidget):
    def __init__(self):
        super().__init__(3,5)

# Compound table widget
class compoundTable(QTableWidget):
    def __init__(self):
        super().__init__(1, 5)


class TestWindow(QWidget):
    def __init__(self):
        super().__init__()

        # setting title
        self.setWindowTitle('Reflectometry of Quantum Materials (RQM)')

        # setting geometry
        self.setGeometry(200,80,1000,600)

        # calling method
        # creating a combo box widget
        self.layers = QComboBox(self)

        # setting geometry of layers
        self.layers.setGeometry(200, 150, 120, 30)
        self.layers.move(31, 131)
        self.layers.resize(98, 30)

        # this is where the layers list will be initialized based on input sample info
        my_list = ['Substrate', 'Layer 1', 'Layer 2', 'Layer 3']

        self.layers.addItems(my_list)
        self.layers.activated.connect(self.layerProperties)

        self._createButtons()

        self._createElementTable()
        

    def _createButtons(self):
        # button used to add or remove layers
        self.addlayer_button = QPushButton('Add Layer', self)
        self.addlayer_button.resize(100, 32)
        self.addlayer_button.move(30, 40)
        self.addlayer_button.clicked.connect(self._addLayer)

        self.copylayer_button = QPushButton('Copy Layer', self)
        self.copylayer_button.resize(100, 32)
        self.copylayer_button.move(30, 70)
        self.copylayer_button.clicked.connect(self._copyLayer)

        self.removelayer_button = QPushButton('Remove Layer', self)
        self.removelayer_button.resize(100, 32)
        self.removelayer_button.move(30, 100)
        self.removelayer_button.clicked.connect(self._removeLayer)

        self.struct_button = QPushButton('Structural', self)
        self.struct_button.resize(100, 32)
        self.struct_button.move(150, 5)

        self.poly_button = QPushButton('Element Variation', self)
        self.poly_button.resize(100, 32)
        self.poly_button.move(250, 5)

        self.mag_button = QPushButton('Magnetic', self)
        self.mag_button.resize(100, 32)
        self.mag_button.move(350, 5)

        self.elementMode_button = QPushButton('Element', self)
        self.elementMode_button.resize(100, 32)
        self.elementMode_button.move(150, 175)

        self.compoundMode_button = QPushButton('Compound', self)
        self.compoundMode_button.resize(100, 32)
        self.compoundMode_button.move(250, 175)

    def _createElementTable(self):
        # create table
        self.elementTable = QTableWidget(self)
        self.elementTable.move(150, 40)
        self.elementTable.resize(660, 125)

        self.elementTable.setRowCount(3)
        self.elementTable.setColumnCount(6)

        self.elementTable.setHorizontalHeaderLabels(
            ['Element', 'Thickness', 'Density', 'Roughness', 'Linked Roughness', 'Scattering Factor'])

    def _addLayer(self):
        num = self.layers.count()
        if num == 0:
            self.layers.addItem('Substrate')
        else:
            self.layers.addItem('Layer ' + str(num))

    def _removeLayer(self):
        num = self.layers.count()

        if num != 0:
            self.layers.removeItem(num-1)

    def _copyLayer(self):
        num = self.layers.count()
        if num == 0:
            self.layers.addItem('Substrate')
        else:
            self.layers.addItem('Layer ' + str(num))


    def layerProperties(self, index):
            print('Activated: ' + str(index))


if __name__ == '__main__':
    app = QApplication(sys.argv)
    demo = ReflectometryApp()
    demo.show()
    sys.exit(app.exec_())