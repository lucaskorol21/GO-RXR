from PyQt5.QtWidgets import *
from PyQt5 import QtCore, QtGui
from PyQt5.QtCore import Qt
import sys

# Element table widget
class elementTable(QTableWidget):
    def __init__(self):
        super().__init__(3,5)

# Compound table widget
class compoundTable(QTableWidget):
    def __init__(self):
        super().__init__(1, 5)

class MainWindow(QMainWindow):
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
        self.layers.move(31, 101)
        self.layers.resize(98, 30)

        # this is where the layers list will be initialized based on input sample info
        my_list = ['Substrate', 'Layer 1', 'Layer 2', 'Layer 3']

        self.layers.addItems(my_list)
        self.layers.activated.connect(self.layerProperties)

        self.createButtons()

        # create table
        self.elementTable = QTableWidget(self)
        self.elementTable.move(150, 40)
        self.elementTable.resize(660,125)

        self.elementTable.setRowCount(3)
        self.elementTable.setColumnCount(6)

        self.elementTable.setHorizontalHeaderLabels(['Element', 'Thickness', 'Density', 'Roughness', 'Linked Roughness', 'Scattering Factor'])


    def createButtons(self):
        # button used to add or remove layers
        self.addlayer_button = QPushButton('Add Layer', self)
        self.addlayer_button.resize(100, 32)
        self.addlayer_button.move(30, 40)
        self.addlayer_button.clicked.connect(self._addLayer)

        self.removelayer_button = QPushButton('Remove Layer', self)
        self.removelayer_button.resize(100, 32)
        self.removelayer_button.move(30, 70)
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


    def layerProperties(self, index):
            print('Activated: ' + str(index))


if __name__ == '__main__':
    app = QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())