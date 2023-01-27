import numpy as np
import h5py
import pyqtgraph as pg
from PyQt5.QtWidgets import *
import sys


def hdf5ToDict(hform):

    mydict = dict()
    for key in hform.keys():
        mydict[key] = dict()
        mydict[key]['Data'] = list(hform[key])

        for attrskey,val in hform[key].attrs.items():
            mydict[key][attrskey] = val

    return mydict

def ReadDataHDF5(fname):
    """
    Purpose: Reads in the experimental and simulated data from hdf5 file and then plots spectrum chosen by user
    :param fname: File name
    :return:
    """

    f = h5py.File(fname, 'r')
    experiment = f['Experimental_data']
    simulation = f['Simulated_data']

    RS = experiment['Reflectivity_Scan']
    Rsim = simulation['Reflectivity_Scan']

    ES = experiment['Energy_Scan']
    Esim = simulation['Energy_Scan']

    # Collects data information to print to terminal
    data = list()
    data_dict = dict()
    sim_dict = dict()

    for Rkey in RS.keys():
        mydata = RS[Rkey]
        data_dict[Rkey] = RS[Rkey]
        sim_dict[Rkey] = Rsim[Rkey]
        Rdat = [int(mydata.attrs['DatasetNumber']),'Reflectivity', Rkey]
        data.append(Rdat)

    for Ekey in ES.keys():
        mydata = ES[Ekey]
        data_dict[Ekey] = ES[Ekey]
        sim_dict[Ekey] = Esim[Ekey]
        Edat = [int(mydata.attrs['DatasetNumber']),'Energy', Ekey]
        data.append(Edat)

    # Sorts data in appropriate order
    data = np.array(data)
    sort_idx = np.argsort(data[:,0].astype(int))
    data = data[sort_idx]



    data_dict = hdf5ToDict(data_dict)
    sim_dict = hdf5ToDict(sim_dict)
    f.close()

    return data, data_dict, sim_dict

class DataViewerApp(QMainWindow):
    def __init__(self):
        super().__init__()

        self.data_dict = dict()
        self.energyScans = list()
        self.reflectivityScans = list()

        self.data_dict = dict()
        self.energyScans2 = list()
        self.reflectivityScans2 = list()
        self.reflBool = True



        # set the title
        self.setWindowTitle('Data Viewer')

        # set the geometry of the window
        self.setGeometry(180, 60, 1400, 800)

        # load file

        menuBar = self.menuBar()

        # saving and loading area
        fileMenu = QMenu("&File", self)
        # load a file
        self.loadFile = QAction("&Load", self)
        self.loadFile.triggered.connect(self._loadFile)
        fileMenu.addAction(self.loadFile)

        # load a file
        self.loadFile2 = QAction("&Load 2nd File", self)
        self.loadFile2.triggered.connect(self._loadFile2)
        fileMenu.addAction(self.loadFile2)

        menuBar.addMenu(fileMenu)

        pagelayout = QHBoxLayout()
        buttonLayout = QVBoxLayout()
        # plotting space
        self.plottingSpace = pg.PlotWidget()
        self.plottingSpace.setBackground('w')
        self.plottingSpace.addLegend()

        # separate in terms of energy and reflectometry scans
        self.energyButton = QPushButton('Energy Scan')
        self.energyButton.clicked.connect(self._energyButton)
        self.energyButton.setFixedWidth(300)
        self.energyButton.setFixedHeight(50)

        self.reflButton = QPushButton('Reflectivity Scan')
        self.reflButton.clicked.connect(self._reflectivityButton)
        self.reflButton.setFixedWidth(300)
        self.reflButton.setFixedHeight(50)

        self.energyButton.setStyleSheet('background: grey')
        self.reflButton.setStyleSheet('background: grey')

        self.scans = QComboBox()  # stores the scans
        self.scans.activated.connect(self._plot_scans)

        self.scans2 = QComboBox()  # stores the scans
        self.scans2.activated.connect(self._plot_scans)

        buttonLayout.addWidget(self.reflButton)
        buttonLayout.addWidget(self.energyButton)
        buttonLayout.addWidget(self.scans)
        buttonLayout.addWidget(self.scans2)
        buttonLayout.addStretch(1)

        pagelayout.addLayout(buttonLayout)
        pagelayout.addWidget(self.plottingSpace)

        widget = QWidget()
        widget.setLayout(pagelayout)
        self.setCentralWidget(widget)
        # plot the data

    def _plot_scans(self):
        self.plottingSpace.clear()
        name = self.scans.currentText()
        name2 = self.scans2.currentText()
        if name != None and name != '':
            if self.reflBool:  # reflectivity scan
                pol = self.data_dict[name]['Polarization']
                qz = self.data_dict[name]['Data'][0]
                R = self.data_dict[name]['Data'][2]
                if pol == 'AL' or pol == 'AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    qz = qz[rm_idx]
                    R = R[rm_idx]


                self.plottingSpace.plot(qz, R, pen=pg.mkPen((2, 3), width=2), name=name)
                if pol == 'AL' or pol == 'AC':
                    self.plottingSpace.setLogMode(False, False)
                else:
                    self.plottingSpace.setLogMode(False, True)

            else:  # energy scan
                E = self.data_dict[name]['Data'][3]
                R = self.data_dict[name]['Data'][2]
                self.plottingSpace.setLogMode(False, False)
                self.plottingSpace.plot(E, R, pen=pg.mkPen((2, 3), width=2), name=name)

        if name2 != None and name2 != '':
            if self.reflBool:  # reflectivity scan
                pol = self.data_dict2[name2]['Polarization']
                qz = self.data_dict2[name2]['Data'][0]
                R = self.data_dict2[name2]['Data'][2]
                if pol == 'AL' or pol == 'AC':
                    rm_idx = [i for i in range(len(R)) if R[i] < 4 and R[i] > -4]
                    qz = qz[rm_idx]
                    R = R[rm_idx]


                self.plottingSpace.plot(qz, R, pen=pg.mkPen((0, 3), width=2), name=name)
                if pol == 'AL' or pol == 'AC':
                    self.plottingSpace.setLogMode(False, False)
                else:
                    self.plottingSpace.setLogMode(False, True)

            else:  # energy scan
                E = self.data_dict2[name2]['Data'][3]
                R = self.data_dict2[name2]['Data'][2]
                self.plottingSpace.setLogMode(False, False)
                self.plottingSpace.plot(E, R, pen=pg.mkPen((0, 3), width=2), name=name)

    def _energyButton(self):
        self.energyButton.setStyleSheet('background: cyan')
        self.reflButton.setStyleSheet('background: grey')

        self.scans.blockSignals(True)
        self.scans.clear()
        self.scans.addItems(self.energyScans)
        self.scans.blockSignals(False)

        self.scans2.blockSignals(True)
        self.scans2.clear()
        self.scans2.addItems(self.energyScans2)
        self.scans2.blockSignals(False)

        self.reflBool = False

    def _reflectivityButton(self):
        self.energyButton.setStyleSheet('background: grey')
        self.reflButton.setStyleSheet('background: cyan')

        self.scans.blockSignals(True)
        self.scans.clear()
        self.scans.addItems(self.reflectivityScans)
        self.scans.blockSignals(False)

        self.scans2.blockSignals(True)
        self.scans2.clear()
        self.scans2.addItems(self.reflectivityScans2)
        self.scans2.blockSignals(False)

        self.reflBool = True

    def _loadFile(self):
        # Loading the file
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')  # retrieve file name
        if fname.endswith('.h5'):
            data, data_dict, sim_dict = ReadDataHDF5(fname)
            self.data_dict = data_dict
            self.energyScans = list()  # reset scan names
            self.reflectivityScans = list()
            for key in list(self.data_dict.keys()):
                if 'Angle' in list(self.data_dict[key].keys()):
                    self.energyScans.append(key)
                else:
                    self.reflectivityScans.append(key)
            self.energyButton.setStyleSheet('background: grey')
            self.reflButton.setStyleSheet('background: grey')
        else:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("File Error")
            messageBox.setText("Application can only handle hdf5 file types.")
            messageBox.exec()

    def _loadFile2(self):
        # Loading the file
        fname, _ = QFileDialog.getOpenFileName(self, 'Open File')  # retrieve file name
        if fname.endswith('.h5'):
            data2, data_dict2, sim_dict2 = ReadDataHDF5(fname)
            self.data_dict2 = data_dict2
            self.energyScans2 = list()  # reset scan names
            self.reflectivityScans2 = list()
            for key in list(self.data_dict2.keys()):
                if 'Angle' in list(self.data_dict2[key].keys()):
                    self.energyScans2.append(key)
                else:
                    self.reflectivityScans2.append(key)
            self.energyButton.setStyleSheet('background: grey')
            self.reflButton.setStyleSheet('background: grey')
        else:
            messageBox = QMessageBox()
            messageBox.setWindowTitle("File Error")
            messageBox.setText("Application can only handle hdf5 file types.")
            messageBox.exec()



if __name__ == "__main__":
    app = QApplication(sys.argv)
    demo = DataViewerApp()

    demo.show()

    sys.exit(app.exec_())