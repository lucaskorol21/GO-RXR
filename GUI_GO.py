from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt
import sys

class TableWidget(QTableWidget):
    def __init__(self):
        super().__init__(1, 5)
        self.setHorizontalHeaderLabels(['Formula','Thickness (A)','Density (g/cm^3)','Roughness (A)', 'Linked Roughness'])
        self.verticalHeader().setDefaultSectionSize(25)
        self.horizontalHeader().setDefaultSectionSize(150)
    def _addRow(self):
        rowCount = self.rowCount()
        self.insertRow(rowCount)
    def _removeRow(self):
        if self.rowCount() > 0:
            self.removeRow(self.rowCount()-1)
    def _copyRow(self):
        self.insertRow(self.rowCount())
        rowCount = self.rowCount()
        columnCount = self.columnCount()

        for j in range(columnCount):
            if not self.item(rowCount-2,j) is None:
                self.setItem(rowCount-1,j,QTableWidgetItem(self.item(rowCount-2,j).text()))

class App(QWidget):
    def __init__(self):
        super().__init__()
        self.resize(1600, 600)

        mainLayout = QHBoxLayout()

        table = TableWidget()
        mainLayout.addWidget(table)


        buttonLayout = QVBoxLayout()

        button_new = QPushButton("New")
        button_new.clicked.connect(table._addRow)
        buttonLayout.addWidget(button_new)

        button_remove = QPushButton("Remove")
        button_remove.clicked.connect(table._removeRow)
        buttonLayout.addWidget(button_remove)

        button_copy = QPushButton("Copy")
        button_copy.clicked.connect(table._copyRow)
        buttonLayout.addWidget(button_copy, alignment=Qt.AlignTop)

        mainLayout.addLayout(buttonLayout)
        self.setLayout(mainLayout)
if __name__ == '__main__':
    app = QApplication(sys.argv)
    app.setStyleSheet('QPushButton{font-size: 20px; width: 200px; height: 50px}')
    demo = App()
    demo.show()
    sys.exit(app.exec_())