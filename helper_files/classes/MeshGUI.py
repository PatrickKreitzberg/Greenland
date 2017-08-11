import pyqtgraph as pg
from pyqtgraph.Qt import QtGui


class MeshGUI(QtGui.QMainWindow):
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        self.cw = QtGui.QWidget()
        self.setCentralWidget(self.cw)
        self.horLayout = QtGui.QHBoxLayout()
        self.cw.setLayout(self.horLayout)
        self.verLayout = QtGui.QVBoxLayout()
        self.pw = pg.PlotWidget()



