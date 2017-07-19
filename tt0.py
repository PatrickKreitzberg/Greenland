import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import PyQt4
import h5py

from pylab import sqrt
from helper_files.cm import *


allData = h5py.File('./data/AllDataSets.h5', 'r')
vx = allData['VX'][:]
vy = allData['VY'][:]
v = sqrt(vx**2 + vy**2)
allData.close()
velVals = getCM('velocity').map(v)

app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.setWindowTitle('GREENLAND')    # MAIN WINDOW
cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
mw.setCentralWidget(cw)
l = QtGui.QGridLayout()            # CENTRAL WIDGET LAYOUT (layout of the central widget)
cw.setLayout(l)
p = pg.PlotWidget()
ii = pg.ImageItem(velVals)
ii.setOpts(axisOrder='row-major')
p.addItem(ii)
p.invertY(True)
p.setAspectLocked(True)
l.addWidget(p)
print ii.save('./velII.png')
mw.show()




## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()