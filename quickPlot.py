# Found the way to load a QImage into the label on GraphicsLayoutWidget.
#
# Cheers,
# Vasilije

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import h5py

app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.setWindowTitle('GREENLAND')    # MAIN WINDOW
cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
mw.setCentralWidget(cw)
l = QtGui.QGridLayout()            # CENTRAL WIDGET LAYOUT (layout of the central widget)
cw.setLayout(l)
iiContainer = QtGui.QStackedWidget()    # STACKED WIDGET (inside the layout)

allData = h5py.File('./data/AllDataSets.h5', 'a')
smb = allData['smb'][:]
# allData.__delitem__('smb')
# allData.create_dataset('smb_norm', data=smb)
# allData.create_dataset('smb', data=np.flipud(smb))
allData.close()

pw = pg.PlotWidget()
ii = pg.ImageItem(smb)
ii.setOpts(axisOrder='row-major')
pw.addItem(ii)
pw.invertY(True)
pw.setAspectLocked(True)
iiContainer.addWidget(pw)
l.addWidget(iiContainer)




mw.show()


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        pg.QtGui.QApplication.instance().exec_( )
