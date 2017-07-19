# Found the way to load a QImage into the label on GraphicsLayoutWidget.
#
# Cheers,
# Vasilije

import pyqtgraph
from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets


app = pyqtgraph.QtGui.QApplication([])
win = pyqtgraph.GraphicsLayoutWidget()
qi = pyqtgraph.QtGui.QImage()
label = pyqtgraph.QtGui.QLabel(win)

qi.load('v.jpg')
label.setPixmap(pyqtgraph.QtGui.QPixmap.fromImage(qi))

win.show()


if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        pyqtgraph.QtGui.QApplication.instance().exec_()