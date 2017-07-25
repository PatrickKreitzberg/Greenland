# Found the way to load a QImage into the label on GraphicsLayoutWidget.
#
# Cheers,
# Vasilije

import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import h5py
import pyqtgraph.examples

#
# app = QtGui.QApplication([])
# mw = QtGui.QMainWindow()
# mw.setWindowTitle('GREENLAND')    # MAIN WINDOW
# cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
# mw.setCentralWidget(cw)
# l = QtGui.QGridLayout()            # CENTRAL WIDGET LAYOUT (layout of the central widget)
# cw.setLayout(l)
# iiContainer = QtGui.QStackedWidget()    # STACKED WIDGET (inside the layout)
#
# allData = h5py.File('./data/AllDataSets.h5', 'a')
# smb = allData['smb'][:]
# # allData.__delitem__('smb')
# # allData.create_dataset('smb_norm', data=smb)
# # allData.create_dataset('smb', data=np.flipud(smb))
# allData.close()
#
# pw = pg.PlotWidget()
# ii = pg.ImageItem(smb)
# ii.setOpts(axisOrder='row-major')
# pw.addItem(ii)
# pw.invertY(True)
# pw.setAspectLocked(True)
# iiContainer.addWidget(pw)
# l.addWidget(iiContainer)
#
#
#
#
# mw.show()
#
#
# if __name__ == '__main__':
#     import sys
#     if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
#         pg.QtGui.QApplication.instance().exec_( )





#QtGui.QApplication.setGraphicsSystem('raster')

#mw.resize(800,800)

# win = pg.GraphicsWindow(title="Basic plotting examples")
# win.resize(1000,600)
# win.setWindowTitle('pyqtgraph example: Plotting')

# Enable antialiasing for prettier plots

# x = np.arange(1000)
# y = np.random.normal(size=(3, 1000))
# plotWidget = pg.plot(title="Three plot curves")
# for i in range(3):
#     plotWidget.plot(x, y[i], pen=(i,3))  ## setting pen=(i,3) automaticaly creates three different-colored pens
#
# plotPen2       = pg.mkPen(color=(255, 255, 255), width=2)
# ax = pg.AxisItem('left', pen=plotPen2)
# ax.setLabel('bleh')
#
# plotWidget.getPlotItem().getViewBox().addItem(ax)

from pyqtgraph.Qt import QtGui, QtCore
import numpy as np
import pyqtgraph as pg
app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
win = pg.PlotWidget()
x = np.arange(0, 2*np.pi, 0.1)
y = np.sin(x)
l = pg.PlotDataItem(x=x,y=y)
win.addItem(l)
cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
lay = QtGui.QGridLayout()
cw.setLayout(lay)
lay.addWidget(win)
l.rotate(90)
win.getPlotItem().invertX(True)

ax = win.getPlotItem().getAxis('left')
ax.setLabel('Y')
ax = win.getPlotItem().getAxis('bottom')
ax.setLabel('X')
mw.setCentralWidget(cw)
mw.show()
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()

#
# p2 = win.addPlot(title="Multiple curves")
# p2.plot(np.random.normal(size=100), pen=(255,0,0), name="Red curve")
# p2.plot(np.random.normal(size=110)+5, pen=(0,255,0), name="Green curve")
# p2.plot(np.random.normal(size=120)+10, pen=(0,0,255), name="Blue curve")
#
# p3 = win.addPlot(title="Drawing with points")
# p3.plot(np.random.normal(size=100), pen=(200,200,200), symbolBrush=(255,0,0), symbolPen='w')
#
#
# win.nextRow()
#
# p4 = win.addPlot(title="Parametric, grid enabled")
# x = np.cos(np.linspace(0, 2*np.pi, 1000))
# y = np.sin(np.linspace(0, 4*np.pi, 1000))
# p4.plot(x, y)
# p4.showGrid(x=True, y=True)
#
# p5 = win.addPlot(title="Scatter plot, axis labels, log scale")
# x = np.random.normal(size=1000) * 1e-5
# y = x*1000 + 0.005 * np.random.normal(size=1000)
# y -= y.min()-1.0
# mask = x > 1e-15
# x = x[mask]
# y = y[mask]
# p5.plot(x, y, pen=None, symbol='t', symbolPen=None, symbolSize=10, symbolBrush=(100, 100, 255, 50))
# p5.setLabel('left', "Y Axis", units='A')
# p5.setLabel('bottom', "Y Axis", units='s')
# p5.setLogMode(x=True, y=False)
#
# p6 = win.addPlot(title="Updating plot")
# curve = p6.plot(pen='y')
# data = np.random.normal(size=(10,1000))
# ptr = 0
# def update():
#     global curve, data, ptr, p6
#     curve.setData(data[ptr%10])
#     if ptr == 0:
#         p6.enableAutoRange('xy', False)  ## stop auto-scaling after the first data set is plotted
#     ptr += 1
# timer = QtCore.QTimer()
# timer.timeout.connect(update)
# timer.start(50)
#
#
# win.nextRow()
#
# p7 = win.addPlot(title="Filled plot, axis disabled")
# y = np.sin(np.linspace(0, 10, 1000)) + np.random.normal(size=1000, scale=0.1)
# p7.plot(y, fillLevel=-0.3, brush=(50,50,200,100))
# p7.showAxis('bottom', False)
#
#
# x2 = np.linspace(-100, 100, 1000)
# data2 = np.sin(x2) / x2
# p8 = win.addPlot(title="Region Selection")
# p8.plot(data2, pen=(255,255,255,200))
# lr = pg.LinearRegionItem([400,700])
# lr.setZValue(-10)
# p8.addItem(lr)
#
# p9 = win.addPlot(title="Zoom on selected region")
# p9.plot(data2)
# def updatePlot():
#     p9.setXRange(*lr.getRegion(), padding=0)
# def updateRegion():
#     lr.setRegion(p9.getViewBox().viewRange()[0])
# lr.sigRegionChanged.connect(updatePlot)
# p9.sigXRangeChanged.connect(updateRegion)
# updatePlot()

## Start Qt event loop unless running in interactive mode or using pyside.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()


