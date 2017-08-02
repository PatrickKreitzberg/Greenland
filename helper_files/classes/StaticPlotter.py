from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import pyqtgraph as pg
from pyqtgraph import LegendItem
import PyQt4
from pylab import plot,show,ion,subplots
from dolfin import project
from MyLegend import *
from ..dataset_objects import *
from ..gui import *

def runStaticPlot():
    staticPlotWindow = QtGui.QMainWindow(mw)
    staticPlotWindow.setWindowTitle('Static Plotter')
    dummyWidget = QtGui.QWidget()
    staticPlotWindow.setCentralWidget(dummyWidget)
    layout = QtGui.QVBoxLayout()
    dummyWidget.setLayout(layout)
    plt1 = pg.PlotWidget()
    plt2 = pg.PlotWidget()
    plt3 = pg.PlotWidget()
    layout.addWidget(plt1)
    layout.addWidget(plt2)
    layout.addWidget(plt3)

    surfPlt  = plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=surface.pen)
    bedPlt   = plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bed.pen)
    thickPlt = plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=thickness.pen)


    velocityPlt = plt2.getPlotItem().plot(velocity.distanceData, velocity.pathData, pen=velocity.pen)
    smbPlt = plt3.getPlotItem().plot(smb.distanceData, smb.pathData, pen=smb.pen)
    legend1 = pg.LegendItem(offset=(-1,1))
    legend1.addItem(surfPlt, 'Surface(m)')
    legend1.addItem(bedPlt, 'Bed(m)')
    legend1.addItem(thickPlt, 'Thickness(m)')
    legend1.setParentItem(plt1.getPlotItem())

    legend2 = pg.LegendItem(offset=(-1, 1))
    legend2.addItem(velocityPlt, 'Velocity(m/yr)')
    legend2.setParentItem(plt2.getPlotItem())

    legend3 = pg.LegendItem(offset=(-1, 1))
    legend3.addItem(smbPlt, 'SMB(m)')
    legend3.setParentItem(plt3.getPlotItem())
    staticPlotWindow.show()