from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import pyqtgraph as pg
from pyqtgraph import LegendItem
import PyQt4
from pylab import plot,show,ion,subplots, sqrt
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
    legend1 = pg.LegendItem(offset=(-1, 1))
    legend2 = pg.LegendItem(offset=(-1, 1))
    legend3 = pg.LegendItem(offset=(-1, 1))

    if surfaceCheck.checkState() == 2:
        surfPlt  = plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=surface.pen)
        legend1.addItem(surfPlt, 'Surface(m)')

    if bedCheck.checkState() == 2:
        bedPlt   = plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bed.pen)
        legend1.addItem(bedPlt, 'Bed(m)')

    if thicknessCheck.checkState() == 2:
        thickPlt = plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=thickness.pen)
        legend1.addItem(thickPlt, 'Thickness(m)')

    if velocityCheck.checkState() == 2:
        velocityPlt = plt2.getPlotItem().plot(velocity.distanceData, velocity.pathData, pen=velocity.pen)
        dv = velocity.pathData[1:] - velocity.pathData[:-1]
        dv = sqrt(dv**2)
        dvPlt = plt1.getPlotItem().plot(velocity.distanceData[:-1], dv, pen=smb.pen)
        legend2.addItem(velocityPlt, 'Velocity(m/yr)')

    if smbCheck.checkState() == 2:
        smbPlt = plt3.getPlotItem().plot(smb.distanceData, smb.pathData, pen=smb.pen)
        legend3.addItem(smbPlt, 'SMB(m)')

    legend1.setParentItem(plt1.getPlotItem())

    legend2.setParentItem(plt2.getPlotItem())

    legend3.setParentItem(plt3.getPlotItem())
    staticPlotWindow.show()