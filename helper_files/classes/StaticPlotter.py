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
    # plt3 = pg.PlotWidget()
    layout.addWidget(plt1)
    layout.addWidget(plt2)
    surfPlt  = plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=whitePlotPen)
    bedPlt   = plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bluePlotPen)
    thickPlt = plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=greenPlotPen)
    legend1 = pg.LegendItem()
    legend1.addItem(surfPlt, 'Surface')
    legend1.addItem(bedPlt, 'Bed')
    legend1.addItem(thickPlt, 'Thickness')
    plt1.addItem(legend1)
    staticPlotWindow.show()