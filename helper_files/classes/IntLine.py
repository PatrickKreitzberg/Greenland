import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
from Marker import *

class intLine():
    def __init__(self, ox, oy, pen, marker):
        self.PlotDataItem = pg.PlotDataItem(ox, oy, pen=whitePlotPen)
        self.marker = marker
        self.PlotDataItem.curve.setClickable(True)
        self.PlotDataItem.curve.opts['mouseWidth'] = 20


    def setFunction(self, fun):
        self.PlotDataItem.sigClicked.connect(fun)



    #def __del__(self):



