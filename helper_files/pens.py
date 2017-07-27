import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

'''
Creates the pens for plotting.
'''

#####################################################
####         CREATE PENS                         ####
#####################################################
plotPen = QtGui.QPen()
plotPen.setWidth(2)
plotPen2       = pg.mkPen(color=(255, 255, 255), width=2)
plotPen3       = pg.mkPen(color=(0,     0,   0), width=2)
plotPen4       = pg.mkPen(color=(101,   4,   4), width=2)
greyPlotPen    = pg.mkPen(color=(200, 200, 200), width=2)
redPlotPen     = pg.mkPen(color=(100,   0,   0), width=2)
bluePlotPen    = pg.mkPen(color=(  0,   0, 255), width=2)
greenPlotPen   = pg.mkPen(color=( 76, 153,   0), width=2)
purplePlotPen  = pg.mkPen(color=(102,   0, 204), width=2)
orangePlotPen  = pg.mkPen(color=(255, 128,   0), width=2)