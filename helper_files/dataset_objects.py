import numpy as np
from classes.dataset import dataset
from gui import *
from pens import *
from cm import *
from classes.cbanchor import *


#####################################################
####         CREATE INITIAL DATA SET(S)          ####
#####################################################
print 'Creating data sets and loading their data/plots'

vpts = [] #holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates

velocity = dataset('velocity', bpLegend, greenPlotPen, map=True)
bp.addItem(velocity.pathPlotItem)
iiContainer.addWidget(velocity.plotWidget)
iiContainer.setCurrentWidget(velocity.plotWidget)
velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=[0,10018], yRange=[0,17946])

smb = dataset('smb', bpLegend, redPlotPen, map=True)
bp.addItem(smb.pathPlotItem)
iiContainer.addWidget(smb.plotWidget)

bed = dataset('bed', bpLegend, bluePlotPen, map=True)
bp.addItem(bed.pathPlotItem)
iiContainer.addWidget(bed.plotWidget)

surface = dataset('surface', bpLegend, greyPlotPen, map=True)
bp.addItem(surface.pathPlotItem)
iiContainer.addWidget(surface.plotWidget)

velocityWidth = dataset('velocitywidth', bpLegend, purplePlotPen)
bp.addItem(velocityWidth.pathPlotItem)

thickness = dataset('thickness', bpLegend, orangePlotPen, map=True)#, map=True)
bp.addItem(thickness.pathPlotItem)
iiContainer.addWidget(thickness.plotWidget)

velocity.pathPlotItem.clear()
surface.pathPlotItem.clear()
smb.pathPlotItem.clear()
bed.pathPlotItem.clear()


colorMap  = getCM('velocity')
colorBar  = getColorBar('velocity', colorMap)
cba = cbAnchor()
# cba.hideAxis('left')
# cba.hideAxis('bottom')
cba.addItem(colorBar)
velocity.plotWidget.addItem(cba)
cba.setFixedWidth(100)
cba.setFixedHeight(400)
cba.setAspectLocked(True)
cba.getViewBox().setRange(xRange=[0,50], yRange=[0,400])

# cba.invertX(True)
cba.invertY(True)
cba.setParentItem(velocity.plotWidget.getPlotItem())
# cba.getViewBox().setMouseEnabled(x=False, y=False)
cba.anchor(itemPos=(1,0), parentPos=(1,0), offset=(-10,10))
print 'CBA Mouse Enabled: ', cba.getViewBox().mouseEnabled()



print 'Done loading data'