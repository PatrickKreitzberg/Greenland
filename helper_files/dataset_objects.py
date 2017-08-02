import numpy as np
from classes.Dataset import Dataset
from gui import *
from pens import *
from colorMaps import *
from classes.ColorBarAnchorWidget import *
from constants import *


#####################################################
####         CREATE INITIAL DATA SET(S)          ####
#####################################################
print 'Creating data sets and loading their data/plots'

vpts = [] #holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates

velocity = Dataset('velocity', greenPlotPen, map=True)
# bp.addItem(velocity.pathPlotItem)
iiContainer.addWidget(velocity.plotWidget)
iiContainer.setCurrentWidget(velocity.plotWidget)
velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=[0,10018], yRange=[0,17946], padding=0.1)

smb = Dataset('smb', redPlotPen, map=True)
# bp.addItem(smb.pathPlotItem)
iiContainer.addWidget(smb.plotWidget)

bed = Dataset('bed', bluePlotPen, map=True)
# bp.addItem(bed.pathPlotItem)
iiContainer.addWidget(bed.plotWidget)

surface = Dataset('surface', greyPlotPen, map=True)
# bp.addItem(surface.pathPlotItem)
iiContainer.addWidget(surface.plotWidget)

velocityWidth = Dataset('velocitywidth', purplePlotPen)
# bp.addItem(velocityWidth.pathPlotItem)

thickness = Dataset('thickness', orangePlotPen, map=True)
print 'thickness range: ', np.amin(thickness.data), np.amax(thickness.data)
# bp.addItem(thickness.pathPlotItem)
iiContainer.addWidget(thickness.plotWidget)

velocity.pathPlotItem.clear()
surface.pathPlotItem.clear()
smb.pathPlotItem.clear()
bed.pathPlotItem.clear()





print 'Done loading data'