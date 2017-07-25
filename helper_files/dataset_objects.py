import numpy as np
from classes.dataset import dataset
from gui import *
from pens import *


#####################################################
####         CREATE INITIAL DATA SET(S)          ####
#####################################################
vpts = [] #holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates

velocity = dataset('velocity', bpLegend, greenPlotPen, map=True)
bp.addItem(velocity.pathPlotItem)
iiContainer.addWidget(velocity.plotWidget)
iiContainer.setCurrentWidget(velocity.plotWidget)

smb = dataset('smb', bpLegend, redPlotPen, map=True)
print 'SMB: ', np.amin(smb.data), np.amax(smb.data) #-11493.3860928 6060.80339304
bp.addItem(smb.pathPlotItem)
iiContainer.addWidget(smb.plotWidget)

surface = dataset('surface', bpLegend, greyPlotPen, map=True)
bp.addItem(surface.pathPlotItem)
iiContainer.addWidget(surface.plotWidget)

bed = dataset('bed', bpLegend, bluePlotPen, map=True)
bp.addItem(bed.pathPlotItem)
iiContainer.addWidget(bed.plotWidget)

velocityWidth = dataset('velocitywidth', bpLegend, purplePlotPen)