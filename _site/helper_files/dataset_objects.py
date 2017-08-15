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

dataSet = h5py.File('./data/GreenlandInBedCoord.h5', 'r')
map['x1'] = len(dataSet['bed'][:][0])
map['y1'] = len(dataSet['bed'][:])
map['proj_x1'] = dataSet['x'][:][-1]
map['proj_y1'] = dataSet['y'][:][-1]
dataSet.close()

vpts = []     # holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates
intLines = [] # holds integration lines

velocity = Dataset('velocity', greenPlotPen, draw=True)
iiContainer.addWidget(velocity.plotWidget)
iiContainer.setCurrentWidget(velocity.plotWidget)
velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=[0, 10018], yRange=[0, 17964], padding=0.1)

smb = Dataset('smb', redPlotPen, draw=True)
iiContainer.addWidget(smb.plotWidget)

bed = Dataset('bed', bluePlotPen, draw=True)
iiContainer.addWidget(bed.plotWidget)

surface = Dataset('surface', greyPlotPen, draw=True)
iiContainer.addWidget(surface.plotWidget)

velocityWidth = Dataset('velocitywidth', purplePlotPen)

thickness = Dataset('thickness', orangePlotPen, draw=True)
iiContainer.addWidget(thickness.plotWidget)

velocity.pathPlotItem.clear()
surface.pathPlotItem.clear()
smb.pathPlotItem.clear()
bed.pathPlotItem.clear()

print 'Done loading data'