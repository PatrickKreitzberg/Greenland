import numpy as np
from classes.Dataset import Dataset
from gui import *
from pens import *
from colorMaps import *
from classes.Colorbar_Anchor_Widget import *
from constants import *
import time


#####################################################
####         CREATE INITIAL DATA SET(S)          ####
#####################################################
print 'Creating data sets and loading their data/plots'
dataT0 = time.time()
dataSet = h5py.File('./data/GreenlandInBedCoord.h5', 'r')
map['x1'] = len(dataSet['bed'][:][0])
map['y1'] = len(dataSet['bed'][:])
map['proj_x1'] = dataSet['x'][:][-1]
map['proj_y1'] = dataSet['y'][:][-1]
dataSet.close()

vpts = []     # holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates
intLines = [] # holds integration lines
velT0 = time.time()
velocity = Dataset('velocity', greenPlotPen, draw=True)
iiContainer.addWidget(velocity.plotWidget)
iiContainer.setCurrentWidget(velocity.plotWidget)
velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=[0, 10018], yRange=[0, 17964], padding=0.1)
print 'Loaded velocity object in ', time.time() - velT0, 'seconds'

smbT0 = time.time()
smb = Dataset('smb', redPlotPen, draw=True)
iiContainer.addWidget(smb.plotWidget)
print 'Loaded SMB in', time.time() - smbT0,' seconds'

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

print 'Done loading data in', time.time()-dataT0, 'seconds'