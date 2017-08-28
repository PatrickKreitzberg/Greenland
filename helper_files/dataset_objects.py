import numpy as np
from classes.Dataset import Dataset
from gui import *
from pens import *
from colorMaps import *
from classes.ColorBarAnchorWidget import *
from constants import *
import time
from multiprocessing import *
from scipy.interpolate import RectBivariateSpline

'''
#####################################################
####         CREATE INITIAL DATA SET(S) WITH PARALLEL          ####
#####################################################
print 'Creating data sets and loading their data/plots'
dataT0 = time.time()
dataSet = h5py.File('./data/GreenlandInBedCoord.h5', 'r')
map['x1'] = len(dataSet['bed'][:][0])
map['y1'] = len(dataSet['bed'][:])
map['proj_x1'] = dataSet['x'][:][-1]
map['proj_y1'] = dataSet['y'][:][-1]
dataSet.close()

out_q = Queue()
jobs = []
resultdict = {}
datasets = {}
dataFile = h5py.File('./data/GreenlandInBedCoord.h5', 'r')
def initDS(name, pen, drw, datasets, dataFile):
    print 'initDS'
    datasets[name] = None
    print 'ds keys', datasets.keys()
    ds = Dataset(name, pen, dataFile, draw=drw)
    datasets[name] = ds

    return


print 'starting processes'
p0 = Process(target=initDS, args=('smb', redPlotPen, True, datasets, dataFile))
# print 'process smb'
jobs.append(p0)
p0.start()

p1 = Process(target=initDS, args=('velocity', greenPlotPen, True, datasets, dataFile))
# print 'process velocity'
jobs.append(p1)
p1.start()

p2 = Process(target=initDS, args=('bed', bluePlotPen, True, datasets, dataFile))
# print 'process velocity'
jobs.append(p2)
p2.start()


p3 = Process(target=initDS, args=('surface', greyPlotPen, True, datasets, dataFile))
# print 'process surface'
jobs.append(p3)
p3.start()

p4 = Process(target=initDS, args=('thickness', orangePlotPen, True, datasets, dataFile))
# print 'process thickness'
jobs.append(p4)
p4.start()

# for j in jobs:
#     print 'for j in '
#     resultdict.update(out_q.get())

for j in jobs:
    print 'for j in 2'
    j.join()

print datasets.keys()
dataFile.close()
# vpts = []     # holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates
# intLines = [] # holds integration lines
# velT0 = time.time()
# velocity = datasets['velocity']#Dataset('velocity', greenPlotPen, draw=True)
# iiContainer.addWidget(velocity.plotWidget)
# iiContainer.setCurrentWidget(velocity.plotWidget)
# velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=[0, 10018], yRange=[0, 17964], padding=0.1)
# print 'Loaded velocity object in ', time.time() - velT0, 'seconds'
#
# smbT0 = time.time()
# smb = datasets['smb']#Dataset('smb', redPlotPen, draw=True)
# iiContainer.addWidget(smb.plotWidget)
# print 'Loaded SMB in', time.time() - smbT0,' seconds'
#
# bed = datasets['bed']#Dataset('bed', bluePlotPen, draw=True)
# iiContainer.addWidget(bed.plotWidget)
#
# surface = datasets['surface']#Dataset('surface', greyPlotPen, draw=True)
# iiContainer.addWidget(surface.plotWidget)
#
# velocityWidth = datasets['velocitywidth']#Dataset('velocitywidth', purplePlotPen)
#
# thickness = datasets['thickness']#Dataset('thickness', orangePlotPen, draw=True)
# iiContainer.addWidget(thickness.plotWidget)
#
# velocity.pathPlotItem.clear()
# surface.pathPlotItem.clear()
# smb.pathPlotItem.clear()
# bed.pathPlotItem.clear()

print 'Done loading data in', time.time()-dataT0, 'seconds'


'''
#####################################################
####         CREATE INITIAL DATA SET(S) NO PARRALLEL         ####
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
print 'loading interps'
int0 = time.time()
def getInterpolators(ds):
    ds.setInterpolator()

#
#
#   https://stackoverflow.com/questions/19828612/python-multiprocessing-setting-class-attribute-value
#
#
jobs = []
p0 = Process(target=getInterpolators, args=(smb,))
# print 'process smb'
jobs.append(p0)
p0.start()



p1 = Process(target=getInterpolators, args=(velocity,))
# print 'process velocity'
jobs.append(p1)
p1.start()

p2 = Process(target=getInterpolators, args=(bed,))
# print 'process velocity'
jobs.append(p2)
p2.start()


p3 = Process(target=getInterpolators, args=(surface,))
# print 'process surface'
jobs.append(p3)
p3.start()

p4 = Process(target=getInterpolators, args=(thickness,))
# print 'process thickness'
jobs.append(p4)
p4.start()

# for j in jobs:
#     print 'for j in '
#     resultdict.update(out_q.get())

for j in jobs:
    print 'for j in 2'
    j.join()
print 'fart:', smb.interp
print 'loaded interps', time.time() - int0

print 'Done loading data in', time.time()-dataT0, 'seconds'
