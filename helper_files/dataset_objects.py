# LOCAL IMPORTS
from classes.Dataset import Dataset
from gui import *
from pens import *
from colorMaps import *
from classes.ColorBarAnchorWidget import *
from constants import *

# EXTERNAL IMPORTS
import numpy as np
import threading
import time
from multiprocessing import *
from scipy.interpolate import RectBivariateSpline
import pyqtgraph as pg

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

smbT0 = time.time()
smb = Dataset('smb', redPlotPen, draw=True)
# iiContainer.addWidget(smb.plotWidget)

bed = Dataset('bed', bluePlotPen, draw=True)
# iiContainer.addWidget(bed.plotWidget)

surface = Dataset('surface', greyPlotPen, draw=True)
# iiContainer.addWidget(surface.plotWidget)

velocityWidth = Dataset('velocitywidth', purplePlotPen)

thickness = Dataset('thickness', orangePlotPen, draw=True)
# iiContainer.addWidget(thickness.plotWidget)

velocity.pathPlotItem.clear()
surface.pathPlotItem.clear()
smb.pathPlotItem.clear()
bed.pathPlotItem.clear()


#######################################################
####         CREATE INTERPOLATORS IN PARRALLEL     ####
#######################################################



print 'loading interps'
int0 = time.time()


def getInterpolators(ds, out_q):
    # ds.setInterpolator()
    outdict = {}
    if ds.name == 'velocity':
        interp = RectBivariateSpline(ds.bed_xarray, ds.bed_yarray, np.flipud(ds.data).transpose())
        vxInterp = RectBivariateSpline(ds.bed_xarray, ds.bed_yarray, np.flipud(ds.vx).transpose())
        vyInterp = RectBivariateSpline(ds.bed_xarray, ds.bed_yarray, np.flipud(ds.vy).transpose())
        outdict[ds.name] =                    interp
        outdict[str('vx' + ds.name)] =        vxInterp
        outdict[str('vy' + ds.name)] = vyInterp
    else:
        interp = RectBivariateSpline(ds.bed_xarray, ds.bed_yarray, np.flipud(ds.data).transpose())
        outdict[ds.name] = interp
    out_q.put(outdict)
    return

#
#
#   https://stackoverflow.com/questions/19828612/python-multiprocessing-setting-class-attribute-value
#
#

out_q = Queue()
jobs = []
p0 = Process(target=getInterpolators, args=(smb, out_q))
# print 'process smb'
jobs.append(p0)
p0.start()


p1 = Process(target=getInterpolators, args=(velocity, out_q))
# print 'process velocity'
jobs.append(p1)
p1.start()

p2 = Process(target=getInterpolators, args=(bed, out_q))
# print 'process velocity'
jobs.append(p2)
p2.start()


p3 = Process(target=getInterpolators, args=(surface, out_q))
# print 'process surface'
jobs.append(p3)
p3.start()

p4 = Process(target=getInterpolators, args=(thickness, out_q))
# print 'process thickness'
jobs.append(p4)
p4.start()
resultdict = {}

for j in jobs:
    resultdict.update(out_q.get())

for j in jobs:
    j.join()

bed.interp = resultdict['bed']
smb.interp = resultdict['smb']
thickness.interp = resultdict['thickness']
surface.interp = resultdict['surface']
velocity.interp = resultdict['velocity']
velocity.vxInterp = resultdict['vxvelocity']
velocity.vyInterp = resultdict['vyvelocity']

print 'Created ALL interpolators in ', time.time() - int0, 'seconds'



#######################################################
#### CREATE COLORMAPS IN PARRALLEL INSIDE A THREAD ####
#######################################################

cmT0 = time.time()
print 'Creating colormaps in parallel'
def createColormaps(ds, out_q):
    print ''
    pg.PlotWidget()
    # od = {}
    # od[ds.name] = ds.p
    # out_q.put(od)


def colorMapThread():
    out_q = Queue()
    jobs = []
    p0 = Process(target=createColormaps, args=(smb, out_q))
    # print 'process smb'
    jobs.append(p0)
    p0.start()


    p1 = Process(target=createColormaps, args=(velocity, out_q))
    # print 'process velocity'
    jobs.append(p1)
    p1.start()

    p2 = Process(target=createColormaps, args=(bed, out_q))
    # print 'process velocity'
    jobs.append(p2)
    p2.start()


    p3 = Process(target=createColormaps, args=(surface, out_q))
    # print 'process surface'
    jobs.append(p3)
    p3.start()

    p4 = Process(target=createColormaps, args=(thickness, out_q))
    # print 'process thickness'
    jobs.append(p4)
    p4.start()

    resultdict = {}

    # for j in jobs:
    #     resultdict.update(out_q.get())

    for j in jobs:
        j.join()
def crcm(ds, cmdata):
    colorData = cmdata[ds.name][:]
    ds.imageItem.setImage(colorData)
    ds.plotWidget.addItem(ds.imageItem)
    ds.plotWidget.setAspectLocked(True)
    ds.plotWidget.invertY(True)
    ds.colorMap = getCM(ds.name)
    ds.colorBar = getColorBar(ds.name, ds.colorMap)
    ds.colorBarAnchorWidget = ColorBarAnchorWidget()
    ds.colorBarAnchorWidget.hideAxis('left')
    ds.colorBarAnchorWidget.hideAxis('bottom')
    ds.colorBarAnchorWidget.addItem(ds.colorBar)
    ds.plotWidget.addItem(ds.colorBarAnchorWidget)
    ds.colorBarAnchorWidget.setFixedWidth(158)
    ds.colorBarAnchorWidget.setFixedHeight(250)
    ds.colorBarAnchorWidget.setAspectLocked(True)
    ds.colorBarAnchorWidget.getViewBox().setRange(xRange=[-64.0, 114], yRange=[-15, 247], padding=0.0)
    ds.colorBarAnchorWidget.invertY(True)
    ds.colorBarAnchorWidget.setParentItem(ds.plotWidget.getPlotItem())
    ds.colorBarAnchorWidget.getViewBox().setMouseEnabled(x=False, y=False)
    ds.colorBarAnchorWidget.anchor(itemPos=(1, 0), parentPos=(1, 0), offset=(-10, -10))

def createCMInThread():
    cmdata = h5py.File(fileDict['color'], 'r')
    # crcm(thickness, cmdata)
    # crcm(smb, cmdata)
    # crcm(bed, cmdata)
    # crcm(thickness, cmdata)
    jobs = []
    p0 = Process(target=crcm, args=(smb, cmdata))
    # print 'process smb'
    jobs.append(p0)
    p0.start()

    p1 = Process(target=crcm, args=(velocity, cmdata))
    # print 'process velocity'
    jobs.append(p1)
    p1.start()

    p2 = Process(target=crcm, args=(bed, cmdata))
    # print 'process velocity'
    jobs.append(p2)
    p2.start()

    p3 = Process(target=crcm, args=(surface, cmdata))
    # print 'process surface'
    jobs.append(p3)
    p3.start()

    p4 = Process(target=crcm, args=(thickness, cmdata))
    # print 'process thickness'
    jobs.append(p4)
    p4.start()

    for j in jobs:
        j.join()

    cmdata.close()
    # thickness.createColorMap()
    iiContainer.addWidget(thickness.plotWidget)
    # surface.createColorMap()
    # smb.createColorMap()
    # bed.createColorMap()

# cmThread = threading.Thread(target=createCMInThread)
# cmThread.start()
createCMInThread()
# colorMapThread()
print 'Created ALL Colormaps in', time.time()-cmT0, 'seconds'


print 'Done loading data in', time.time()-dataT0, 'seconds'
