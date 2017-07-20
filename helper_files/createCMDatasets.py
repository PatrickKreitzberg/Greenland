import pyqtgraph as pg
from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import PyQt4
import h5py
from pylab import sqrt
from cm import *

allData = h5py.File('../data/AllDataSets.h5', 'r')
vx = allData['VX'][:]
vy = allData['VY'][:]
v = sqrt(vx**2 + vy**2)
bed = allData['bed'][:]
surface = allData['surface'][:]
allData.close()


velVals = getCM('velocity').mapToByte(v)
velVals = velVals.astype(np.uint8)

bedVals = getCM('bed').mapToByte(bed)
bedVals = bedVals.astype(np.uint8)

surfaceVals = getCM('surface').mapToByte(surface)
surfaceVals = surfaceVals.astype(np.uint8)

outF = h5py.File('../data/dataCMValues.h5','w') #719
outF.create_dataset('velocity', data=velVals)
outF.create_dataset('bed', data=bedVals)
outF.create_dataset('surface', data=surfaceVals)
outF.close()
