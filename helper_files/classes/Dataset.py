import time
import h5py
from helper_files.math_functions import *
from helper_files.colorMaps import *
from helper_files.classes.ColorBarAnchorWidget import ColorBarAnchorWidget
from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline
import numpy as np
from ..gui import *
from ..constants import map

class Dataset():
    def __init__(self, name, pen, draw=False, dataFileName='./data/GreenlandInBedCoord.h5', dataCMFileName='./data/dataCMValues.h5', subSample=spatialRes):
        '''
        names: bed, surface, SMB_rec
        dataFileName is the name of the hdf5 file with all the data in it
        name is name of the data inside the dataFileName
        interpParamter is parameter to send to mathFunctions.getInterpolators()
        pen is the pen for the bottom plot legend
        '''
        spatialRes = subSample
        subSample = 1
        interpSS = 6
        self.name = name
        self.pen = pen
        self.dataCMFileName = dataCMFileName
        if self.name == 'velocity':
            self.data, self.vx, self.vy = self.setData(dataFileName, name)
            bed_xarray = linspace(map['proj_x0'], map['proj_x1'], map['x1'], endpoint=True)
            bed_yarray = linspace(map['proj_y1'], map['proj_y0'], map['y1'], endpoint=True)
            t0 = time.time()
            self.interp   = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.data).transpose())
            self.vxInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vx).transpose())
            self.vyInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.vy).transpose())
            self.createColorMap()
            # print "interp took ", time.time() - t0
        elif self.name == 'velocitywidth':
            self.data = None
        else:
            # subSample=10
            self.data = self.setData(dataFileName, name)
            t0 = time.time()
            map['x1'] = len(self.data[0])
            map['y1'] = len(self.data)
            bed_xarray = linspace(map['proj_x0'], map['proj_x1'], map['x1'], endpoint=True)
            bed_yarray = linspace(map['proj_y1'], map['proj_y0'], map['y1'], endpoint=True)
            self.interp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(self.data).transpose())
            self.createColorMap()
            # print "interp took ", time.time() - t0
        # if draw:
        #     self.colorData = self.setColorData(dataCMFileName, name)
        #
        #     # Setup imageitem
        #     self.imageItem    = pg.ImageItem(self.colorData)
        #     self.imageItem.setOpts(axisOrder='row-major')
        #
        #     # Setup plotWidget
        #     self.plotWidget   = pg.PlotWidget()      # velW
        #     self.plotWidget.addItem(self.imageItem)
        #     self.plotWidget.setAspectLocked(True)
        #     self.plotWidget.invertY(True)
        #     self.colorMap  = getCM(name)
        #     self.colorBar  = getColorBar(name, self.colorMap)
        #
        #     self.colorBarAnchorWidget = ColorBarAnchorWidget()
        #     self.colorBarAnchorWidget.hideAxis('left')
        #     self.colorBarAnchorWidget.hideAxis('bottom')
        #     self.colorBarAnchorWidget.addItem(self.colorBar)
        #
        #     self.plotWidget.addItem(self.colorBarAnchorWidget)
        #     self.colorBarAnchorWidget.setFixedWidth(158)
        #     self.colorBarAnchorWidget.setFixedHeight(292)
        #     self.colorBarAnchorWidget.setAspectLocked(True)
        #     self.colorBarAnchorWidget.getViewBox().setRange(xRange=[-44.0,114], yRange=[-15,247], padding=0.0)
        #     self.colorBarAnchorWidget.invertY(True)
        #     self.colorBarAnchorWidget.setParentItem(self.plotWidget.getPlotItem())
        #     self.colorBarAnchorWidget.getViewBox().setMouseEnabled(x=False, y=False)
        #     self.colorBarAnchorWidget.anchor(itemPos=(1,0), parentPos=(1,0), offset=(-10,-10))

        self.pathPlotItem = pg.PlotDataItem([0,0], pen=self.pen)  # bpSurf
        # self.legendItem   = bpLegend.addItem(self.pathPlotItem, name)      # bplSMB
        self.pathData     = None        # nsmb nv etc.
        self.distanceData = None    # x data for plots.  Which is distance in proj coordinates

    def createColorMap(self):
        # fn = './data/' + self.name + 'CM.h5'
        self.colorMapFile = h5py.File(self.dataCMFileName, 'r')
        self.colorData = self.colorMapFile[self.name][:]

        # Setup imageitem
        self.imageItem = pg.ImageItem(self.colorData)
        self.imageItem.setOpts(axisOrder='row-major')

        # Setup plotWidget
        self.plotWidget = pg.PlotWidget()  # velW
        self.plotWidget.addItem(self.imageItem)
        self.plotWidget.setAspectLocked(True)
        self.plotWidget.invertY(True)
        self.colorMap = getCM(self.name)
        self.colorBar = getColorBar(self.name, self.colorMap)

        self.colorBarAnchorWidget = ColorBarAnchorWidget()
        self.colorBarAnchorWidget.hideAxis('left')
        self.colorBarAnchorWidget.hideAxis('bottom')
        self.colorBarAnchorWidget.addItem(self.colorBar)

        self.plotWidget.addItem(self.colorBarAnchorWidget)
        self.colorBarAnchorWidget.setFixedWidth(158)
        self.colorBarAnchorWidget.setFixedHeight(250)
        self.colorBarAnchorWidget.setAspectLocked(True)
        self.colorBarAnchorWidget.getViewBox().setRange(xRange=[-64.0, 114], yRange=[-15, 247], padding=0.0)
        self.colorBarAnchorWidget.invertY(True)
        self.colorBarAnchorWidget.setParentItem(self.plotWidget.getPlotItem())
        self.colorBarAnchorWidget.getViewBox().setMouseEnabled(x=False, y=False)
        self.colorBarAnchorWidget.anchor(itemPos=(1, 0), parentPos=(1, 0), offset=(-10, -10))
        self.colorMapFile.close()

    def setInterpolator(self, subSample):
        bed_xarray = linspace(map['proj_x0'], map['proj_x1'], map['x1'], endpoint=True)
        bed_yarray = linspace(map['proj_y1'], map['proj_y0'], map['y1'], endpoint=True)
        if self.name == 'velocity':
            self.interp = RectBivariateSpline(bed_xarray, bed_yarray,
                                              np.flipud(self.data).transpose())
        else:
            self.interp = RectBivariateSpline(bed_xarray, bed_yarray,
                                              np.flipud(self.data).transpose())

    def setData(self, dataFileName, dataDictName):
        if dataDictName == 'velocity':
            datFile = h5py.File(dataFileName, 'r')
            vx = datFile['VX'][:]
            vy = datFile['VY'][:]
            data = sqrt(vx ** 2 + vy ** 2)
            datFile.close()
            return data, vx, vy
        else:
            datFile = h5py.File(dataFileName, 'r')

            if dataDictName in datFile.keys():
                data = datFile[dataDictName][:]
                datFile.close()
                return data
            else:
                print 'ERROR Dataset not found in ', dataFileName
                datFile.close()
                return 1

    def setColorData(self, dataCMFileName, dataDictName):
        datFile = h5py.File(dataCMFileName, 'r')
        if dataDictName in datFile.keys():
            data = datFile[dataDictName][:]
            datFile.close()
            return data
        else:
            print 'ERROR: ', str(dataDictName + 'CM'), ' Dataset not found in ', dataCMFileName
            datFile.close()
            return 1


