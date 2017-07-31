import time
import h5py
from helper_files.math_functions import *
from helper_files.cm import *
from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline
import numpy as np
from ..gui import *

class dataset():
    def __init__(self, name, bpLegend, pen, map=False, dataFileName='./data/GreenlandInBedCoord.h5', dataCMFileName='./data/dataCMValues.h5', subSample=5):
        '''
        names: bed, surface, SMB_rec
        dataFileName is the name of the hdf5 file with all the data in it
        name is name of the data inside the dataFileName
        interpParamter is parameter to send to mathFunctions.getInterpolators()
        pen is the pen for the bottom plot legend
        '''

        bed_x0 = -637925  # first x
        bed_x1 = 864625  # last x

        bed_y0 = -657675
        bed_y1 = -3349425

        bed_xarray = linspace(bed_x0, bed_x1, 10018,
                              endpoint=True)  # FIXME should maybe be one less point?? Prob not because +150 on one side, -150 on the other
        bed_yarray = linspace(bed_y1, bed_y0, 17946, endpoint=True)

        self.name = name
        if self.name == 'velocity':
            self.data, self.vx, self.vy = self.setData(dataFileName, name)
            # self.vxInterp, self.vyInterp = getInterpolators(self.vx, dataDictName, self.vy)
            t0 = time.time()
            self.interp = RectBivariateSpline(bed_xarray[::subSample], bed_yarray[::subSample], np.flipud(self.data[::subSample, ::subSample].transpose()))
            print "interp took ", time.time() - t0
        elif self.name == 'velocitywidth':
            self.data = None
        else:
            self.data = self.setData(dataFileName, name)
            t0 = time.time()
            self.interp = RectBivariateSpline(bed_xarray[::subSample], bed_yarray[::subSample], np.flipud(self.data[::subSample, ::subSample].transpose()))
            print "interp took ", time.time() - t0
        if map:
            self.colorData = self.setColorData(dataCMFileName, name)
            self.colorMap  = getCM(name)
            self.colorBar  = getColorBar(name, self.colorMap)

            # Setup imageitem
            self.imageItem    = pg.ImageItem(self.colorData)
            self.imageItem.setOpts(axisOrder='row-major')

            # Setup plotWidget
            self.plotWidget   = pg.PlotWidget()      # velW
            self.plotWidget.addItem(self.imageItem)
            self.plotWidget.setAspectLocked(True)
            self.plotWidget.invertY(True)
            self.plotWidget.addItem(self.colorBar)
            self.plotWidget.getPlotItem().getViewBox().sigRangeChanged.connect(self.vbRangeChange)
            self.cbScale = 1.1

            # self.colorBar.setParentItem(self.plotWidget.getPlotItem().getViewBox())
            # self.colorBar.anchor(itemPos=(1,0), parentPos=(1,0), offset=(-10,10))
            # self.colorBar.autoAnchor(pos=(1,1))
            # self.colorBar.anchor(itemPos=(0,0), parentPos=(0,0), offset=(-50,-50))
        self.pathPlotItem = pg.PlotDataItem([0,0], pen=pen)  # bpSurf
        self.legendItem   = bpLegend.addItem(self.pathPlotItem, name)      # bplSMB
        self.pathData     = None        # nsmb nv etc.
        self.distanceData = None    # x data for plots.  Which is distance in proj coordinates

    def vbRangeChange(self, e):
        rng = e.state['viewRange'] #[[xmin, xmax], [ymin, ymax]]
        yr = rng[1][1] - rng[1][0]
        xr = rng[0][1] - rng[0][0]
        # print 'INFO'
        # if self.cbScale == 2:
        #     print '1'
        #     self.cbScale = 0.5
        #     self.colorBar.scale(self.cbScale, self.cbScale)
        # else:
        #     print '2'
        #     self.cbScale = 2
        #     self.colorBar.scale(self.cbScale, self.cbScale)

        # print self.colorBar.childrenShape()

        # self.colorBar.sca
        # self.colorBar.mapToParent(self.colorBar.boundingRegion())
        # print yr, xr
        # print self.colorBar.width(), self.colorBar.height()
        # # self.colorBar.setFixedHeight(yr/5)
        #
        #
        # self.colorBar.scale(1.1, 1.1)
        # print self.colorBar.sceneBoundingRect()
        # print self.colorBar.scene().width()


        # self.colorBar.updateGeometry()
        # self.colorBar.setMaximumHeight(yr/5)
        # self.colorBar.setMinimumHeight(yr/5)
        # self.colorBar.setMaximumWidth(xr/6)
        # self.colorBar.setMinimumWidth(xr/6)
        # self.colorBar.update()
        # self.colorBar.sca


    def setData(self, dataFileName, dataDictName):
        if dataDictName == 'velocity':
            datFile = h5py.File(dataFileName, 'r')
            vx = datFile['VX'][:]
            vy = datFile['VY'][:]
            data = sqrt(vx ** 2 + vy ** 2)
            datFile.close()
            return data, vx, vy
        else:
            datFile   = h5py.File(dataFileName, 'r')
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