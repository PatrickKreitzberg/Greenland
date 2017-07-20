import h5py
from mathFunctions import *
from cm import *

class dataset():
    def __init__(self, dataDictName, bpLegend, pen, map=False, dataFileName='./data/AllDataSets.h5', dataCMFileName='./data/dataCMValues.h5',):
        '''
        dataDictNames: bed, surface, SMB_rec
        dataFileName is the name of the hdf5 file with all the data in it
        dataDictName is name of the data inside the dataFileName
        interpParamter is parameter to send to mathFunctions.getInterpolators()
        pen is the pen for the bottom plot legend
        '''
        if dataDictName == 'velocity':
            self.data, self.vx, self.vy = self.setData(dataFileName, dataDictName)
            self.vxInterp, self.vyInterp = getInterpolators(self.vx, dataDictName, self.vy)
        elif dataDictName == 'smb':
            self.data = self.setData(dataFileName, 'SMB_rec')
            # self.interp = getInterpolators(self.data, dataDictName)  # bedInterp
        else:
            self.data = self.setData(dataFileName, dataDictName)
            # self.interp = getInterpolators(self.data, dataDictName)  # bedInterp
        if map:
            self.colorData    = self.setColorData(dataCMFileName, dataDictName)
            self.colorMap     = getCM(dataDictName)
            self.colorBar     = getColorBar(dataDictName, self.colorMap)
            self.imageItem    = pg.ImageItem(self.colorData)
            self.imageItem.setOpts(axisOrder='row-major')

            # Setup plotWidget
            self.plotWidget   = pg.PlotWidget()      # velW
            self.plotWidget.addItem(self.imageItem)
            self.plotWidget.setAspectLocked(True)
            self.plotWidget.invertY(True)

        self.pathPlotItem = pg.PlotDataItem([0,0],  pen=pen)  # bpSurf
        self.legendItem   = bpLegend.addItem(self.pathPlotItem, dataDictName)      # bplSMB
        self.pathData     = None        # nsmb nv etc.


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
            print 'Found color data'
            return data
        else:
            print 'ERROR: ', str(dataDictName + 'CM'), ' Dataset not found in ', dataCMFileName
            datFile.close()
            return 1


# dataFileName='../data/AllDataSets.h5'
# datFile = h5py.File(dataFileName, 'r')
# print datFile.keys()
# datFile.close()

vel_x0 = -638000  # first x coordinate
vel_x1 = 864550  # last x coordinate
vel_y0 = -657600  # first y coordinate
vel_y1 = -3349350  # last y coordinate

print (vel_x1 - vel_x0)/150 + 1
print min(10,11)
# vel_xarray = linspace(vel_x0, vel_x1, 10018, endpoint=True)
# vel_yarray = linspace(vel_y1, vel_y0, 17946, endpoint=True)