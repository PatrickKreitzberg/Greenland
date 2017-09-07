import time
startTime = time.time()
import fenics as fc
import sys
from scipy.interpolate import interp1d
from pyqtgraph.Qt import QtCore  # , QtWidgets
from scipy.integrate import ode
from helper_files.classes.Dataset import Dataset
from helper_files.classes.Marker import *
from helper_files.colorMaps import *
from helper_files.gui import *
from helper_files.math_functions import *
from helper_files.pens import *
from helper_files.profile_driver import runModel
from helper_files.dataset_objects import *
from helper_files.data_functions import *
from helper_files.gui_functions import *
from helper_files.constants import *
from helper_files.mesh_functions import *
from helper_files.classes.Instructions import *
from helper_files.classes.ModelGui import ModelGUI

##################################################

#        BRANCHING   BRANCHING   BRANCHING   BRANCHING   BRANCHING

##################################################


'''
https://stackoverflow.com/questions/38065570/pyqtgraph-is-it-possible-to-have-a-imageview-without-histogram

RACMO:

Latitude of the origin:                             90
Longitude of the origin(central meridian):         -45
Standard parallel:                                  70
Ellipsoid:WGS84
Datum:WGS84
'''

print 'Loading...'
#####################################################
####           CONSTANTS/GLOBAL VARIABLES       #####
#####################################################

botPlot = False  # if bottom plot has been populated or not
inBed = False
inSMB = False
inSurface = False
clickedCurve = False
dataCMFileName = './data/dataCMValues.h5'
dataFileName   = './data/GreenlandInBedCoord.h5'
mapList.currentIndexChanged.connect(changeMap)
integrateLine = None

#####################################################
####                    GUI                     #####
#####################################################

'''
    GUI is created by importing the gui.py file.  
        RIGHTSIDE OF GUI
    buttonBoxWidget: Generic QWidget which holds the buttons
    buttonBox: QLayout widget which controls the buttons/their placement
    dataCheckContainer: QWidget which holds the checkboxes which choose what to display on bot plot
        LEFTSIDE OF GUI
    lsw: Generic QWidget which holds the two plots
    leftSide: QLayout widget which controls the plots/their placement
'''


def hiResInterpolators():
    velocity.setInterpolator(1)
    surface.setInterpolator(1)
    smb.setInterpolator(1)
    bed.setInterpolator(1)
    thickness.setInterpolator(1)
    model_res_lineEdit.setText('150')
    hiResButton.setEnabled(False)

def clearPoints():
    del vpts[:]
    textOut.setText('')
    modelButton.setEnabled(False)
    cProfButton.setEnabled(False)
    velocity.pathPlotItem.clear()
    surface.pathPlotItem.clear()
    smb.pathPlotItem.clear()
    bed.pathPlotItem.clear()
    thickness.pathPlotItem.clear()


def createModelGUI():
    m = ModelGUI(mw)


# def runModelButt():
#     dr = float(model_res_lineEdit.text()) # dr = 150
#     if len(vpts) > 0:
#
#         interpolateData(True)
#
#         ###########################################
#         ###  INTERPOLATE DATA SO EVEN INTERVAL  ###
#         ###########################################
#         '''
#         Interpolating data at an even interval so the data points align with the
#         invterval mesh.
#         '''
#         thickness1dInterp = interp1d(thickness.distanceData, thickness.pathData)
#         bed1dInterp       = interp1d(bed.distanceData,       bed.pathData)
#         surface1dInterp   = interp1d(surface.distanceData,   surface.pathData)
#         smb1dInterp       = interp1d(smb.distanceData,       smb.pathData)
#         velocity1dInterp  = interp1d(velocity.distanceData,  velocity.pathData)
#
#         # N is the number of total data points including the last
#         # Data points on interval [0, N*dr] inclusive on both ends
#
#         N = int(np.floor(bed.distanceData[-1]/float(dr))) # length of path / resolution
#
#         x = np.arange(0, (N+1)*dr, dr) # start point, end point, number of segments. END POINT NOT INCLUDED!
#         print 'N, dr, dr*N', N, dr, dr*N
#         mesh = fc.IntervalMesh(N, 0, dr * N)  # number of cells, start point, end point
#
#         thicknessModelData = thickness1dInterp(x)
#         bedModelData       = bed1dInterp(x)
#         surfaceModelData   = surface1dInterp(x)
#         smbModelData       = smb1dInterp(x)
#         velocityModelData  = velocity1dInterp(x)
#
#         THICKLIMIT = 10.  # Ice is never less than this thick
#         H = surfaceModelData - bedModelData
#         surfaceModelData[H <= THICKLIMIT] = bedModelData[H <= THICKLIMIT]
#
#         #FIXME the intervalMesh is consistantly 150 between each datapoint this not true for the data being sent
#         hdf_name = '.data/latest_profile.h5'
#         hfile = fc.HDF5File(mesh.mpi_comm(), hdf_name, "w")
#         V = fc.FunctionSpace(mesh,"CG",1)
#
#         functThickness = fc.Function(V, name="Thickness")
#         functBed       = fc.Function(V, name="Bed")
#         functSurface   = fc.Function(V, name="Surface")
#         functSMB       = fc.Function(V, name='SMB')
#         functVelocity  = fc.Function(V, name='Velocity')
#
#         surface.pathPlotItem.setData(x, surfaceModelData)
#         pg.QtGui.QApplication.processEvents()
#
#         functThickness.vector()[:] = thicknessModelData
#         functBed.vector()[:]       = bedModelData
#         functSurface.vector()[:]   = surfaceModelData
#         functSMB.vector()[:]       = smbModelData
#         functVelocity.vector()[:]  = velocityModelData
#
#         hfile.write(functThickness.vector(), "/thickness")
#         hfile.write(functBed.vector(),       "/bed")
#         hfile.write(functSurface.vector(),   "/surface")
#         hfile.write(functSMB.vector(),       "/smb")
#         hfile.write(functVelocity.vector(),  "/velocity")
#         hfile.write(mesh, "/mesh")
#         hfile.close()
#         runModel(hdf_name)

def showInstructions():
    Instructions(mw)


# proxy = pg.SignalProxy(bp.getPlotItem().scene().sigMouseMoved, rateLimit=60, slot=mouseMovedBP)


vptCur = None

#####################################################
####         CONNECT BUTTONS TO FUNCTIONS       #####
#####################################################

clearButton.clicked.connect(clearPoints)
calcWidthButton.clicked.connect(cwLoop)
cProfButton.clicked.connect(calcBP) #FIXME should change names so calcProf isn't the integration function
cRegionButton.clicked.connect(intLine)
cVelArrowsButton.clicked.connect(arrows)
modelButton.clicked.connect(createModelGUI)
hiResButton.clicked.connect(hiResInterpolators)
meshButton.clicked.connect(meshGui)
instructionButton.clicked.connect(showInstructions)



velocity.imageItem.hoverEvent = mouseMoved
# bed.imageItem.hoverEvent = mouseMoved
# surface.imageItem.hoverEvent = mouseMoved
# thickness.imageItem.hoverEvent = mouseMoved


velocity.imageItem.mouseClickEvent = mouseClick

# proxy = pg.SignalProxy(bp.getPlotItem().scene().sigMouseMoved, rateLimit=60, slot=mouseMovedBP)

#fixme velocity.imageItem.mouseClickEvent = mouseClick # Default map is vel FIXME

iiContainer.keyboardGrabber()
iiContainer.keyPressEvent = ky
iiContainer.keyReleaseEvent = ky


print 'Loaded in: ', time.time() - startTime, ' seconds!'

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()

'''
TO DO 

    GENERAL
- Work on velocity width alg.
- TRY AND GET MESH WORKING!
- Shift+click Marker to integrate from there

    MODEL
- Align extrapolated points with fenics mesh
- Allow modeling down integrated line


NEW
    Model screwed up when set to 150 resolution when it worked at 300
    
    Could change so where the glacier ends is always on right side of plot for the model
    


'''