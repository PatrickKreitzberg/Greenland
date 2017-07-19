import fenics as fc
import numpy as np
from scipy.integrate import ode
from pylab import sqrt,linspace,array,argmax, meshgrid, argmin
from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import pyqtgraph as pg
import PyQt4
from helper_files.mathFunctions import *
import time
import h5py
from helper_files.pltPoints import *
from helper_files.cm import *
from helper_files.profile_driver import runModel
import sys
sys.path.append("/home/pat/")
'''
https://stackoverflow.com/questions/38065570/pyqtgraph-is-it-possible-to-have-a-imageview-without-histogram

'''
#
clickedCurve = False

# Import all data from hdf5 file
print 'Loading data'
allData = h5py.File('/home/pat/research/Data/AllDataSets.h5', 'r')
bed = allData['bed'][:]
surfaceVals = allData['surface'][:]
print 'Surface min, max: ', np.amin(surfaceVals), np.amax(surfaceVals)
vx = allData['VX'][:]
vy = allData['VY'][:]
v = sqrt(vx**2 + vy**2)
smb = allData['smb'][:]
allData.close()
print 'v interp'
# vin01 = [vx, 'velocity', vy]
# vxInterp, vyInterp = getInterpolators(vx, 'velocity', vy)

vpts = [] #holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates
print 'Done'


#####################################################
####           SET COLOR MAPS                    ####
#####################################################

print 'Processing data and colormap\nAdding plot'
velVals = getCM('velocity').map(v)
bedVals = getBedCM().map(bed)
surfCMVals = getSurfaceCM().map(surfaceVals)
print 'Done'

#####################################################
####           CONSTANTS                         ####
#####################################################
currentMap = 0  # selects which data map to show [velocity, bed, surface]
dr = 150         # spatial resolution for bottom plot interpolation in meters
botPlot = False  # if bottom plot has been populated or not


#####################################################
####           BUILD GUI                         ####
#####################################################

print 'Building GUI'
#Main window
app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.setWindowTitle('GREENLAND')    # MAIN WINDOW


#dumby widget
cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
mw.setCentralWidget(cw)
l = QtGui.QGridLayout()            # CENTRAL WIDGET LAYOUT (layout of the central widget)
cw.setLayout(l)
iiContainer = QtGui.QStackedWidget()    # STACKED WIDGET (inside the layout)
bedVisible = False
velVisible = True
l.setRowStretch(0,2)
velW = pg.PlotWidget()
iiContainer.addWidget(velW)
iiVel = pg.ImageItem(velVals)        # IMAGE ITEM (inside a plot widget)

iiVel.setOpts(axisOrder='row-major')
velW.addItem(iiVel)
vcb = getVelocityColorbar()
velW.addItem(vcb)
velW.invertY(True)
velW.setAspectLocked(True)

# Time to load
# without bedw/surfaceW: 65 seconds
# with bedw/surfacew:    65 seconds

bedW = pg.PlotWidget()
surfaceW = pg.PlotWidget()
iiContainer.addWidget(bedW)
iiContainer.addWidget(surfaceW)

iiBed = pg.ImageItem(bedVals)        # IMAGE ITEM (inside a plot widget)
iiBed.setOpts(axisOrder='row-major')
bcb = getBedColorbar()
bedW.addItem(iiBed)
bedW.addItem(bcb)

iiSurface = pg.ImageItem(surfCMVals)
iiSurface.setOpts(axisOrder='row-major')
scb = getSurfaceColorbar()
surfaceW.addItem(iiSurface)
surfaceW.addItem(scb)

bedW.invertY(True)
surfaceW.invertY(True)

bedW.setAspectLocked(True)
surfaceW.setAspectLocked(True)

#####################################################
####      SIDE WIDGET WITH BUTTONS               ####
#####################################################

bbW = QtGui.QWidget()
buttonBox = QtGui.QVBoxLayout()
bbW.setLayout(buttonBox)

mapList = QtGui.QComboBox()
showBedButton = QtGui.QPushButton('Show Bed Data')
showVelButton = QtGui.QPushButton('Show Velocity Data')
maps = ['Velocity', 'Bed', 'Surface']
mapList.addItems(maps)

clearButton = QtGui.QPushButton('Clear Points')
calcWidthButton = QtGui.QPushButton('Calculate Velocity Width')
intButton = QtGui.QPushButton('Integrate')
cProfButton = QtGui.QPushButton('CalcProfile') #FIXME should automatically get profile
cRegionButton = QtGui.QPushButton('Region')
cVelArrowsButton = QtGui.QPushButton('Arrows')
modelButton = QtGui.QPushButton('Run Model')

buttonBox.addWidget(mapList)
buttonBox.addWidget(clearButton)
buttonBox.addWidget(calcWidthButton)
buttonBox.addWidget(intButton)
buttonBox.addWidget(cProfButton)
buttonBox.addWidget(cRegionButton)
buttonBox.addWidget(cVelArrowsButton)
buttonBox.addWidget(modelButton)
velTF = QtGui.QTextBrowser()
velTF.setText('X\tY\tVelocity\tXp\tYp')
buttonBox.addWidget(velTF)

print 'Done'
#####################################################
####         CREATE PENS                         ####
#####################################################
plotPen = QtGui.QPen()
plotPen.setWidth(2)
plotPen2     = pg.mkPen(color=(255, 255, 255), width=2)
plotPen3     = pg.mkPen(color=(0,     0,   0), width=2)
plotPen4     = pg.mkPen(color=(101,   4,   4), width=2)
greyPlotPen  = pg.mkPen(color=(200, 200, 200), width=2)
redPlotPen   = pg.mkPen(color=(100,   0,   0), width=2)
bluePlotPen  = pg.mkPen(color=(  0,   0, 255), width=2)
greenPlotPen = pg.mkPen(color=( 76, 153,   0), width=2)



#####################################################
####         CREATE BOTTOM PLOT                  ####
#####################################################
bp = pg.PlotWidget()
bpLegend = bp.getPlotItem().addLegend()
dataLen = 0

bpv, bpSMB, bpvw, bpSurf  = pg.PlotDataItem([0,0],  pen=greyPlotPen), pg.PlotDataItem([0,0], pen=redPlotPen), pg.PlotDataItem([0,0], pen=bluePlotPen), pg.PlotDataItem([0,0], pen=greenPlotPen)
bp.addItem(bpv)
bp.addItem(bpSMB)
bp.addItem(bpvw)
bp.addItem(bpSurf)
bplV = bpLegend.addItem(bpv, 'Velocity')
bplSMB = bpLegend.addItem(bpSMB, 'SMB')
bplVW = bpLegend.addItem(bpvw, 'Velocity Width')
bplSurf = bpLegend.addItem(bpSurf, 'Surface Ele.')

nv, nsmb, nvelWidth, nsurf = None, None, None, None

l.addWidget(iiContainer, 0,0,5,1)
l.addWidget(bbW,         0,1,1,1)
l.addWidget(bp,          5,0,1,1)
l.expandingDirections()

# show the window
mw.show()
iiContainer.setMinimumHeight(2*(mw.height()//3))

# botPlotObjects = [dr, bpLegend, dataLen, bplV, bplSMB, bplVW, nv, nsmb, nvel, nsurf, botPlot]

#####################################################
####         CREATE INTERPOLATORS                ####
#####################################################
vxInterp, vyInterp = getInterpolators(vx, 'velocity', vy)
bedInterp          = getInterpolators(bed, 'bed')
smbInterp          = getInterpolators(smb, 'smb')
surfaceInterp      = getInterpolators(surfaceVals, 'surface')

#####################################################
####         SEPERATELY LOAD SURFACE             ####
#####################################################


allData = h5py.File('/home/pat/research/Data/AllDataSets.h5', 'r')
bed = allData['bed'][:]
surfaceVals = allData['surface'][:]
print 'Surface min, max: ', np.amin(surfaceVals), np.amax(surfaceVals)
vx = allData['VX'][:]
vy = allData['VY'][:]
v = sqrt(vx**2 + vy**2)
smb = allData['SMB_rec'][:]
allData.close()


#####################################################
####         SEPERATELY LOAD BED                 ####
#####################################################

allData = h5py.File('/home/pat/research/Data/AllDataSets.h5', 'r')
bed = allData['bed'][:]
surfaceVals = allData['surface'][:]
print 'Surface min, max: ', np.amin(surfaceVals), np.amax(surfaceVals)
vx = allData['VX'][:]
vy = allData['VY'][:]
v = sqrt(vx**2 + vy**2)
smb = allData['SMB_rec'][:]
allData.close()


#####################################################
####         SEPERATELY LOAD SMB                 ####
#####################################################


allData = h5py.File('/home/pat/research/Data/AllDataSets.h5', 'r')
bed = allData['bed'][:]
surfaceVals = allData['surface'][:]
print 'Surface min, max: ', np.amin(surfaceVals), np.amax(surfaceVals)
vx = allData['VX'][:]
vy = allData['VY'][:]
v = sqrt(vx**2 + vy**2)
smb = allData['SMB_rec'][:]
allData.close()




def changeMap(index):
    '''
    Called when data drop down menu is changed.
    :param index:
    :return:
    '''
    global currentMap
    if index == 0 and currentMap != 0:
        print "Show Vel"
        currentMap = 0
        iiContainer.setCurrentWidget(velW)
        iiVel.mouseClickEvent = mouseClick
    elif index == 1 and currentMap != 1:
        currentMap = 1
        print "SHow bed"
        iiContainer.setCurrentWidget(bedW)
        iiBed.mouseClickEvent = mouseClick
    elif index == 2 and currentMap != 2:
        currentMap = 2
        iiContainer.setCurrentWidget(surfaceW)


mapList.currentIndexChanged.connect(changeMap)


def cwLoop(e):
    '''
    Calls the calculate width function in order to get the width along a profile
    :param e:
    :return:
    '''
    calcVelWidth(vpts[0].getX(), vpts[0].getY(), vpts[1].getX(), vpts[1].getY(), True)

    for i in range(1, len(vpts)):
        # if not vpts[i][2]:
        calcVelWidth(vpts[i - 1].getX(), vpts[i - 1].getY(), vpts[i].getX(), vpts[i].getY(), True)

def calcVelWidth(x0, y0, x1, y1, draw):
    '''
    Calculates the width of the ice stream
    :param x0:
    :param y0:
    :param x1:
    :param y1:
    :param draw:
    :return:
    '''

    # input is in map coordinates
    #
    #    This is with interpolation
    #
    theta = np.arctan2(float(y1 - y0), float(x1-x0))
    #Rotation matrix:
    # rotMatrix = np.matrix([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    # cos    -sin    x    =    x*cos + y*-sin    = y*-sin
    # sin     cos    y    =    x*sin + y*cos    = y*cos
    tx1, ty1 = projCoord(x1,y1)
    vx0 = vxInterp(tx1,ty1, grid=False)
    vy0 = vyInterp(tx1,ty1, grid=False)
    v0 = sqrt(vx0**2 + vy0**2)
    dv = [[0, 0, 0], [0, 0, 0]]  # x, y, dv for left and right
    endPoints = [[0, 0], [0, 0]]  # end points [left[x,y], right[x,y]]
    # print 'v0 ', vx0, vy0, v0
    for i in range(2):
        dr = 0
        currentVelocity = 10
        startEndRatio = 0
        vOld = v0
        if i == 0:
            dis = 1
        else:
            dis = -1
        # print 'min([int(v0%100),8]) ', min([int(v0%100),8])
        while currentVelocity > 5 and startEndRatio <= min([int(v0%100),8]):
            dr += 1
            tx, ty = projCoord(x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta)))  # Line perpindicular to flow
            vxd = vxInterp(tx, ty, grid=False)
            vyd = vyInterp(tx, ty, grid=False)
            currentVelocity = sqrt(vxd ** 2 + vyd ** 2)
            if np.abs(currentVelocity - vOld) > dv[i][2]:
                dv[i][0], dv[i][1] = x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta))
                dv[i][2] = np.abs(currentVelocity - vOld)
            startEndRatio = v0/currentVelocity
            # print currentVelocity, startEndRatio, v0
            vOld = currentVelocity

        if currentVelocity < 5:
            # plotting line
            endPoints[i][0], endPoints[i][1] = dv[i][0], dv[i][1]#mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
        else:
            endPoints[i][0], endPoints[i][1] = mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
    # print endPoints
    # print dv
    if draw:
        iiContainer.currentWidget().addItem(pg.PlotDataItem([endPoints[0][0], endPoints[1][0]], [endPoints[0][1], endPoints[1][1]], connect='all', pen=plotPen2))

        # circle plotting
        d = (0.5)*sqrt((endPoints[0][0]-endPoints[1][0])**2 + (endPoints[0][1]-endPoints[1][1])**2)
        cax, cay = endPoints[1][0] + (d * -np.sin(theta)), endPoints[1][1] + (d * np.cos(theta))
        xc, yc = circArr(cax, cay)
        iiContainer.currentWidget().addItem(pg.PlotDataItem(xc, yc, connect='all', pen=plotPen3))
        iiContainer.currentWidget().addItem(pg.PlotDataItem([cax], [cay], pen=plotPen3))
    return endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]


def getProfile(t,y):
    return np.array([t * (-vxInterp([y[0]], [y[1]], grid=False)), t * (-vyInterp([y[0]], [y[1]], grid=False))])

integrateLine = None

def _intMesh0(t,y):
    return np.array([t * (vxInterp([y[0]], [y[1]], grid=False)), t * (vyInterp([y[0]], [y[1]], grid=False))])

def _intMesh(x,y, traceIt):
    # x0p, y0p = projCoord(x,y)
    y0 = np.array([x, y])
    t0, t1, dt = 0, 400, .1
    r = ode(_intMesh0).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    ix = 5
    iy = 5
    ox = []
    oy = []
    while r.successful() and sqrt(ix**2 + iy**2) > 0.5:
        # print(r.t+dt, r.integrate(r.t+dt))
        ai = r.integrate(r.t + dt)
        ix = vxInterp([ai[0]], [ai[1]], grid=False)
        iy = vyInterp([ai[0]], [ai[1]], grid=False)
        xi, yi = mapCoord(ai[0], ai[1])
        # print 'velocity: ', sqrt(ix**2 + iy**2)
        ox.append(np.real(xi))
        oy.append(np.real(yi))
    iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen3))
    if traceIt == 1:
        _t(ox,oy, 1)

    # r2 = ode(getProfile).set_integrator('zvode', method='bdf')
    # r2.set_initial_value(y0, t0)
    # ix = 5
    # iy = 5
    # ox = []
    # oy = []
    # while r2.successful() and sqrt(ix ** 2 + iy ** 2) > 5:
    #     # print(r.t+dt, r.integrate(r.t+dt))
    #     print 'looping'
    #     ai = r2.integrate(r2.t + dt)
    #     ix = interpX([ai[0]], [ai[1]], grid=False)
    #     iy = interpY([ai[0]], [ai[1]], grid=False)
    #     xi, yi = mapCoord(ai[0], ai[1])
    #     print 'velocity: ', sqrt(ix ** 2 + iy** 2)
    #     ox.append(np.real(xi))
    #     oy.append(np.real(yi))
    # iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen2))
    # return ox,oy



def _t(ox, oy, again):
    '''
        Draw lines perpindicular to integration line.  Repeat to trace edge.
    '''

    ipts = [] # points along edge boundry
    ds = 0.5
    mp = float(oy[-1] - oy[-2]) / float(ox[-1]-ox[-2])     # slope of user input line
    theta = np.arctan2(float(oy[-1] - oy[-2]), float(ox[-1]-ox[-2]))
    # nsampOut = 200
    # lwOut = 100
    # yline = linspace(-lwOut,lwOut, nsampOut, endpoint=True) #linspace perpindicular to the line of the 2 pts

    #Rotation matrix:
    rotMatrix = np.matrix([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    xa = [ox[-1]]
    ya = [oy[-1]]
    d = [-ds,ds]
    '''
        d = [-ds,ds]  the ds goes to the left it the line goes from down to up.
    '''

    for i in range(2):
        #rotate coordinates
        t = rotMatrix*np.matrix([[0.0], [d[i]]])
        #transform coordinates into projected coordinates
        # tx,ty = projCoord(ox[-1]+t[0,0], oy[-1]+t[1,0])
        tx, ty = ox[-1] + t[0, 0], oy[-1] + t[1, 0]
        xa.append(tx)
        ya.append(ty)
    if again%2==0:
        iiContainer.currentWidget().addItem(pg.PlotDataItem(xa, ya, pen=plotPen3))
    else:
        iiContainer.currentWidget().addItem(pg.PlotDataItem(xa, ya, pen=plotPen3))
    # _intMesh(xa[-1], ya[-1])
    ipts.append([xa[0],ya[0]])
    ipts.append([xa[1], ya[1]])
    '''
    
        Lets find all the points on the edge then integrate after that.
        
    '''

    '''
        Lets go perp then integrate a step or two and repeat
    '''


    for i in range(1):
        ds = 0.5
        mp = float(ipts[-1][1] - ipts[-2][1]) / float(ipts[-1][0] - ipts[-2][0])  # slope of user input line
        theta = np.arctan2(float(ipts[-1][1] - ipts[-2][1]), float(ipts[-1][0] - ipts[-2][0]))

        # Rotation matrix:
        rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])

        '''
            d = [-ds,ds]  the ds goes to the left it the line goes from down to up.
        '''
        arx = [ipts[-1][0]]
        ary = [ipts[-1][1]]
        iv = 0
        cnt0=0
        while iv < 250:# and cnt0 < 10:
            cnt0 += 1
            # rotate coordinates
            t = rotMatrix * np.matrix([[0.0], [-cnt0*ds]]) #-cnt0*ds means it heads back inland

            # transform coordinates into projected coordinates
            tx,ty = projCoord(ipts[-1][0] + t[0, 0], ipts[-1][1] + t[1, 0])
            ivx = vxInterp([tx], [ty], grid=False)
            ivy = vyInterp([tx], [ty], grid=False)
            iv = sqrt(ivx**2 + ivy**2)
            arx.append(ipts[-1][0] + t[0, 0])
            ary.append(ipts[-1][1] + t[1, 0])
        iiContainer.currentWidget().addItem(pg.PlotDataItem(arx, ary, pen=plotPen3))

        y0 = np.array([arx[-1],ary[-1]])
        t0, t1, dt = 0, 500, 100
        r = ode(getProfile).set_integrator('zvode', method='bdf')
        r.set_initial_value(y0, t0)
        ox = []
        oy = []
        iv = 1
        while r.successful() and r.t < t1:
            # print(r.t+dt, r.integrate(r.t+dt))
            print 'int int int'
            ai = r.integrate(r.t + dt)

            xi, yi = mapCoord(np.real(ai[0]), np.real(ai[1]))

            print xi,yi
            # print 'velocity: ', sqrt(ix**2 + iy**2)
            ox.append(xi)
            oy.append(yi)
        iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen4))


    #     del arx[:]
    #     del ary[:]
    #
    #     arx = [ipts[-1][0]]
    #     ary = [ipts[-1][1]]
    #     iv = 5
    #     cnt1 = 0
    #     while iv > .1 and cnt1 < 10:
    #         cnt1 += 1
    #         # rotate coordinates
    #         print 'second loop'
    #         t = rotMatrix * np.matrix([[0.0], [cnt1 * ds]])  # -i*ds means it heads back inland
    #
    #         # transform coordinates into projected coordinates
    #         tx, ty = projCoord(ipts[-1][0] + t[0, 0], ipts[-1][1] + t[1, 0])
    #         ivx = interpX([tx], [ty], grid=False)
    #         ivy = interpY([tx], [ty], grid=False)
    #         iv = sqrt(ivx ** 2 + ivy ** 2)
    #         arx.append(ipts[-1][0] + t[0, 0])
    #         ary.append(ipts[-1][1] + t[1, 0])
    #     iiContainer.currentWidget().addItem(pg.PlotDataItem(arx, ary, pen=plotPen4))



    # pxa, pya = projCoord(xa[-1], ya[-1])
    # y0 = np.array([pxa, pya])
    # t0, t1, dt = 0, 400, .1
    # r = ode(getProfile).set_integrator('zvode', method='bdf')
    # r.set_initial_value(y0, t0)
    # ix = 5
    # iy = 5
    # ox = []
    # oy = []
    # while r.successful() and sqrt(ix ** 2 + iy ** 2) > 5:
    #     # print(r.t+dt, r.integrate(r.t+dt))
    #     ai = r.integrate(r.t + dt)
    #     ix = interpX([ai[0]], [ai[1]], grid=False)
    #     iy = interpY([ai[0]], [ai[1]], grid=False)
    #     print 'looping ', sqrt(ix ** 2 + iy ** 2)
    #     xi, yi = mapCoord(ai[0], ai[1])
    #     # print 'velocity: ', sqrt(ix**2 + iy**2)
    #     ox.append(np.real(xi))
    #     oy.append(np.real(yi))
    # iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen3))
    #
    # if again < 5:
    #     _t(xa, ya, again + 1)






def trace(ox,oy):
    '''
        This function is supposed to trace along the edge of an outgoing velocity stream.
        This is to help calculate the ice paths in a region.

        ox, oy are the coordinates of the line going up towards the edge so that ox[-1] is the last point which should
        be very close to the edge.

    '''
    ds = 0.5
    mp = float(oy[-1] - oy[-2]) / float(ox[-1]-ox[-2])     # slope of user input line
    theta = np.arctan2(float(oy[-1] - oy[-2]), float(ox[-1]-ox[-2]))
    # nsampOut = 200
    # lwOut = 100
    # yline = linspace(-lwOut,lwOut, nsampOut, endpoint=True) #linspace perpindicular to the line of the 2 pts

    #Rotation matrix:
    rotMatrix = np.matrix([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    xa = [ox[-1]]
    ya = [oy[-1]]
    d = [-ds,ds]
    '''
    d = [-ds,ds]  the ds goes to the left it the line goes from down to up.
    '''

    for i in range(2):
        #rotate coordinates
        t = rotMatrix*np.matrix([[0.0], [d[i]]])

        #transform coordinates into projected coordinates
        # tx,ty = projCoord(ox[-1]+t[0,0], oy[-1]+t[1,0])
        tx, ty = ox[-1] + t[0, 0], oy[-1] + t[1, 0]
        xa.append(tx)
        ya.append(ty)
    iiContainer.currentWidget().addItem(pg.PlotDataItem(xa, ya, pen=plotPen3))


    '''
        Go back and forth along a line ~parallel to the integrated line. 
        Look for edge and create a new line from the edge.
        
        xa, ya are the two points in the line that is perpindicular to the integrated line.
        
    '''
    ds = 0.2
    mp = float(ya[-1] - ya[-2]) / float(xa[-1] - xa[-2])  # slope of user input line #FIXME head both directions
    # theta = np.arctan(mp)# do arctan2

    # Rotation matrix:
    rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    '''
        While loop to find the best point to integrate from
        
        THIS IS NEGATIVE DS
    '''
    found = False
    ivx = 5 # vx interpolated at new spot
    ivy = 5 # vy interpolated at new spot
    i = 0
    fx0 = [] # points along the line ~parallel to the integrated line
    fy0 = [] #
    while sqrt(ivx**2 + ivy**2) > 0.0 and i < 10:
        i += 1
        # rotate coordinates
        t = rotMatrix * np.matrix([[0.0], [-ds*i]])

        # transform coordinates into projected coordinates
        tx, ty = projCoord(xa[-1]+t[0,0], ya[-1]+t[1,0])
        ivx = vxInterp([tx], [ty], grid=False)
        ivy = vyInterp([tx], [ty], grid=False)
        if sqrt(ivx**2 + ivy**2) == 0.0 and i>1:
            found = True
        fx0.append(xa[-1]+t[0,0])
        fy0.append(ya[-1]+t[1,0])
    iiContainer.currentWidget().addItem(pg.PlotDataItem(fx0, fy0, pen=plotPen2))

    if not found:
        ivx = 5  # vx interpolated at new spot
        ivy = 5  # vy interpolated at new spot
        i = 0
        fx0 = []  # points along the line ~parallel to the integrated line
        fy0 = []  #
        while sqrt(ivx ** 2 + ivy ** 2) > 0.0 and i < 10:
            i += 1
            # rotate coordinates
            t = rotMatrix * np.matrix([[0.0], [ds * i]]) # DS IS POSITIVE!!!!!!!!!!!!!!!!!

            # transform coordinates into projected coordinates
            tx, ty = projCoord(xa[-1] + t[0, 0], ya[-1] + t[1, 0])
            ivx = vxInterp([tx], [ty], grid=False)
            ivy = vyInterp([tx], [ty], grid=False)
            if sqrt(ivx ** 2 + ivy ** 2) == 0.0 and i > 1:
                found = True
            fx0.append(xa[-1] + t[0, 0])
            fy0.append(ya[-1] + t[1, 0])

    # x0p, y0p = projCoord(x, y)


    # INT MESH
    y0 = np.array([fx0[-2], fy0[-2]])
    t0, t1, dt = 0, 400, .1
    r = ode(_intMesh0).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    ix = 5
    iy = 5
    ox = []
    oy = []
    while r.successful() and sqrt(ix ** 2 + iy ** 2) > 5:
        # print(r.t+dt, r.integrate(r.t+dt))
        ai = r.integrate(r.t + dt)
        ix = vxInterp([ai[0]], [ai[1]], grid=False)
        iy = vyInterp([ai[0]], [ai[1]], grid=False)
        xi, yi = mapCoord(ai[0], ai[1])
        # print 'velocity: ', sqrt(ix**2 + iy**2)
        ox.append(np.real(xi))
        oy.append(np.real(yi))
    iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen3))



def intMesh():
    rng = iiVel.getViewBox().viewRange()
    t1 = time.time()
    c = 0
    print 'view range: ', rng
    i0, i1 = np.amax([0,int(np.floor(rng[0][0]))]), np.amin([int(np.ceil(rng[0][1])), 10018])
    j0,j1 = np.amax([0, int(np.floor(rng[1][0]))]), np.amin([int(np.ceil(rng[1][1])), 17946])

    #this will scan horizontally until we hit a large velocity
    for i in range(i0, i1, int((i1-i0)//5)):
        for j in range(j0, j1, int((j1-j0)//5)):
            if np.abs(v[j][i]) >= 100.0:
                x,y = projCoord(i,j)
                # ox, ox = _intMesh(x,y)
                _intMesh(x, y, 1)
                # t(ox, ox)
                j = j1
                i = i1
                break
    t2 = time.time()
    print 'count: ', c
    print i0 , i1
    print j0, j1
    print 'time: ', t2-t1



def intLine(e):
    xp0 = vpts[-1].getX()
    yp0 = vpts[-1].getY()
    xp1 = vpts[-2].getX()
    yp1 = vpts[-2].getY()
    xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False) # for vpts[-1]


    theta = np.arctan2((yril - yrir),(xril - xrir))
    d = sqrt((yril - yrir)**2 + (xril - xrir)**2)
    dr = 10
    l = linspace(-d/2,d/2,dr, endpoint=True)

    rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    for i in range(len(l)):
        '''
            Moves along velocity width line. 
            Good spot for parallel programming.
        '''
        print l[i]
        trot = rotMatrix * np.matrix([[l[i]], [0.0]])
        x0p, y0p = projCoord((xp1 + trot[0,0]), (yp1 + trot[1, 0]))
        print 'x,y ', x0p, '  ----  ', y0p
        y0 = np.array([x0p, y0p])
        print y0
        t0, t1, dt = 0, 80, .1

        r = ode(getProfile).set_integrator('zvode', method='bdf')
        r.set_initial_value(y0, t0)
        print "Printing integration points"


        ox = []
        oy = []
        ivx, ivy = 5,5
        while r.successful() and sqrt(ivx**2 + ivy**2) > 5:
            ai = r.integrate(r.t + dt)
            xi, yi = mapCoord(ai[0], ai[1])
            ivx = vxInterp([ai[0]], [ai[1]], grid=False)
            ivy = vyInterp([ai[0]], [ai[1]], grid=False)
            # print 'xi, iy: ', xi, yi
            ox.append(np.real(xi))
            oy.append(np.real(yi))

        r2 = ode(_intMesh0).set_integrator('zvode', method='bdf')
        r2.set_initial_value(y0, 0)
        print "Printing integration points"
        ox2 = []
        oy2 = []
        ivx, ivy = 5, 5
        while r2.successful() and sqrt(ivx**2 + ivy**2) > 5:
            ai2 = r2.integrate(r2.t + dt)
            xi, yi = mapCoord(ai2[0], ai2[1])
            ivx = vxInterp([ai2[0]], [ai2[1]], grid=False)
            ivy = vyInterp([ai2[0]], [ai2[1]], grid=False)
            # print 'xi, iy: ', xi, yi
            ox2.append(np.real(xi))
            oy2.append(np.real(yi))

        # print 'ox, oy: ', ox, oy
        iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen2))
        iiContainer.currentWidget().addItem(pg.PlotDataItem(ox2, oy2, pen=plotPen2))


def calcProf(e):
    global integrateLine
    x0p, y0p = projCoord(vpts[-1].getX(), vpts[-1].getY())
    y0 = np.array([x0p, y0p])
    t0, t1, dt = 0, 80, .1
    r = ode(getProfile).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    print "Printing integration points"
    ox = [vpts[-1].getX()]
    oy = [vpts[-1].getY()]
    while r.successful() and r.t < t1:
        # print(r.t+dt, r.integrate(r.t+dt))
        ai = r.integrate(r.t+dt)
        xi, yi = mapCoord(ai[0], ai[1])
        # print 'xi, iy: ', xi, yi
        ox.append(np.real(xi))
        oy.append(np.real(yi))

    print 'ox, oy: ', ox, oy
    integrateLine = pg.PlotDataItem(ox,oy,pen=plotPen2)
    iiContainer.currentWidget().addItem(integrateLine)

def interpolateData(botPlotBool):
    '''
    Calculate the data for bottom plot or to run the model.

    If botPlotBool, calculate all the data.  Else, calculate just bed/surface.

    :return:
    '''
    global dr, bpLegend, dataLen, bplV, bplSMB, bplVW, bplSurf, nv, nsmb, nvelWidth, nsurf, botPlot
    botPlot = True
    velValues = []
    xValues = []
    smbValues = []
    surfValues = []
    bedValues = []
    linePoints = [0]
    vwValues = []
    graphX = []
    nv, nsmb, nvelWidth, nsurf = None, None, None, None,

    for i in range(1, len(vpts)):
        '''
        This part compares neighbor points to each other. 
        '''
        theta = np.arctan2(float(vpts[i].getY() - vpts[i - 1].getY()), float(vpts[i].getX() - vpts[i - 1].getX()))
        distance = (sqrt((vpts[i].getY() - vpts[i - 1].getY()) ** 2 + (vpts[i].getX() - vpts[i - 1].getX()) ** 2))
        remainder = distance%dr
        xline = linspace(0, distance, distance * (dr / 150), endpoint=True)  # * const makes it every 150/const meters
        linePoints.append(distance * (dr / 150) + linePoints[-1])

        # Rotation matrix:
        rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        px = []  # px, py are projected coordinates used to get values from the interpolator.  Projected meaning IRL values
        py = []
        for j in range(len(xline)):
            # rotate coordinates
            t = rotMatrix * np.matrix([[xline[j]], [0.0]])
            # FIXME probably more elegant way to do this
            # transform coordinates into projected coordinates
            tx, ty = projCoord(vpts[i - 1].getX() + t[0, 0], vpts[i - 1].getY() + t[1, 0])
            px.append(tx)
            py.append(ty)
            if len(px) > 1:
                graphX.append(graphX[-1] + sqrt((px[-1]-px[-2])**2 + (py[-1]-py[-2])**2))
            elif len(graphX) == 0:
                graphX.append(0)
            else:
                graphX.append(graphX[-1] + sqrt((px[-1]) ** 2 + (py[-1]) ** 2))

        ########################################
        ##    CALCULATE SURFACE ELEVATION     ##
        ########################################
        localSurface = surfaceInterp(px, py, grid=False)
        surfValues.append(localSurface)

        ########################################
        ##         CALCULATE BED              ##
        ########################################
        localBed = bedInterp(px, py, grid=False)
        bedValues.append(localBed)

        ########################################
        ##   COMPILE DATA                     ##
        ########################################
        nbed  = np.array(bedValues[0])
        nsurf = np.array(surfValues[0])
        for i in range(1, len(bedValues)):
            nsurf = np.append(nsurf, surfValues[i])
            nbed  = np.append(nbed,  bedValues[i])

        if botPlotBool:
            ########################################
            ##        CALCULATE VELOCITY          ##
            ########################################
            vxd = vxInterp(px, py, grid=False)
            vyd = vyInterp(px, py, grid=False)  # 1D array
            vi = sqrt(vxd ** 2 + vyd ** 2)
            xValues.append(xline)
            velValues.append(vi)
            print 'velocity len ', len(vi)
            print 'px len ', len(px)
            ########################################
            ##     CALCULATE VELOCITY WIDTH       ##
            ########################################

            vwd = []
            for i in range(len(px)):
                xp0, yp0 = mapCoord(px[i - 1], py[i - 1])
                xp1, yp1 = mapCoord(px[i], py[i])
                xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False)
                vwd.append(sqrt((xril - xrir) ** 2 + (yril - yrir) ** 2))
            vwValues.append(vwd)

            ########################################
            ##   CALCULATE SURFACE MASS-BALANCE   ##
            ########################################
            localSMB = smbInterp(px, py, grid=False)
            smbValues.append(localSMB)

            ########################################
            ##   COMPILE DATA                     ##
            ########################################
            nv        = np.array(velValues[0])
            nsmb      = np.array(smbValues[0])
            nvelWidth = np.array(vwValues[0])

            for i in range(1, len(velValues)):
                nv        = np.append(nv,        velValues[i])
                nsmb      = np.append(nsmb,      smbValues[i])
                nvelWidth = np.append(nvelWidth, vwValues[i])

    if botPlotBool:
        return nbed, nsurf, nv, nsmb, nvelWidth, linePoints, np.array(graphX)
    else:
        return nbed, nsurf  # , linePoints


def calcBP():
    '''
    Calculate the data for bottom plot then populate the plot.
    :return:
    '''
    #FIXME why have nv, nsmb, etc here?
    global dr, bpLegend, dataLen, bplV, bplSMB, bplSurf, bplVW, nv, nsmb, nvelWidth, nsurf, botPlot
    if len(vpts) > 0:
        print 'plotting'
        nbed, nsurf, nv, nsmb, nvelWidth, linePoints, graphX = interpolateData(True)
        bpv.setData(graphX   , nv) #FIXME shouldnt have to use :-2
        bpSMB.setData(graphX , nsmb)
        bpvw.setData(graphX  , nvelWidth) #velocity width
        bpSurf.setData(graphX, nsurf)
        pg.QtGui.QApplication.processEvents()

def arrows():
    ''''
    Plot velocity vectors
    '''
    print 'arrows'
    print iiVel.getViewBox().viewRange() # [y0, y1] [x0, x1]
    rngx, rngy = iiVel.getViewBox().viewRange()
    print rngx
    print rngy
    print 'Starting arrows'
    for i in range(int(rngx[0]),int(rngx[1])):
        for j in range(int(rngy[0]), int(rngy[1])):
            vdir = [vx[j][i], -vy[j][i]]
            vmag = sqrt(vx[j][i]**2 +  vy[j][i]**2)
            # theta = np.arctan2(vdir[1], vdir[0])
            iiContainer.currentWidget().addItem(pg.PlotDataItem([(i+0.5),(i+0.5 + vdir[0]/(1.5*vmag))], [(j+0.5),(j+0.5 + vdir[1]/(1.5*vmag))], pen=plotPen2))
            iiContainer.currentWidget().addItem(pg.PlotDataItem([i+0.5], [j+0.5]), pen=(255,255,255), symbolBrush=(255,0,0), symbolPen='w')
    print 'finished arrows'



def mouseMovedBP(evt):
    global bpLegend, dataLen, nv, nsmb, nvelWidth, bplV, bplSMB, bplV, bplSurf
    print 'mouse moved bottom plot'
    print 'bplV is ' , bplV is None
    if botPlot:
        print 'is not not not none'
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if bp.getPlotItem().sceneBoundingRect().contains(pos):
            mousePoint = bp.getPlotItem().vb.mapSceneToView(pos)
            index = int(mousePoint.x())
            if index > 0 and index < dataLen:
                bplV.setText('Velocity' + str(nv[index]))
                bplSMB.setText('SMB ' + str(nsmb[index]))
                bplVW.setText('Velocity Width ' + str(nvelWidth[index]))
                bplSurf.setTxt('Surface ele ' + str(nsurf[index]))
    # else:
        # print 'It is NONE'

def clearPoints():
    del vpts[:]


def runModelButt():
    global dr
    if len(vpts) > 0:
        nbed, nsurf = interpolateData(False)
        print 'bed data:'
        print nbed
        print 'surface data:'
        print nsurf
        THICKLIMIT = 10.  # Ice is never less than this thick
        H = nsurf - nbed
        nsurf[H <= THICKLIMIT] = nbed[H <= THICKLIMIT]
        N = len(nbed)
        mesh = fc.IntervalMesh(N-1, 0, dr*(N-1))
        #FIXME the intervalMesh is consistantly 150 between each datapoint this not true for the data being sent
        hdf_name = '/home/pat/research/latest_profile.h5'
        hfile = fc.HDF5File(mesh.mpi_comm(), hdf_name, "w")
        V = fc.FunctionSpace(mesh,"CG",1)
        functBed     = fc.Function(V, name="Bed")
        functSurface = fc.Function(V, name="Surface")
        functBed.vector()[:]     = nbed
        functSurface.vector()[:] = nsurf
        hfile.write(functBed.vector(), "/bed")
        hfile.write(functSurface.vector(), "/surface")
        hfile.write(mesh, "/mesh")
        hfile.close()
        runModel(hdf_name)



# proxy = pg.SignalProxy(bp.getPlotItem().scene().sigMouseMoved, rateLimit=60, slot=mouseMovedBP)

vptSel = False
vptCur = None
clearButton.clicked.connect(clearPoints)
calcWidthButton.clicked.connect(cwLoop)
intButton.clicked.connect(calcProf)
cProfButton.clicked.connect(calcBP) #FIXME should change names so calcProf isn't the integration function
cRegionButton.clicked.connect(intLine)
cVelArrowsButton.clicked.connect(arrows)
modelButton.clicked.connect(runModelButt)



def mouseClick(e):
    global vptSel, vptCur, integrateLine
    first = False
    if len(vpts) == 0:
        first = True
        x = e.pos().x()
        y = e.pos().y()
        print 'map coords ', x,y
        print 'proj coords ', projCoord(x,y)
        x0p, y0p = projCoord(x,y)
        print 'x,y ', x0p, '  ----  ', y0p
        y0 = np.array([x0p, y0p])
        print y0
        t0, t1, dt = 0, 80, .1
        r = ode(getProfile).set_integrator('zvode', method='bdf')
        r.set_initial_value(y0, t0)
        ox = []
        oy = []
        ind = 0
        while r.successful() and ind < 2:
            ai = r.integrate(r.t + dt)
            ind += 1
            xi, yi = mapCoord(np.real(ai[0]), np.real(ai[1]))
            ox.append(xi)
            oy.append(yi)
        iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen3))

        endPoints = [[0,0], [0,0]]
        endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1] = calcVelWidth(x, y, ox[-1], oy[-1], False)
        d = (0.5) * sqrt((endPoints[0][0] - endPoints[1][0]) ** 2 + (endPoints[0][1] - endPoints[1][1]) ** 2)
        theta = np.arctan2(float(y - oy[-1]), float(x - ox[-1]))
        x, y = endPoints[0][0] + (d * -np.sin(theta)), endPoints[0][1] + (d * np.cos(theta))
        print 'x,y ', x,y


    if shift:
        intLine(e)
    else:
        if not vptSel:
            for pt in vpts:
                if pt.checkClicked(e.pos()):
                    vptSel = True
                    vptCur = pt
                    print 'Found it! at ', pt.getPos()
            if not vptSel:
                if not first:
                    x = e.pos().x()
                    y = e.pos().y()
                px, py = projCoord(x,y)
                vxd = vxInterp([px], [py], grid=False)
                vyd = vyInterp([px], [py], grid=False)
                sur = surfaceInterp([px], [py], grid=False)
                velocity = sqrt(vxd**2 + vyd**2)
                vpts.append(vpt(x, y, velocity, velW)) #in projected coordinates, like 600000
                txt = 'x: ' + str(x) + '\t y: ' + str(y) + '\t v: ' + str(velocity)+ '\t surface ele: ' + str(sur)
                velTF.append(txt)
                iiContainer.currentWidget().addItem(vpts[-1].getCross()[0])
                iiContainer.currentWidget().addItem(vpts[-1].getCross()[1])
                if len(vpts) > 1:
                    xa = [vpts[-1].getX(), vpts[-2].getX()]
                    ya = [vpts[-1].getY(), vpts[-2].getY()]
                    vpts[-2].setLine(pg.PlotDataItem(xa,ya,connect='all'))
                    iiContainer.currentWidget().addItem(vpts[-2].getLine())#,pen=plotPen)
        else:
            vptSel = False


mouseCoordinates = QtGui.QLabel('x:\ty:')
buttonBox.addWidget(mouseCoordinates)

def mouseMoved(e):
    global vptSel, vptCur, integrateLine, currentMap
    # print 'mm 1'
    if e.isExit() is False:
        # print 'mm 2'
        if vptSel and integrateLine is not None:
            cData = integrateLine.curve.getData()
            imin = curveDistance(e.pos().x(), e.pos().y(), cData)
            if imin != -1:
                x = cData[0][imin]
                y = cData[1][imin]
                vptCur.updateCross(x,y)
            else:
                vptCur.updateCross(e.pos().x(), e.pos().y())
        # else:
        x = int(np.floor(e.pos().x()))
        y = int(np.floor(e.pos().y()))
        if np.abs(x) <= 10018 and np.abs(y) <= 17946:
            # print 'mm 3'
            if currentMap == 0:
                # print 'mm 30'
                mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y) + '\t v: ' + str(v[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))])+ '\t bed: ' + str(bed[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))])  + '\t surface: ' + str(surfaceVals[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))]))
            elif currentMap == 1:
                # print 'mm 31'
                mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y) + '\t bed: ' + str(bed[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))]))
            elif currentMap == 2:
                # print 'mm 32'
                mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y) + '\t surface: ' + str(surfaceVals[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))]))



def ky(e):
    print 'key pressed ', e.key()
    print 'type: -' + str(e.type()) + '-'
    #6 is press, 7 is release
    if e.type() == 6:
        print 'shift true'
        shift = True
    else:
        shift = False
        print 'shift false'
    # 16777249 is shift

shift = False
iiVel.hoverEvent = mouseMoved
iiBed.hoverEvent = mouseMoved
iiSurface.hoverEvent = mouseMoved
# proxy = pg.SignalProxy(bp.getPlotItem().scene().sigMouseMoved, rateLimit=60, slot=mouseMovedBP)

iiVel.mouseClickEvent = mouseClick # Default map is vel FIXME

# iiContainer.keyboardGrabber()
# iiContainer.keyPressEvent = ky
# iiContainer.keyReleaseEvent = ky


print 'Done Loading'

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()



'''

Monday June 12: 
Added constant cursor coordinates with velocity readout

*Jesse's width search is strange.  Shouldn't look over such a large area.

Tuesday June 13:
Did velocity width WITHOUT interpolation.  Interpolation makes the dv so small that the program 
has troubles finding large drops.  Could be improved by maybe makine 3 lines instead of 1


Wednesday June 14:
used np.flipup on the 

'''