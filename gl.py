import time
startTime = time.time()
import fenics as fc
import sys

from pyqtgraph.Qt import QtCore  # , QtWidgets
from scipy.integrate import ode
from helper_files.classes.dataset import dataset
from helper_files.classes.pltPoints import *
from helper_files.cm import *
from helper_files.gui import *
from helper_files.math_functions import *
from helper_files.pens import *
from helper_files.profile_driver import runModel
from helper_files.dataset_objects import *
from helper_files.data_functions import *
from helper_files.gui_functions import *

# from helper_files.gui_functions import calcProf
sys.path.append("/home/pat/")
'''
https://stackoverflow.com/questions/38065570/pyqtgraph-is-it-possible-to-have-a-imageview-without-histogram

RACMO:

Latitude of the origin:     90
Longitude of the origin(central meridian):         -45
Standard parallel:          70
Ellipsoid:WGS84
Datum:WGS84
'''

print 'Loading...'
#####################################################
####           CONSTANTS/GLOBAL VARIABLES       ####
#####################################################

# dr = 150         # spatial resolution for bottom plot interpolation in meters
botPlot = False  # if bottom plot has been populated or not
inBed = False
inSMB = False
inSurface = False
clickedCurve = False
dataCMFileName = './data/dataCMValues.h5'
dataFileName   = './data/AllDataSets.h5'

#
# #####################################################
# ####         CREATE INITIAL DATA SET(S)          ####
# #####################################################
#
# velocity = dataset('velocity', bpLegend, greyPlotPen, map=True)
# bp.addItem(velocity.pathPlotItem)
# iiContainer.addWidget(velocity.plotWidget)
# iiContainer.setCurrentWidget(velocity.plotWidget)
#
# smb = dataset('smb', bpLegend, redPlotPen, map=True)
# print 'SMB: ', np.amin(smb.data), np.amax(smb.data) #-11493.3860928 6060.80339304
# bp.addItem(smb.pathPlotItem)
# iiContainer.addWidget(smb.plotWidget)
#
# surface = dataset('surface', bpLegend, greenPlotPen, map=True)
# bp.addItem(surface.pathPlotItem)
# iiContainer.addWidget(surface.plotWidget)
#
# bed = dataset('bed', bpLegend, bluePlotPen, map=True)
# bp.addItem(bed.pathPlotItem)
# iiContainer.addWidget(bed.plotWidget)
#
# velocityWidth = dataset('velocitywidth', bpLegend, purplePlotPen)

#
# def changeMap(index):
#     '''
#     Called when data-set drop down menu is changed.
#     :param index:
#     :return:
#     '''
#
#     # print iiContainer.currentWidget()
#     vr = iiContainer.currentWidget().getPlotItem().getViewBox().viewRange()
#     global currentMap
#     if index == 0 and currentMap != 0:
#         #velocity
#         currentMap = 0
#         iiContainer.setCurrentWidget(velocity.plotWidget)
#         velocity.imageItem.mouseClickEvent = mouseClick
#         velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])
#     elif index == 1 and currentMap != 1:
#         #bed
#         currentMap = 1
#         iiContainer.setCurrentWidget(bed.plotWidget)
#         bed.imageItem.mouseClickEvent = mouseClick
#         bed.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])
#     elif index == 2 and currentMap != 2:
#         #surface
#         currentMap = 2
#         iiContainer.setCurrentWidget(surface.plotWidget)
#         surface.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])
#     elif index == 3 and currentMap != 3:
#         #SMB
#         currentMap = 3
#         iiContainer.setCurrentWidget(smb.plotWidget)
#         smb.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])


mapList.currentIndexChanged.connect(changeMap)


# def cwLoop(e):
#     '''
#     Calls the calculate width function in order to get the width along a profile
#     :param e:
#     :return:
#     '''
#     calcVelWidth(vpts[0].getX(), vpts[0].getY(), vpts[1].getX(), vpts[1].getY(), True)
#
#     for i in range(1, len(vpts)):
#         # if not vpts[i][2]:
#         calcVelWidth(vpts[i - 1].getX(), vpts[i - 1].getY(), vpts[i].getX(), vpts[i].getY(), True)
#
# def calcVelWidth(x0, y0, x1, y1, draw):
#     '''
#     Calculates the width of the ice stream at one point, (x1, y1).  (x0, y0) is there
#     to give an idea of where the velocity width begins and ends which should be on a
#     line which is perpindicular to the line from (x1, y1) to (x0, y0).
#     :param x0:
#     :param y0:
#     :param x1:
#     :param y1:
#     :param draw:
#     :return:
#     '''
#
#     # input is in map coordinates
#     #
#     #    This is with interpolation
#     #
#     theta = np.arctan2(float(y1 - y0), float(x1-x0))
#     #Rotation matrix:
#     # rotMatrix = np.matrix([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
#     # cos    -sin    x    =    x*cos + y*-sin    = y*-sin
#     # sin     cos    y    =    x*sin + y*cos    = y*cos
#     vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', x1, y1, d2=velocity.vy)
#     tx1, ty1 = projCoord(x1,y1)
#
#     vx0 = vxInterp(tx1,ty1, grid=False)
#     vy0 = vyInterp(tx1,ty1, grid=False)
#     v0 = sqrt(vx0**2 + vy0**2)
#
#     dv = [[0, 0, 0], [0, 0, 0]]  # x, y, dv for left and right
#     endPoints = [[0, 0], [0, 0]]  # end points [left[x,y], right[x,y]]
#     # print 'v0 ', vx0, vy0, v0
#     for i in range(2):
#         dr = 0
#         currentVelocity = 10
#         startEndRatio = 0
#         vOld = v0
#         if i == 0:
#             dis = 1
#         else:
#             dis = -1
#         # print 'min([int(v0%100),8]) ', min([int(v0%100),8])
#         while currentVelocity > 5 and startEndRatio <= min([int(v0%100),8]):
#             dr += 1
#             vxInterp, vyInterp = getInterpolators(velocity.vx, velocity.name, (x1 + (dr*dis * -np.sin(theta))), (y1 + (dr*dis * np.cos(theta))), d2=velocity.vy)
#             tx, ty = projCoord(x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta)))  # Line perpindicular to flow
#             vxd = vxInterp(tx, ty, grid=False)
#             vyd = vyInterp(tx, ty, grid=False)
#             currentVelocity = sqrt(vxd ** 2 + vyd ** 2)
#             if np.abs(currentVelocity - vOld) > dv[i][2]:
#                 dv[i][0], dv[i][1] = x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta))
#                 dv[i][2] = np.abs(currentVelocity - vOld)
#             startEndRatio = v0/currentVelocity
#             # print currentVelocity, startEndRatio, v0
#             vOld = currentVelocity
#         if currentVelocity < 5:
#             # plotting line
#             endPoints[i][0], endPoints[i][1] = dv[i][0], dv[i][1]#mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
#         else:
#             endPoints[i][0], endPoints[i][1] = mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
#     if draw:
#         iiContainer.currentWidget().addItem(pg.PlotDataItem([endPoints[0][0], endPoints[1][0]], [endPoints[0][1], endPoints[1][1]], connect='all', pen=plotPen2))
#
#         # circle plotting
#         d = (0.5)*sqrt((endPoints[0][0]-endPoints[1][0])**2 + (endPoints[0][1]-endPoints[1][1])**2)
#         cax, cay = endPoints[1][0] + (d * -np.sin(theta)), endPoints[1][1] + (d * np.cos(theta))
#         xc, yc = circArr(cax, cay)
#         iiContainer.currentWidget().addItem(pg.PlotDataItem(xc, yc, connect='all', pen=plotPen3))
#         iiContainer.currentWidget().addItem(pg.PlotDataItem([cax], [cay], pen=plotPen3))
#     return endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]


# def getProfile(t,y):
#     '''
#     Prints a line
#     :param t:
#     :param y:
#     :return:
#     '''
#     # print 'getPro ', np.real(y[0]), ' ', np.real(y[1])
#     mx, my = mapCoord(np.real(y[0]), np.real(y[1]))
#     vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', np.floor(mx), np.floor(my), d2=velocity.vy)
#     return np.array([t * (-vxInterp([y[0]], [y[1]], grid=False)), t * (-vyInterp([y[0]], [y[1]], grid=False))])

integrateLine = None

#
# def intLine(e):
#     '''
#     Shows the line which follows the path in the velocity field
#     :param e:
#     :return:
#     '''
#     xp0 = vpts[-1].getX()
#     yp0 = vpts[-1].getY()
#     xp1 = vpts[-2].getX()
#     yp1 = vpts[-2].getY()
#     xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False) # for vpts[-1]
#
#     theta = np.arctan2((yril - yrir),(xril - xrir))
#     d = sqrt((yril - yrir)**2 + (xril - xrir)**2)
#     dr = 10
#     l = linspace(-d/2,d/2,dr, endpoint=True)
#
#     rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])
#
#     for i in range(len(l)):
#         '''
#             Moves along velocity width line.
#             Good spot for parallel programming.
#         '''
#         print l[i]
#         trot = rotMatrix * np.matrix([[l[i]], [0.0]])
#         x0p, y0p = projCoord((xp1 + trot[0,0]), (yp1 + trot[1, 0]))
#         print 'x,y ', x0p, '  ----  ', y0p
#         y0 = np.array([x0p, y0p])
#         print y0
#         t0, t1, dt = 0, 80, .1
#
#         r = ode(getProfile).set_integrator('zvode', method='bdf')
#         r.set_initial_value(y0, t0)
#         print "Printing integration points"
#
#         ox = []
#         oy = []
#         ivx, ivy = 5,5
#         while r.successful() and sqrt(ivx**2 + ivy**2) > 5:
#             ai = r.integrate(r.t + dt)
#             xi, yi = mapCoord(ai[0], ai[1])
#             ivx = vxInterp([ai[0]], [ai[1]], grid=False)
#             ivy = vyInterp([ai[0]], [ai[1]], grid=False)
#             # print 'xi, iy: ', xi, yi
#             ox.append(np.real(xi))
#             oy.append(np.real(yi))
#
#         r2 = ode(_intMesh0).set_integrator('zvode', method='bdf')
#         r2.set_initial_value(y0, 0)
#         print "Printing integration points"
#         ox2 = []
#         oy2 = []
#         ivx, ivy = 5, 5
#         while r2.successful() and sqrt(ivx**2 + ivy**2) > 5:
#             ai2 = r2.integrate(r2.t + dt)
#             xi, yi = mapCoord(ai2[0], ai2[1])
#             ivx = vxInterp([ai2[0]], [ai2[1]], grid=False)
#             ivy = vyInterp([ai2[0]], [ai2[1]], grid=False)
#             # print 'xi, iy: ', xi, yi
#             ox2.append(np.real(xi))
#             oy2.append(np.real(yi))
#
#         # print 'ox, oy: ', ox, oy
#         iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen2))
#         iiContainer.currentWidget().addItem(pg.PlotDataItem(ox2, oy2, pen=plotPen2))


# def calcProf(e):
#     '''
#     Calculates the line that shows velocity flow inwards (negative direction).
#     :param e:
#     :return:
#     '''
#     global integrateLine
#     x0p, y0p = projCoord(vpts[-1].getX(), vpts[-1].getY())
#     y0 = np.array([x0p, y0p])
#     t0, t1, dt = 0, 80, .1
#     r = ode(getProfile).set_integrator('zvode', method='bdf')
#     r.set_initial_value(y0, t0)
#     print "Printing integration points"
#     ox = [vpts[-1].getX()]
#     oy = [vpts[-1].getY()]
#     while r.successful() and r.t < t1:
#         # print(r.t+dt, r.integrate(r.t+dt))
#         ai = r.integrate(r.t+dt)
#         xi, yi = mapCoord(ai[0], ai[1])
#         # print 'xi, iy: ', xi, yi
#         ox.append(np.real(xi))
#         oy.append(np.real(yi))
#
#     print 'ox, oy: ', ox, oy
#     integrateLine = pg.PlotDataItem(ox,oy,pen=plotPen2)
#     iiContainer.currentWidget().addItem(integrateLine)

# def interpolateData(botPlotBool):
#     '''
#     Calculate the data for bottom plot or to run the model.
#
#     If botPlotBool, calculate all the data.  Else, calculate just bed/surface.
#
#     :return:
#     '''
#     global dr, bpLegend, dataLen, botPlot
#     botPlot = True
#     velValues = []
#     xValues = []
#     smbValues = []
#     surfValues = []
#     bedValues = []
#     linePoints = [0]
#     vwValues = []
#     graphX = []
#
#     ########################################
#     ##    GATHER LOCAL INTERPOLATORS      ##
#     ########################################
#     mxx = max(pt.x for pt in vpts)
#     mxy = max(pt.y for pt in vpts)
#     mix = min(pt.x for pt in vpts)
#     miy = min(pt.y for pt in vpts)
#
#     surfaceInterp      = getInterpolators(surface.data, surface.name,   mix, miy, x1=mxx, y1=mxy)
#     bedInterp          = getInterpolators(bed.data,     bed.name,       mix, miy, x1=mxx, y1=mxy)
#     vxInterp, vyInterp = getInterpolators(velocity.vx,  velocity.name,  mix, miy, x1=mxx, y1=mxy, d2=velocity.vy)
#     smbInterp          = getInterpolators(smb.data,     smb.name,       mix, miy, x1=mxx, y1=mxy)
#
#     # Find a distance ~150m which gets close to dividing the distance between first 2 spots
#
#     for i in range(1, len(vpts)):
#         '''
#         This part compares neighbor points to each other.
#         '''
#         theta = np.arctan2(float(vpts[i].getY() - vpts[i - 1].getY()), float(vpts[i].getX() - vpts[i - 1].getX()))
#         # pvx0, pvy0 = projCoord(vpts[i-1].x, vpts[i-1].y)
#         # pvx1, pvy1 = projCoord(vpts[i].x, vpts[i].y)
#         # distance = (sqrt((vpts[i].getY() - vpts[i - 1].getY()) ** 2 + (vpts[i].getX() - vpts[i - 1].getX()) ** 2))
#         distance = sqrt((vpts[i-1].x-vpts[i].x)**2 + (vpts[i-1].y-vpts[i].y)**2)
#         remainder = distance%dr
#         # Xline needs to be in map coordinates because tx,ty are in map coordinates
#         xline = linspace(0, distance, distance, endpoint=True)  # * const makes it every 150/const meters
#         '''
#         #FIXME NEED TO CHANGE SO IT LINES UP WITH FENICS MESH
#         '''
#
#         yline = linspace(0, distance, distance/150, endpoint=False)  # * const makes it every 150/const meters
#         # print 'xline: '
#         # print xline
#         # print 'yline: '
#         # print yline
#         linePoints.append(distance * (1 / dr) + linePoints[-1])
#
#         # Rotation matrix:
#         rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])
#         px = []  # px, py are projected coordinates used to get values from the interpolator.  Projected meaning IRL values
#         py = []
#         for j in range(len(xline)):
#             # rotate coordinates
#             t = rotMatrix * np.matrix([[xline[j]], [0.0]])
#             # FIXME probably more elegant way to do this
#             # transform coordinates into projected coordinates
#             tx, ty = projCoord(vpts[i - 1].getX() + t[0, 0], vpts[i - 1].getY() + t[1, 0])
#             px.append(tx)
#             py.append(ty)
#             if len(px) > 1:
#                 graphX.append(graphX[len(graphX)-1] + sqrt((px[-1]-px[-2])**2 + (py[-1]-py[-2])**2))
#                 # print 'dist: ', graphX[-1] + sqrt((px[-1]-px[-2])**2 + (py[-1]-py[-2])**2)
#             elif len(graphX) == 0:
#                 graphX.append(0)
#             else:
#                 graphX.append(graphX[len(graphX)-1])# + sqrt((px[-1]) ** 2 + (py[-1]) ** 2))
#             #     print 'dis3: ', graphX[-1] + sqrt((px[-1]) ** 2 + (py[-1]) ** 2)
#         print 'px ', px
#         print 'py ', py
#         print 'gx ', graphX
#
#         ########################################
#         ##    CALCULATE SURFACE ELEVATION     ##
#         ########################################
#         localSurface = surfaceInterp(px, py, grid=False)
#         surfValues.append(localSurface)
#
#         ########################################
#         ##         CALCULATE BED              ##
#         ########################################
#         localBed = bedInterp(px, py, grid=False)
#         bedValues.append(localBed)
#
#         ########################################
#         ##   COMPILE DATA                     ##
#         ########################################
#         bed.pathData     = np.array(bedValues[0])
#         surface.pathData = np.array(surfValues[0])
#         for i in range(1, len(bedValues)):
#             bed.pathData     = np.append(bed.pathData, bedValues[i])
#             surface.pathData = np.append(surface.pathData, surfValues[i])
#
#
#         if botPlotBool:
#             ########################################
#             ##        CALCULATE VELOCITY          ##
#             ########################################
#             vxd = vxInterp(px, py, grid=False)
#             vyd = vyInterp(px, py, grid=False)  # 1D array
#             vi = sqrt(vxd ** 2 + vyd ** 2)
#             xValues.append(xline)
#             velValues.append(vi)
#             print 'velocity len ', len(vi)
#             print 'px len ', len(px)
#             ########################################
#             ##     CALCULATE VELOCITY WIDTH       ##
#             ########################################
#
#             vwd = []
#             for i in range(len(px)):
#                 xp0, yp0 = mapCoord(px[i - 1], py[i - 1])
#                 xp1, yp1 = mapCoord(px[i], py[i])
#                 xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False)
#                 vwd.append(sqrt((xril - xrir) ** 2 + (yril - yrir) ** 2))
#             vwValues.append(vwd)
#
#             ########################################
#             ##   CALCULATE SURFACE MASS-BALANCE   ##
#             ########################################
#             localSMB = smbInterp(px, py, grid=False)
#             smbValues.append(localSMB)
#
#             ########################################
#             ##   COMPILE DATA                     ##
#             ########################################
#             velocity.pathData      = np.array(velValues[0])
#             smb.pathData           = np.array(smbValues[0])
#             velocityWidth.pathData = np.array(vwValues[0])
#
#             for i in range(1, len(velValues)):
#                 velocity.pathData      = np.append(velocity.pathData,        velValues[i])
#                 smb.pathData           = np.append(smb.pathData,      smbValues[i])
#                 velocityWidth.pathData = np.append(velocityWidth.pathData, vwValues[i])
#
#
#     print 'graphx[-1]: ', graphX[len(graphX)-1]
#     dist = 0
#     for i in range(len(vpts)-1):
#         print 'calc distance...'
#         xd0, yd0 = projCoord(vpts[i].x, vpts[i].y)
#         xd1, yd1 = projCoord(vpts[i+1].x, vpts[i+1].y)
#         print xd0, yd0, xd1, yd1
#         dist += sqrt(((xd1 - xd0) ** 2 + (yd1 - yd0) ** 2))
#     print 'dist: ', dist
#     surface.distanceData = graphX
#     bed.distanceData = graphX
#     if botPlotBool:
#         velocity.distanceData = graphX
#         smb.distanceData = graphX
#         velocityWidth.distanceData = graphX
#         return linePoints, np.array(graphX)
    # else:
    #     return nbed, nsurf  # , linePoints


# def calcBP():
#     '''
#     Calculate the data for bottom plot then populate the plot.
#     :return:
#     '''
#     global dr, bpLegend, dataLen, botPlot
#     if len(vpts) > 0:
#         print 'plotting'
#         # nbed, nsurf, nv, nsmb, nvelWidth, linePoints, graphX = interpolateData(True)
#         linePoints, graphX = interpolateData(True)
#         # velocity.pathPlotItem.setData(velocity.distanceData          , velocity.pathData)
#         # smb.pathPlotItem.setData(smb.distanceData                    , smb.pathData)
#         # velocityWidth.pathPlotItem.setData(velocityWidth.distanceData, velocityWidth.pathData)
#         surface.pathPlotItem.setData(surface.distanceData              , surface.pathData)
#         # bed.pathPlotItem.setData(bed.distanceData                    , bed.pathData)
#         pg.QtGui.QApplication.processEvents()


# def arrows():
#     ''''
#     Plot velocity vectors to see the velocity field and its borders
#     '''
#     print 'arrows'
#     print velocity.imageItem.getViewBox().viewRange() # [y0, y1] [x0, x1]
#     rngx, rngy = velocity.imageItem.getViewBox().viewRange()
#     print rngx
#     print rngy
#     print 'Starting arrows'
#     for i in range(int(rngx[0]),int(rngx[1])):
#         for j in range(int(rngy[0]), int(rngy[1])):
#             vdir = [velocity.vx[j][i], -velocity.vy[j][i]]
#             vmag = sqrt(velocity.vx[j][i]**2 +  velocity.vy[j][i]**2)
#             # theta = np.arctan2(vdir[1], vdir[0])
#             iiContainer.currentWidget().addItem(pg.PlotDataItem([(i+0.5),(i+0.5 + vdir[0]/(1.5*vmag))], [(j+0.5),(j+0.5 + vdir[1]/(1.5*vmag))], pen=plotPen2))
#             iiContainer.currentWidget().addItem(pg.PlotDataItem([i+0.5], [j+0.5]), pen=(255,255,255), symbolBrush=(255,0,0), symbolPen='w')
#     print 'finished arrows'



# def mouseMovedBP(evt):
#     global bpLegend, dataLen
#     print 'mouse moved bottom plot'
#     if botPlot:
#         print 'is not not not none'
#         pos = evt[0]  ## using signal proxy turns original arguments into a tuple
#         if bp.getPlotItem().sceneBoundingRect().contains(pos):
#             mousePoint = bp.getPlotItem().vb.mapSceneToView(pos)
#             index = int(mousePoint.x())
#             if index > 0 and index < len(velocity.pathData):
#                 velocity.legendItem.setText('Velocity' + str(velocity.pathData[index]))
#                 smb.legendItem.setText('SMB ' + str(smb.pathData[index]))
#                 velocityWidth.legendItem.setText('Velocity Width ' + str(velocityWidth.pathData[index]))
#                 surface.legendItem.setTxt('Surface ele ' + str(surface.pathData[index]))
#     # else:
#         # print 'It is NONE'

def clearPoints():
    del vpts[:]
    velocity.pathPlotItem.clear()
    surface.pathPlotItem.clear()
    smb.pathPlotItem.clear()
    bed.pathPlotItem.clear()


def runModelButt():
    global dr
    if len(vpts) > 0:

        interpolateData(True)
        # print 'bed data:'
        # print bed.pathData
        # print 'surface data:'
        # print surface.pathData

        THICKLIMIT = 10.  # Ice is never less than this thick
        H = surface.pathData - bed.pathData
        surface.pathData[H <= THICKLIMIT] = bed.pathData[H <= THICKLIMIT]
        N = len(bed.pathData)

        mesh = fc.IntervalMesh(N-1, 0, dr*(N-1))

        #FIXME the intervalMesh is consistantly 150 between each datapoint this not true for the data being sent
        hdf_name = '/home/pat/research/latest_profile.h5'
        hfile = fc.HDF5File(mesh.mpi_comm(), hdf_name, "w")
        V = fc.FunctionSpace(mesh,"CG",1)

        functBed      = fc.Function(V, name="Bed")
        functSurface  = fc.Function(V, name="Surface")
        functSMB      = fc.Function(V, name='SMB')
        functVelocity = fc.Function(V, name='Velocity')

        functBed.vector()[:]      = bed.pathData
        functSurface.vector()[:]  = surface.pathData
        functSMB.vector()[:]      = smb.pathData
        functVelocity.vector()[:] = velocity.pathData

        hfile.write(functBed.vector(), "/bed")
        hfile.write(functSurface.vector(), "/surface")
        hfile.write(functSMB.vector(), '/smb')
        hfile.write(functVelocity.vector(), '/velocity')
        hfile.write(mesh, "/mesh")
        hfile.close()
        runModel(hdf_name)



# proxy = pg.SignalProxy(bp.getPlotItem().scene().sigMouseMoved, rateLimit=60, slot=mouseMovedBP)


vptCur = None
clearButton.clicked.connect(clearPoints)
calcWidthButton.clicked.connect(cwLoop)
intButton.clicked.connect(calcProf)
cProfButton.clicked.connect(calcBP) #FIXME should change names so calcProf isn't the integration function
cRegionButton.clicked.connect(intLine)
cVelArrowsButton.clicked.connect(arrows)
modelButton.clicked.connect(runModelButt)


# def mouseClick(e):
#     global vptSel, vptCur, integrateLine
#     first = False
#     print 'mouseclick ', e.pos().x(), ' ', e.pos().y()
#     if len(vpts) == 0:
#         first = True
#         x = e.pos().x()
#         y = e.pos().y()
#         x0p, y0p = projCoord(x,y)
#         y0 = np.array([x0p, y0p])
#         t0, t1, dt = 0, 80, .1
#         r = ode(getProfile).set_integrator('zvode', method='bdf')
#         r.set_initial_value(y0, t0)
#         ox = []
#         oy = []
#         ind = 0
#         while r.successful() and ind < 2:
#             '''
#             what the fuck was I trying to do here??
#             '''
#             ai = r.integrate(r.t + dt)
#             ind += 1
#             xi, yi = mapCoord(np.real(ai[0]), np.real(ai[1]))
#             ox.append(xi)
#             oy.append(yi)
#         iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen3))
#
#         endPoints = [[0,0], [0,0]]
#         endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1] = calcVelWidth(x, y, ox[-1], oy[-1], False)
#         d = (0.5) * sqrt((endPoints[0][0] - endPoints[1][0]) ** 2 + (endPoints[0][1] - endPoints[1][1]) ** 2)
#         theta = np.arctan2(float(y - oy[-1]), float(x - ox[-1]))
#         x, y = endPoints[0][0] + (d * -np.sin(theta)), endPoints[0][1] + (d * np.cos(theta))
#
#
#     if shift:
#         intLine(e)
#     else:
#         if not vptSel:
#             for pt in vpts:
#                 if pt.checkClicked(e.pos()):
#                     vptSel = True
#                     vptCur = pt
#                     print 'Found it! at ', pt.getPos()
#             if not vptSel:
#                 if not first:
#                     x = e.pos().x()
#                     y = e.pos().y()
#                 px, py = projCoord(x,y)
#                 vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', x, y, d2=velocity.vy)
#                 surfaceInterp = getInterpolators(surface.data, surface.name, x, y)
#                 vxd = vxInterp([px], [py], grid=False)
#                 vyd = vyInterp([px], [py], grid=False)
#                 sur = surfaceInterp([px], [py], grid=False)
#                 v0 = sqrt(vxd**2 + vyd**2)
#                 vpts.append(vpt(x, y, v0, velocity.plotWidget)) # in map coordinates x<10018, y< 17946
#                 x = int(np.floor(x))
#                 y = int(np.floor(y))
#                 txt = 'Point ' + str(len(vpts)-1) + \
#                 ':\n=================\n' + \
#                 'v: ' +      "{:.3f}".format(velocity.data[y][x]) + \
#                 '\nbed: ' +  "{:.3f}".format(bed.data[y][x]) + \
#                 '\nsurf: ' + "{:.3f}".format(surface.data[y][x]) + \
#                 '\nSMB: ' +  "{:.3f}".format(smb.data[y][x]) + '\n\n'
#
#                 textOut.append(txt)
#                 iiContainer.currentWidget().addItem(vpts[-1].getCross()[0])
#                 iiContainer.currentWidget().addItem(vpts[-1].getCross()[1])
#                 if len(vpts) > 1:
#                     xa = [vpts[-1].getX(), vpts[-2].getX()]
#                     ya = [vpts[-1].getY(), vpts[-2].getY()]
#                     vpts[-2].setLine(pg.PlotDataItem(xa,ya,connect='all'))
#                     iiContainer.currentWidget().addItem(vpts[-2].getLine())#,pen=plotPen)
#         else:
#             vptSel = False


# mouseCoordinates = QtGui.QLabel('x:\ty:')
buttonBox.addWidget(mouseCoordinates)

# def mouseMoved(e):
#     global vptSel, vptCur, integrateLine, currentMap
#     # print 'mm 1'
#     if e.isExit() is False:
#         # print 'mm 2'
#         if vptSel and integrateLine is not None:
#             cData = integrateLine.curve.getData()
#             imin = curveDistance(e.pos().x(), e.pos().y(), cData)
#             if imin != -1:
#                 x = cData[0][imin]
#                 y = cData[1][imin]
#                 vptCur.updateCross(x,y)
#             else:
#                 vptCur.updateCross(e.pos().x(), e.pos().y())
#         # else:
#         x = int(np.floor(e.pos().x()))
#         y = int(np.floor(e.pos().y()))
#         if np.abs(x) <= 10018 and np.abs(y) <= 17946:
#             # print 'mm 3'
#             if currentMap == 0:
#                 # Interpolating every spot is too much, causes lag
#                 # vInt = getInterpolators(sqrt(velocity.vx**2 + velocity.vy**2), 'v', x,y)
#                 # px, py = projCoord(x,y)
#                 # vLocal = vInt([px],[py], grid=False)
#                 # print 'mm 30'
#                 t = 'v: ' + str(velocity.data[y][x]) + '\tbed: ' + str(bed.data[y][x]) + '\tsurf: ' + str(surface.data[y][x]) + '\tSMB: ' + str(smb.data[y][x])
#                 # t = t + '\t bed: ' + str(bed[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))])  + '\t surface: ' + str(surfaceVals[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))])
#                 # mouseCoordinates.setText(t)
#
#             elif currentMap == 1:
#                 # print 'mm 31'
#                 mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y) + '\t bed: ' + str(bed[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))]))
#             elif currentMap == 2:
#                 # print 'mm 32'
#                 mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y) + '\t surface: ' + str(surfaceVals[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))]))



# def ky(e):
#     print 'key pressed ', e.key()
#     print 'type: -' + str(e.type()) + '-'
#     #6 is press, 7 is release
#     if e.type() == 6:
#         print 'shift true'
#         shift = True
#     else:
#         shift = False
#         print 'shift false'
#     # 16777249 is shift
# shift = False

#fixme velocity.imageItem.hoverEvent = mouseMoved

velocity.imageItem.hoverEvent = mouseMoved
bed.imageItem.hoverEvent = mouseMoved
surface.imageItem.hoverEvent = mouseMoved


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
- Shift+click vpt to integrate from there

    MODEL
- Align extrapolated points with fenics mesh
- Allow modeling down integrated line


'''