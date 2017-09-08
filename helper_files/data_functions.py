import time
import os
from peakdetect import *
from scipy import sqrt
import numpy as np
import pyqtgraph as pg


# LOCAL IMPORTS
from gui import *
from dataset_objects import *
from math_functions import *
from pens import whitePlotPen


# dr = 150

def getProfile(t,y):
    '''
    Prints a line
    :param t:
    :param y:
    :return:
    '''
    # print 'getPro ', np.real(y[0]), ' ', np.real(y[1])
    # mx, my = mapCoord(np.real(y[0]), np.real(y[1]))
    return np.array([t * (-velocity.vxInterp([y[0]], [y[1]], grid=False)), t * (-velocity.vyInterp([y[0]], [y[1]], grid=False))])


def cwLoop(e):
    '''
    Calls the calculate width function in order to get the width along a profile
    :param e:
    :return:
    '''
    calcVelWidth(markers[0].cx, markers[0].cy, markers[1].cx, markers[1].cy, True)

    for i in range(1, len(markers)):
        # if not vpts[i][2]:
        calcVelWidth(markers[i - 1].cx, markers[i - 1].cy, markers[i].cx, markers[i].cy, True)

def gaussian(x, A, x0, sig):
    return A*math.exp(-(x-x0)**2/(2.0*sig**2))

def fit(p,x):
    return np.sum([gaussian(x, p[i*3],p[i*3+1],p[i*3+2])
                   for i in xrange(len(p)/3)],axis=0)

# def calcVelWidth(x0, y0, x1, y1, draw):
#     '''
#     Calculates the width of the ice stream at one point, (x1, y1).  (x0, y0) is there
#     to give an idea of where the velocity width begins and ends which should be on a
#     line which is perpindicular to the line from (x1, y1) to (x0, y0).
#     :param x0: color coord
#     :param y0: color coord
#     :param x1: color coord
#     :param y1: color coord
#     :param draw:
#     :return:
#     '''
#     staticPlotWindow = QtGui.QMainWindow(mw)
#     staticPlotWindow.setWindowTitle('Static Plotter')
#     dummyWidget = QtGui.QWidget()
#     staticPlotWindow.setCentralWidget(dummyWidget)
#     layout = QtGui.QVBoxLayout()
#     dummyWidget.setLayout(layout)
#     plt1 = pg.PlotWidget()
#     layout.addWidget(plt1)
#
#
#     # input is in map coordinates
#     #
#     # This is with interpolation
#     #
#     theta = np.arctan2(float(y1 - y0), float(x1-x0))
#     rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])
#     ls = linspace(-1,1,2,endpoint=True)
#
#     # lines perp to path
#     tn = rotMatrix * np.matrix([[0.0], [-10]])
#     tp = rotMatrix * np.matrix([[0.0], [10]])
#
#     tx1, ty1 = colorToProj(x1,y1)
#     v0 = velocity.interp(tx1,ty1, grid=False)
#     dv = [[0, 0, 0], [0, 0, 0]]   # x, y, dv for left and right
#     endPoints = [[0, 0], [0, 0]]  # end points [left[x,y], right[x,y]]
#     # print 'v0 ', vx0, vy0, v0
#
#     txn, tyn = colorToProj(x1 + tn[0], y1 + tn[1])  # dis either + or negative
#     txp, typ = colorToProj(x1 + tp[0], y1 + tp[1])
#     vn = velocity.interp(txn, tyn, grid=False)
#     vp = velocity.interp(txp, typ, grid=False)
#     dwidth = 550
#     if vn > v0: #head in negative direction
#         path = linspace(0, dwidth, dwidth/10, endpoint=True)  # would be every 10 meters this in proj coord
#         t = rotMatrix * np.matrix([[0] * len(path), path])  # proj coord
#
#         vn = velocity.interp(t[0]*-150 + tx1, t[1]*-150 + ty1, grid=False)
#
#         peaks = peakdetect(vn[0], path*150, lookahead=5)
#         print 'up'
#         print vn
#         print 'peaks', peaks
#         print 'path', path
#         plt1.getPlotItem().plot(path*150, vn[0])
#         for p in peaks[0]:
#             plt1.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0]))
#         endPoints = [[x1, y1], [x1 - t.item((0, t.shape[1]-1)), y1 - t.item((1, t.shape[1]-1))]]
#         #     px, py = colorToProj(x1, y1)
#         #     # pathx = linspace(px, )
#         #     tx1, ty1, = txn + (dr * -1 * -np.sin(theta)), tyn + (dr * -1 * np.cos(theta))
#     else:
#         path = linspace(0, dwidth, dwidth/10, endpoint=True)  # would be every 10 meters
#         t = rotMatrix * np.matrix([[0] * len(path), path])  # proj coord
#         print 'down'
#
#         vp = velocity.interp(t[0] * 150 + tx1, t[1] * 150 + ty1, grid=False)
#         print vp
#         peaks = peakdetect(vp[0], path*150, lookahead=5)
#         print 'peaks', peaks
#         print 'path', path
#         plt1.getPlotItem().plot(path*150, vp[0])
#         for p in peaks[0]:
#             plt1.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0]))
#         print t
#         print t.shape
#         endPoints = [[x1, y1], [x1 + t.item((0, t.shape[1]-1)), y1 + t.item((1, t.shape[1]-1))]]
#     print endPoints
#     draw = True
#     if draw:
#         iiContainer.currentWidget().addItem(pg.PlotDataItem([endPoints[0][0], endPoints[1][0]], [endPoints[0][1], endPoints[1][1]], connect='all', pen=whitePlotPen))
#         # circle plotting
#         d = (0.5)*sqrt((endPoints[0][0]-endPoints[1][0])**2 + (endPoints[0][1]-endPoints[1][1])**2)
#         cax, cay = endPoints[1][0] + (d * -np.sin(theta)), endPoints[1][1] + (d * np.cos(theta))
#         xc, yc = circArr(cax, cay)
#         iiContainer.currentWidget().addItem(pg.PlotDataItem(xc, yc, connect='all', pen=blackPlotPen))
#         iiContainer.currentWidget().addItem(pg.PlotDataItem([cax], [cay], pen=blackPlotPen))
#     staticPlotWindow.show()
#     return endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]

def calcVelWidth(x0, y0, x1, y1, draw):
    '''
    Calculates the width of the ice stream at one point, (x1, y1).  (x0, y0) is there
    to give an idea of where the velocity width begins and ends which should be on a
    line which is perpindicular to the line from (x1, y1) to (x0, y0).
    :param x0: color coord
    :param y0: color coord
    :param x1: color coord
    :param y1: color coord
    :param draw:
    :return:
    '''

    # input is in map coordinates
    #
    # This is with interpolation
    #
    theta = np.arctan2(float(y1 - y0), float(x1-x0))
    tx1, ty1 = colorToProj(x1,y1)
    v0 = velocity.interp(tx1,ty1, grid=False) # integrate for a bit to find which direction to look for boundary (perp to integration line)
    dv = [[0, 0, 0], [0, 0, 0]]   # x, y, dv for left and right
    endPoints = [[0, 0], [0, 0]]  # end points [left[x,y], right[x,y]]
    # print 'v0 ', vx0, vy0, v0
    vOverTTot = []
    distance = []
    # staticPlotWindow = QtGui.QMainWindow(mw)
    # staticPlotWindow.setWindowTitle('Static Plotter')
    # dummyWidget = QtGui.QWidget()
    # staticPlotWindow.setCentralWidget(dummyWidget)
    # layout = QtGui.QVBoxLayout()
    # dummyWidget.setLayout(layout)
    # plt1 = pg.PlotWidget()
    # layout.addWidget(plt1)
    for i in range(2):
        # one loop for each boundary
        vOverT = []
        dist = []
        dr = 0
        currentVelocity = 10
        startEndRatio = 0
        vOld = v0
        minReached = False
        minIndex = 0
        minVel = 9999
        run = True
        if i == 0:
            dis = 1
        else:
            dis = -1
        # print 'min([int(v0%100),8]) ', min([int(v0%100),8])
        while currentVelocity > 5 and startEndRatio <= min([int(v0%100),8]) and run:
            # this while loop determines when the edge is reached
            dr += 1
            tx, ty = colorToProj(x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta)))  # Line perpindicular to flow
            dist.append(dis*sqrt((dr*dis * -np.sin(theta))**2 + (dr*dis * np.cos(theta))**2))
            currentVelocity = velocity.interp(tx, ty, grid=False)
            if currentVelocity < minVel:
                minVel = currentVelocity
                minIndex = dr - 1
                if minVel < 1/3:
                    minReached = True
            vOverT.append(currentVelocity)
            if np.abs(currentVelocity - vOld) > dv[i][2]:
                dv[i][0], dv[i][1] = x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta))
                dv[i][2] = np.abs(currentVelocity - vOld)
            if currentVelocity !=0:
                startEndRatio = v0/currentVelocity
            if minReached and currentVelocity >= v0:
                #
                # If velocity got much smaller than equally as big should end.  Means there is probably another stream very close.
                #
                run = False
                dist = dist[:minIndex]
                vOverT = vOverT[:minIndex]
            vOld = currentVelocity


        if i == 0:
            for j in range(len(vOverT)):
                vOverTTot.append(vOverT[::-1][j])
                distance.append(dist[::-1][j]*150)

            # vOverTTot.append(vOverT[::-1][:])
        else:
            for j in range(len(vOverT)):
                vOverTTot.append(vOverT[j])
                distance.append(dist[j]*150)
            # vOverTTot.append(vOverT[:])
        if currentVelocity < 5:
            # plotting line
            endPoints[i][0], endPoints[i][1] = dv[i][0], dv[i][1]#mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
        else:
            endPoints[i][0], endPoints[i][1] = colorCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
    lah = 6 #int(700*len(distance)/math.fabs(distance[0]-distance[-1]))
    # print 'look ahead ', lah
    # print len(distance)
    # print 'distance', math.fabs(distance[0]-distance[-1])
    # print 'vou', vOverTTot
    # peaks = peakdetect(vOverTTot, distance, lookahead=lah, delta=np.amax(vOverTTot)/5)
    # print 'peaks', peaks
    # plt1.getPlotItem().plot(distance, vOverTTot)
    # if len(peaks[0]) > 0:
    #     for p in peaks[0]:
    #         plt1.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0]))
    # if len(peaks[1]) > 0:
    #     for p in peaks[1]:
    #         plt1.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0], pen=redPlotPen))
    # staticPlotWindow.show()
    # draw=True
    if draw:
        iiContainer.currentWidget().addItem(pg.PlotDataItem([endPoints[0][0], endPoints[1][0]], [endPoints[0][1], endPoints[1][1]], connect='all', pen=whitePlotPen))
        # circle plotting
        d = (0.5)*sqrt((endPoints[0][0]-endPoints[1][0])**2 + (endPoints[0][1]-endPoints[1][1])**2)
        cax, cay = endPoints[1][0] + (d * -np.sin(theta)), endPoints[1][1] + (d * np.cos(theta))
        xc, yc = circArr(cax, cay)
        iiContainer.currentWidget().addItem(pg.PlotDataItem(xc, yc, connect='all', pen=blackPlotPen))
        iiContainer.currentWidget().addItem(pg.PlotDataItem([cax], [cay], pen=blackPlotPen))
    return endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]


def interpolateData(runModel, dr, dataSetsToPopulate):
    '''
    Populates the data (velocity, thickness, etc.) along the path.
    :parameter dataSetsToPopulate is a dictionary of the datasets to interpolate.
    '''
    # dr = float(model_res_lineEdit.text()) # dr = 150
    velValues = []
    xValues = []
    smbValues = []
    surfValues = []
    bedValues = []
    thickValues = []
    linePoints = [0]
    vwValues = []
    graphX = []

    # d = 0
    # for i in range(1, len(markers)):
    #     # Calculate total path distance (this is in units of 150meters
    #     d += sqrt((markers[i - 1].cx - markers[i].cx) ** 2 + (markers[i - 1].cy - markers[i].cy) ** 2)

    # print 'distance is: ', d*150, 'divide', int(d*(150/dr))

    for i in range(1, len(markers)):
        '''
        This part compares neighbor points to each other. 
        '''

        # For each two points (i and the pt before i) get an array of points between the two,
        # roughly dr distance between the 2
        theta = np.arctan2(float(markers[i].cy - markers[i - 1].cy), float(markers[i].cx - markers[i - 1].cx))
        distance = sqrt((markers[i - 1].cx - markers[i].cx) ** 2 + (markers[i - 1].cy - markers[i].cy) ** 2)
        xline = linspace(0, distance, int(distance*(150/dr)), endpoint=True)  # * 1/dr*150 makes the resolution dr


        # Rotation matrix:
        # Takes points in a line and puts them along the path between the 2 markers
        rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)],
                               [np.sin(theta),      np.cos(theta)]])
        px = []  # px, py are projected coordinates used to get values from the interpolator.  Projected meaning IRL values
        py = []


        for j in range(len(xline)):
            t = rotMatrix * np.matrix([[xline[j]], [0.0]])
            # transform coordinates into projected coordinates
            # (the interpolators are set to use projected coordinates)
            tx, ty = colorToProj(markers[i - 1].cx + t[0, 0], markers[i - 1].cy + t[1, 0])

            # The following was used for debugging, prints marks at every point in tx, ty
            # c = 3
            # iiContainer.currentWidget().addItem(
            #     pg.PlotDataItem([markers[i - 1].cx + t[0, 0] - c,
            #                      markers[i - 1].cx + t[0, 0] + c],
            #                     [markers[i - 1].cy + t[1, 0] - c,
            #                      markers[i - 1].cy + t[1, 0] + c],
            #                     connect='all', pen=whitePlotPen))
            # iiContainer.currentWidget().addItem(
            #     pg.PlotDataItem([markers[i - 1].cx + t[0, 0] - c,
            #                      markers[i - 1].cx + t[0, 0] + c],
            #                     [markers[i - 1].cy + t[1, 0] - c,
            #                      markers[i - 1].cy + t[1, 0] + c],
            #                     connect='all', pen=whitePlotPen))

            px.append(tx)
            py.append(ty)
            if len(px) > 1:
                graphX.append(graphX[len(graphX) - 1] + sqrt((px[-1] - px[-2]) ** 2 + (py[-1] - py[-2]) ** 2))
                # print 'dist: ', graphX[-1] + sqrt((px[-1]-px[-2])**2 + (py[-1]-py[-2])**2)
            elif len(graphX) == 0:
                graphX.append(0)
            else:
                graphX.append(graphX[len(graphX) - 1])  # + sqrt((px[-1]) ** 2 + (py[-1]) ** 2))
                #     print 'dis3: ', graphX[-1] + sqrt((px[-1]) ** 2 + (py[-1]) ** 2)

        ########################################
        ##    CALCULATE SURFACE ELEVATION     ##
        ########################################
        if runModel or 'sur' in dataSetsToPopulate:
            localSurface = surface.interp(px, py, grid=False)
            surfValues.append(localSurface)

        ########################################
        ##         CALCULATE BED              ##
        ########################################
        if runModel or 'bed' in dataSetsToPopulate:
            localBed = bed.interp(px, py, grid=False)
            bedValues.append(localBed)

        ########################################
        ##        CALCULATE VELOCITY          ##
        ########################################
        if runModel or 'vel' in dataSetsToPopulate:
            vi = velocity.interp(px, py, grid=False)
            xValues.append(xline)
            velValues.append(vi)

        ########################################
        ##     CALCULATE VELOCITY WIDTH       ##
        ########################################
        if runModel or 'wth' in dataSetsToPopulate:
            vwd = []
            for i in range(len(px)):
                xp0, yp0 = colorCoord(px[i - 1], py[i - 1])
                xp1, yp1 = colorCoord(px[i], py[i])
                xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False)
                vwd.append(sqrt((xril - xrir) ** 2 + (yril - yrir) ** 2))
            vwValues.append(vwd)

        ########################################
        ##   CALCULATE SURFACE MASS BALANCE   ##
        ########################################
        if runModel or 'smb' in dataSetsToPopulate:
            'init smb'
            localSMB = smb.interp(px, py, grid=False)
            smbValues.append(localSMB)

        ########################################
        ##   CALCULATE THICKNESS              ##
        ########################################
        if runModel or 'thk' in dataSetsToPopulate:
            localThick = thickness.interp(px, py, grid=False)
            thickValues.append(localThick)

        ########################################
        ##   COMPILE DATA                     ##
        ########################################
        if runModel or 'bed' in dataSetsToPopulate:
            bed.pathData = np.array(bedValues[0])

        if runModel or 'sur' in dataSetsToPopulate:
            surface.pathData = np.array(surfValues[0])

        if runModel or 'vel' in dataSetsToPopulate:
            velocity.pathData = np.array(velValues[0])

        if runModel or 'smb' in dataSetsToPopulate:
            smb.pathData = np.array(smbValues[0])

        if runModel or 'wth' in dataSetsToPopulate:
            velocityWidth.pathData = np.array(vwValues[0])

        if runModel or 'thk' in dataSetsToPopulate:
            thickness.pathData = np.array(thickValues[0])

        for i in range(1, len(velValues)):
            if runModel or 'vel' in dataSetsToPopulate:
                velocity.pathData      = np.append(velocity.pathData, velValues[i])

            if runModel or 'smb' in dataSetsToPopulate:
                smb.pathData           = np.append(smb.pathData, smbValues[i])

            if runModel or 'wth' in dataSetsToPopulate:
                velocityWidth.pathData = np.append(velocityWidth.pathData, vwValues[i])

            if runModel or 'bed' in dataSetsToPopulate:
                bed.pathData           = np.append(bed.pathData, bedValues[i])

            if runModel or 'sur' in dataSetsToPopulate:
                surface.pathData       = np.append(surface.pathData, surfValues[i])

            if runModel or 'thk' in dataSetsToPopulate:
                thickness.pathData = np.append(thickness.pathData, thickValues[i])

    if runModel or 'smb' in dataSetsToPopulate:
        smb.pathData = smb.pathData*(1.0/1000.0)*(916.7/1000.0) # millimeters -> meters then water-equivalent to ice-equivalent

    dist = 0
    for i in range(len(markers) - 1):

        xd0, yd0 = colorToProj(markers[i].cx, markers[i].cy)
        xd1, yd1 = colorToProj(markers[i + 1].cx, markers[i + 1].cy)

        dist += sqrt(((xd1 - xd0) ** 2 + (yd1 - yd0) ** 2))

    if runModel or 'thk' in dataSetsToPopulate:
        thickness.distanceData = graphX

    if runModel or 'sur' in dataSetsToPopulate:
        surface.distanceData = graphX

    if runModel or 'bed' in dataSetsToPopulate:
        bed.distanceData = graphX

    if runModel or 'vel' in dataSetsToPopulate:
        velocity.distanceData = graphX

    if runModel or 'smb' in dataSetsToPopulate:
        smb.distanceData = graphX

    if runModel or 'wth' in dataSetsToPopulate:
        velocityWidth.distanceData = graphX

