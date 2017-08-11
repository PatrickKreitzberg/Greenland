from gui import *
from dataset_objects import *
from math_functions import *
import time
import pickle
import os
import scipy.signal as signal
from peakdetect import *



dr = 150

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
    calcVelWidth(vpts[0].cx, vpts[0].cy, vpts[1].cx, vpts[1].cy, True)

    for i in range(1, len(vpts)):
        # if not vpts[i][2]:
        calcVelWidth(vpts[i - 1].cx, vpts[i - 1].cy, vpts[i].cx, vpts[i].cy, True)

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


# def getInterpolators(d1, choice, x0, y0, x1=-99, y1=-99, d2=None):
#     '''
#     Determines the local interpolator and returns it.
#     Interpolates data d1 (and d2 if necessary)
#     If x1,y1 not specified then creates a local interpolator of size 10x10 with
#     the point x0,y0 in the middle.
#
#     :param d1:
#     :param choice:  Which dataset to process
#     :param x0:  IN MAP COORDINATES
#     :param y0:  IN MAP COORDINATES
#     :param x1:
#     :param y1:
#     :param d2:
#     :return:
#     '''
#     # vel_x0 = -638000  # first x coordinate
#     # vel_x1 = 864550  # last x coordinate
#     # vel_y0 = -657600  # first y coordinate
#     # vel_y1 = -3349350  # last y coordinate
#     # vel_xarray = linspace(vel_x0, vel_x1, 10018, endpoint=True)
#     # vel_yarray = linspace(vel_y1, vel_y0, 17946, endpoint=True)
#
#     # Make sure points are in bounds
#     p0 = [x0, y0]
#     p1 = [x1, y1]
#     p0[0], p0[1] = np.floor(p0[0]), np.floor(p0[1])
#
#     if p0[0] < 0:
#         p0[0] = 0
#     if p0[0] > map['x1']:
#         p0[0] = map['x1']
#     if p0[1] < 0:
#         p0[1] = 0
#     if p0[1] > map['y1']:
#         p0[1] = map['y1']
#     if p1[0] == -99:
#         # if the function is sent a single point then create interpolator around the point
#         # in this case a 10x10 interpolator
#         p1[0] = p0[0] - 5
#         p1[1] = p0[1] - 5
#         p0[0] = p0[0] + 5
#         p0[1] = p0[1] + 5
#     else:
#         p1[0], p1[1] = np.floor(p1[0]), np.floor(p1[1])
#         if p1[0] < 0:
#             p1[0] = 0
#         if p1[0] > map['x1']:
#             p1[0] = map['x1']
#         if p1[1] < 0:
#             p1[1] = 0
#         if p1[1] > map['y1']:
#             p1[1] = map['y1']
#
#
#     minSpacing = 10
#     # p1, p2 in map coordinates
#     projx0, projy0 = projCoord(p0[0], p0[1])
#     projx1, projy1 = projCoord(p1[0], p1[1])
#     #FIXME Should there be a minimum dx, dy?
#     dx = 1 + math.fabs(p1[0] - p0[0])
#     dy = 1 + math.fabs(p1[1] - p0[1])
#
#     if p0[1] < p1[1]:
#         vel_yarray = linspace(projy1, projy0, int(dy), endpoint=True)
#         y0 = p0[1]
#         y1 = p1[1] + 1
#     else:
#         vel_yarray = linspace(projy0, projy1, int(dy), endpoint=True)
#         y0 = p1[1]
#         y1 = p0[1] + 1
#     if p0[0] < p1[0]:
#         vel_xarray = linspace(projx0, projx1, int(dx), endpoint=True)
#         x0 = p0[0]
#         x1 = p1[0] + 1
#     else:
#         vel_xarray = linspace(projx1, projx0, int(dx), endpoint=True)
#         x0 = p1[0]
#         x1 = p0[0] + 1
#     y0, y1, x0, x1 = int(y0), int(y1), int(x0), int(x1)
#     possibleChoices = ['velocity', 'bed', 'surface', 'smb', 'thickness']
#     if choice == 'vxvy':
#         return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose()), \
#                RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d2[y0:y1, x0:x1])).transpose())
#     elif choice in possibleChoices:
#         return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
#     # elif choice is 'bed':
#     #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
#     # elif choice is 'surface':
#     #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
#     # elif choice is 'smb':
#     #     #smb and all other RACMO data is 'upside down' compared to the rest of the data
#     #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
#     # elif choice is 'v':
#     #     #just velocity, not vx and vy
#     #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
#     else:
#         print "ERROR: No interpolator selected!  ./helper_files/data_functions.getInterpolators()"
#
#
#
#     # RectBivariateSpline(vel_xarray, vel_yarray, (npk.flipud(vy)).transpose()) #FIXME SHOULD IT BE FLIPPED!?!?!?
#     #FIXME Make sure flipping the transpose works!!!
#     #FIXME changed the interpolator


def interpolateData(runModel):
    '''
    Calculate the data for bottom plot or to run the model.
    If botPlotBool, calculate all the data.  Else, calculate just bed/surface.
    :return:
    '''
    dr = float(model_res_lineEdit.text()) # dr = 150
    velValues = []
    xValues = []
    smbValues = []
    surfValues = []
    bedValues = []
    thickValues = []
    linePoints = [0]
    vwValues = []
    graphX = []
    ########################################
    ##    GATHER LOCAL INTERPOLATORS      ##
    ########################################
    '''
    mxx = max(pt.x for pt in vpts)
    mxy = max(pt.y for pt in vpts)
    mix = min(pt.x for pt in vpts)
    miy = min(pt.y for pt in vpts)

    
    if velocityCheck.checkState() == 2 or runModel:
        # vxInterp, vyInterp = getInterpolators(velocity.vx, velocity.name, mix, miy, x1=mxx, y1=mxy, d2=velocity.vy)
        velInterp = getInterpolators(velocity.data, 'velocity', mix, miy, x1=mxx, y1=mxy)
    if surfaceCheck.checkState() == 2 or runModel:
        surfaceInterp = getInterpolators(surface.data, surface.name, mix, miy, x1=mxx, y1=mxy)
    if bedCheck.checkState() == 2 or runModel:
        bedInterp = getInterpolators(bed.data, bed.name, mix, miy, x1=mxx, y1=mxy)
    if smbCheck.checkState() == 2 or runModel:
        smbInterp = getInterpolators(smb.data, smb.name, mix, miy, x1=mxx, y1=mxy)
    thickInterp   = getInterpolators(thickness.data, thickness.name, mix, miy, x1=mxx, y1=mxy)
    '''

    # thickInterp = getInterpolators(thickness.data, thickness.name, mix, miy, x1=mxx, y1=mxy)

    # Find a distance ~150m which gets close to dividing the distance between first 2 spots
    d = 0
    for i in range(1, len(vpts)):
        d += sqrt((vpts[i - 1].cx - vpts[i].cx) ** 2 + (vpts[i - 1].cy - vpts[i].cy) ** 2)
    # print 'distance is: ', d*150, 'divide', int(d*(150/dr))

    for i in range(1, len(vpts)):
        '''
        This part compares neighbor points to each other. 
        '''

        theta = np.arctan2(float(vpts[i].cy - vpts[i - 1].cy), float(vpts[i].cx - vpts[i - 1].cx))
        # pvx0, pvy0 = projCoord(vpts[i-1].x, vpts[i-1].y)
        # pvx1, pvy1 = projCoord(vpts[i].x, vpts[i].y)
        # distance = (sqrt((vpts[i].getY() - vpts[i - 1].getY()) ** 2 + (vpts[i].getX() - vpts[i - 1].getX()) ** 2))
        distance = sqrt((vpts[i - 1].cx - vpts[i].cx) ** 2 + (vpts[i - 1].cy - vpts[i].cy) ** 2)
        # remainder = distance % dr
        # Xline needs to be in map coordinates because tx,ty are in map coordinates
        xline = linspace(0, distance, int(distance*(150/dr)), endpoint=True)  # * 1/dr*150 makes the resolution dr

        '''
        #FIXME NEED TO CHANGE SO IT LINES UP WITH FENICS MESH
        '''

        # Rotation matrix:
        rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])
        px = []  # px, py are projected coordinates used to get values from the interpolator.  Projected meaning IRL values
        py = []



        for j in range(len(xline)):
            # rotate coordinates
            t = rotMatrix * np.matrix([[xline[j]], [0.0]])
            # transform coordinates into projected coordinates
            tx, ty = colorToProj(vpts[i - 1].cx + t[0, 0], vpts[i - 1].cy + t[1, 0])
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
        if surfaceCheck.checkState() == 2 or runModel:
            # surfaceInterp = getInterpolators(surface.data, surface.name, mix, miy, x1=mxx, y1=mxy)
            localSurface = surface.interp(px, py, grid=False)
            surfValues.append(localSurface)

        ########################################
        ##         CALCULATE BED              ##
        ########################################
        if bedCheck.checkState() == 2 or runModel:
            # bedInterp = getInterpolators(bed.data, bed.name, mix, miy, x1=mxx, y1=mxy)
            localBed = bed.interp(px, py, grid=False)
            bedValues.append(localBed)

        ########################################
        ##        CALCULATE VELOCITY          ##
        ########################################
        if velocityCheck.checkState() == 2 or runModel:
            # velInterp = getInterpolators(velocity.data, 'velocity', mix, miy, x1=mxx, y1=mxy)
            vi = velocity.interp(px, py, grid=False)
            xValues.append(xline)
            velValues.append(vi)

        ########################################
        ##     CALCULATE VELOCITY WIDTH       ##
        ########################################
        if vWidthCheck.checkState() == 2 or runModel:
            vwd = []
            for i in range(len(px)):
                xp0, yp0 = colorCoord(px[i - 1], py[i - 1])
                xp1, yp1 = colorCoord(px[i], py[i])
                xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False)
                vwd.append(sqrt((xril - xrir) ** 2 + (yril - yrir) ** 2))
            vwValues.append(vwd)

        ########################################
        ##   CALCULATE SURFACE MASS-BALANCE   ##
        ########################################
        if smbCheck.checkState() == 2 or runModel:
            # smbInterp = getInterpolators(smb.data, smb.name, mix, miy, x1=mxx, y1=mxy)
            'init smb'
            localSMB = smb.interp(px, py, grid=False)
            smbValues.append(localSMB)

        ########################################
        ##   CALCULATE THICKNESS              ##
        ########################################
        # if smbCheck == 2 or runModel:
        # thickInterp   = getInterpolators(thickness.data, thickness.name, mix, miy, x1=mxx, y1=mxy)
        if thicknessCheck.checkState() == 2 or runModel:
            localThick = thickness.interp(px, py, grid=False)
            thickValues.append(localThick)

        ########################################
        ##   COMPILE DATA                     ##
        ########################################
        if bedCheck.checkState() == 2 or runModel:
            bed.pathData = np.array(bedValues[0])
        if surfaceCheck.checkState() == 2 or runModel:
            surface.pathData = np.array(surfValues[0])
        if velocityCheck.checkState() == 2 or runModel:
            velocity.pathData = np.array(velValues[0])
        if smbCheck.checkState() == 2 or runModel:
            smb.pathData = np.array(smbValues[0])
        if vWidthCheck.checkState() == 2 or runModel:
            velocityWidth.pathData = np.array(vwValues[0])
        if thicknessCheck.checkState() == 2 or runModel:
            thickness.pathData = np.array(thickValues[0])

        for i in range(1, len(velValues)):
            if velocityCheck.checkState() == 2 or runModel:
                velocity.pathData      = np.append(velocity.pathData, velValues[i])
            if smbCheck.checkState() == 2 or runModel:
                smb.pathData           = np.append(smb.pathData, smbValues[i])
            if vWidthCheck.checkState() == 2 or runModel:
                velocityWidth.pathData = np.append(velocityWidth.pathData, vwValues[i])
            if bedCheck.checkState() == 2 or runModel:
                bed.pathData           = np.append(bed.pathData, bedValues[i])
            if surfaceCheck.checkState() == 2 or runModel:
                surface.pathData       = np.append(surface.pathData, surfValues[i])
            if thicknessCheck.checkState() == 2 or runModel:
                thickness.pathData = np.append(thickness.pathData, thickValues[i])
    if smbCheck.checkState() == 2 or runModel:
        smb.pathData = smb.pathData*(1.0/1000.0)*(916.7/1000.0) # millimeters -> meters then water-equivalent to ice-equivalent
    # print 'graphx[-1]: ', graphX[len(graphX) - 1]
    dist = 0
    for i in range(len(vpts) - 1):
        # print 'calc distance...'
        xd0, yd0 = colorToProj(vpts[i].cx, vpts[i].cy)
        xd1, yd1 = colorToProj(vpts[i + 1].cx, vpts[i + 1].cy)
        # print xd0, yd0, xd1, yd1
        dist += sqrt(((xd1 - xd0) ** 2 + (yd1 - yd0) ** 2))
    # print 'dist: ', dist
    if thicknessCheck.checkState() == 2 or runModel:
        thickness.distanceData = graphX
    if surfaceCheck.checkState() == 2 or runModel:
        surface.distanceData = graphX
    if bedCheck.checkState() == 2 or runModel:
        bed.distanceData = graphX
    if velocityCheck.checkState() == 2 or runModel:
        velocity.distanceData = graphX
    if smbCheck.checkState() == 2 or runModel:
        smb.distanceData = graphX
    if vWidthCheck.checkState() == 2 or runModel:
        velocityWidth.distanceData = graphX

