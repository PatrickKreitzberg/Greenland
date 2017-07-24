import fenics as fc
import sys
import time

from pyqtgraph.Qt import QtCore  # , QtWidgets
from scipy.integrate import ode

from helper_files.classes.dataset import dataset
from helper_files.classes.pltPoints import *
from helper_files.cm import *
from helper_files.gui import *
from helper_files.math_functions import *
from helper_files.pens import *
from helper_files.profile_driver import runModel

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
startTime = time.time()
print 'Loading...'
#####################################################
####           CONSTANTS/GLOBAL VARIABLES       ####
#####################################################
currentMap = 0  # selects which data map to show [velocity, bed, surface]
dr = 150         # spatial resolution for bottom plot interpolation in meters
botPlot = False  # if bottom plot has been populated or not
inBed = False
inSMB = False
inSurface = False
clickedCurve = False
dataCMFileName = './data/dataCMValues.h5'
dataFileName   = './data/AllDataSets.h5'

vpts = [] #holds [x,y,v] values, where x,y are the coordinates and v is the velocity magnitude at those coordinates



#####################################################
####         CREATE INITIAL DATA SET(S)          ####
#####################################################

velocity = dataset('velocity', bpLegend, greyPlotPen, map=True)
bp.addItem(velocity.pathPlotItem)
iiContainer.addWidget(velocity.plotWidget)
iiContainer.setCurrentWidget(velocity.plotWidget)

smb = dataset('smb', bpLegend, redPlotPen, map=True)
print 'SMB: ', np.amin(smb.data), np.amax(smb.data) #-11493.3860928 6060.80339304
bp.addItem(smb.pathPlotItem)
iiContainer.addWidget(smb.plotWidget)

surface = dataset('surface', bpLegend, greenPlotPen, map=True)
bp.addItem(surface.pathPlotItem)
iiContainer.addWidget(surface.plotWidget)

bed = dataset('bed', bpLegend, bluePlotPen, map=True)
bp.addItem(bed.pathPlotItem)
iiContainer.addWidget(bed.plotWidget)

velocityWidth = dataset('velocitywidth', bpLegend, purplePlotPen)


def changeMap(index):
    '''
    Called when data drop down menu is changed.
    :param index:
    :return:
    '''

    # print iiContainer.currentWidget()
    vr = iiContainer.currentWidget().getPlotItem().getViewBox().viewRange()
    global currentMap
    if index == 0 and currentMap != 0:
        #velocity
        currentMap = 0
        iiContainer.setCurrentWidget(velocity.plotWidget)
        velocity.imageItem.mouseClickEvent = mouseClick
        velocity.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])
    elif index == 1 and currentMap != 1:
        #bed
        currentMap = 1
        iiContainer.setCurrentWidget(bed.plotWidget)
        bed.imageItem.mouseClickEvent = mouseClick
        bed.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])
    elif index == 2 and currentMap != 2:
        #surface
        currentMap = 2
        iiContainer.setCurrentWidget(surface.plotWidget)
        surface.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])
    elif index == 3 and currentMap != 3:
        #SMB
        currentMap = 3
        iiContainer.setCurrentWidget(smb.plotWidget)
        smb.plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1])


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
    vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', x1, y1, d2=velocity.vy)
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
            vxInterp, vyInterp = getInterpolators(velocity.vx, velocity.name, (x1 + (dr*dis * -np.sin(theta))), (y1 + (dr*dis * np.cos(theta))), d2=velocity.vy)
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
    # print 'getPro ', np.real(y[0]), ' ', np.real(y[1])
    mx, my = mapCoord(np.real(y[0]), np.real(y[1]))
    vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', np.floor(mx), np.floor(my), d2=velocity.vy)
    return np.array([t * (-vxInterp([y[0]], [y[1]], grid=False)), t * (-vyInterp([y[0]], [y[1]], grid=False))])

integrateLine = None

def _intMesh0(t,y):
    vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', [math.floor(y[0]) - 1, math.floor(y[1]) - 1],
                                          p1=[math.ceil(y[0]) + 1, math.ceil(y[1]) + 1], d2=velocity.vy)
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

def intMesh():
    rng = velocity.imageItem.getViewBox().viewRange()
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
    '''
    Shows the line which follows the path in the velocity field
    :param e:
    :return:
    '''
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
    global dr, bpLegend, dataLen, botPlot
    botPlot = True
    velValues = []
    xValues = []
    smbValues = []
    surfValues = []
    bedValues = []
    linePoints = [0]
    vwValues = []
    graphX = []

    ########################################
    ##    GATHER LOCAL INTERPOLATORS      ##
    ########################################
    mxx = max(pt.x for pt in vpts)
    mxy = max(pt.y for pt in vpts)
    mix = min(pt.x for pt in vpts)
    miy = min(pt.y for pt in vpts)

    surfaceInterp      = getInterpolators(surface.data, surface.name,   mix, miy, x1=mxx, y1=mxy)
    bedInterp          = getInterpolators(bed.data,     bed.name,       mix, miy, x1=mxx, y1=mxy)
    vxInterp, vyInterp = getInterpolators(velocity.vx,  velocity.name,  mix, miy, x1=mxx, y1=mxy, d2=velocity.vy)
    smbInterp          = getInterpolators(smb.data,     smb.name,       mix, miy, x1=mxx, y1=mxy)

    for i in range(1, len(vpts)):
        '''
        This part compares neighbor points to each other. 
        '''
        theta = np.arctan2(float(vpts[i].getY() - vpts[i - 1].getY()), float(vpts[i].getX() - vpts[i - 1].getX()))
        pvx0, pvy0 = projCoord(vpts[i-1].x, vpts[i-1].y)
        pvx1, pvy1 = projCoord(vpts[i].x, vpts[i].y)
        # distance = (sqrt((vpts[i].getY() - vpts[i - 1].getY()) ** 2 + (vpts[i].getX() - vpts[i - 1].getX()) ** 2))
        distance = sqrt((pvx1-pvx0)**2 + (pvy1-pvy0)**2)
        remainder = distance%dr
        xline = linspace(0, distance, distance/150, endpoint=True)  # * const makes it every 150/const meters
        '''
        #FIXME NEED TO CHANGE SO IT LINES UP WITH FENICS MESH
        '''
        
        yline = linspace(0, distance, distance/150, endpoint=False)  # * const makes it every 150/const meters
        # print 'xline: '
        # print xline
        # print 'yline: '
        # print yline
        linePoints.append(distance * (1 / dr) + linePoints[-1])

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
                # print 'dist: ', graphX[-1] + sqrt((px[-1]-px[-2])**2 + (py[-1]-py[-2])**2)
            elif len(graphX) == 0:
                graphX.append(0)
            else:
                graphX.append(graphX[-1])# + sqrt((px[-1]) ** 2 + (py[-1]) ** 2))
            #     print 'dis3: ', graphX[-1] + sqrt((px[-1]) ** 2 + (py[-1]) ** 2)

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
        bed.pathData     = np.array(bedValues[0])
        surface.pathData = np.array(surfValues[0])
        for i in range(1, len(bedValues)):
            bed.pathData     = np.append(bed.pathData, bedValues[i])
            surface.pathData = np.append(surface.pathData, surfValues[i])


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
            velocity.pathData      = np.array(velValues[0])
            smb.pathData           = np.array(smbValues[0])
            velocityWidth.pathData = np.array(vwValues[0])

            for i in range(1, len(velValues)):
                velocity.pathData      = np.append(velocity.pathData,        velValues[i])
                smb.pathData           = np.append(smb.pathData,      smbValues[i])
                velocityWidth.pathData = np.append(velocityWidth.pathData, vwValues[i])

    print 'graphx[-1]: ', graphX[-1]
    xd0, yd0 = projCoord(vpts[0].x, vpts[0].y)
    xd1, yd1 = projCoord(vpts[1].x, vpts[1].y)
    xd2, yd2 = projCoord(vpts[2].x, vpts[2].y)
    print 'dist: ', (sqrt(((xd1-xd0)**2 + (yd1-yd0)**2)) + sqrt(((xd2-xd1)**2 + (yd2-yd1)**2)))
    if botPlotBool:
        return linePoints, np.array(graphX)
    # else:
    #     return nbed, nsurf  # , linePoints


def calcBP():
    '''
    Calculate the data for bottom plot then populate the plot.
    :return:
    '''
    global dr, bpLegend, dataLen, botPlot
    if len(vpts) > 0:
        print 'plotting'
        # nbed, nsurf, nv, nsmb, nvelWidth, linePoints, graphX = interpolateData(True)
        linePoints, graphX = interpolateData(True)
        # velocity.pathPlotItem.setData(graphX     , velocity.pathData)
        # smb.pathPlotItem.setData(graphX          , smb.pathData)
        # velocityWidth.pathPlotItem.setData(graphX, velocityWidth.pathData)
        surface.pathPlotItem.setData(graphX      , surface.pathData)
        # bed.pathPlotItem.setData(graphX          , bed.pathData)
        pg.QtGui.QApplication.processEvents()

def arrows():
    ''''
    Plot velocity vectors
    '''
    print 'arrows'
    print velocity.imageItem.getViewBox().viewRange() # [y0, y1] [x0, x1]
    rngx, rngy = velocity.imageItem.getViewBox().viewRange()
    print rngx
    print rngy
    print 'Starting arrows'
    for i in range(int(rngx[0]),int(rngx[1])):
        for j in range(int(rngy[0]), int(rngy[1])):
            vdir = [velocity.vx[j][i], -velocity.vy[j][i]]
            vmag = sqrt(velocity.vx[j][i]**2 +  velocity.vy[j][i]**2)
            # theta = np.arctan2(vdir[1], vdir[0])
            iiContainer.currentWidget().addItem(pg.PlotDataItem([(i+0.5),(i+0.5 + vdir[0]/(1.5*vmag))], [(j+0.5),(j+0.5 + vdir[1]/(1.5*vmag))], pen=plotPen2))
            iiContainer.currentWidget().addItem(pg.PlotDataItem([i+0.5], [j+0.5]), pen=(255,255,255), symbolBrush=(255,0,0), symbolPen='w')
    print 'finished arrows'



def mouseMovedBP(evt):
    global bpLegend, dataLen
    print 'mouse moved bottom plot'
    if botPlot:
        print 'is not not not none'
        pos = evt[0]  ## using signal proxy turns original arguments into a tuple
        if bp.getPlotItem().sceneBoundingRect().contains(pos):
            mousePoint = bp.getPlotItem().vb.mapSceneToView(pos)
            index = int(mousePoint.x())
            if index > 0 and index < len(velocity.pathData):
                velocity.legendItem.setText('Velocity' + str(velocity.pathData[index]))
                smb.legendItem.setText('SMB ' + str(smb.pathData[index]))
                velocityWidth.legendItem.setText('Velocity Width ' + str(velocityWidth.pathData[index]))
                surface.legendItem.setTxt('Surface ele ' + str(surface.pathData[index]))
    # else:
        # print 'It is NONE'

def clearPoints():
    del vpts[:]
    velocity.pathPlotItem.clear()
    surface.pathPlotItem.clear()
    smb.pathPlotItem.clear()
    bed.pathPlotItem.clear()


def runModelButt():
    global dr
    if len(vpts) > 0:
        interpolateData(False)
        print 'bed data:'
        print bed.pathData
        print 'surface data:'
        print surface.pathData
        THICKLIMIT = 10.  # Ice is never less than this thick
        H = surface.pathData - bed.pathData
        surface.pathData[H <= THICKLIMIT] = bed.pathData[H <= THICKLIMIT]
        N = len(bed.pathData)
        mesh = fc.IntervalMesh(N-1, 0, dr*(N-1))
        #FIXME the intervalMesh is consistantly 150 between each datapoint this not true for the data being sent
        hdf_name = '/home/pat/research/latest_profile.h5'
        hfile = fc.HDF5File(mesh.mpi_comm(), hdf_name, "w")
        V = fc.FunctionSpace(mesh,"CG",1)
        functBed     = fc.Function(V, name="Bed")
        functSurface = fc.Function(V, name="Surface")
        functBed.vector()[:]     = bed.pathData
        functSurface.vector()[:] = surface.pathData
        hfile.write(functBed.vector(), "/bed")
        hfile.write(functSurface.vector(), "/surface")
        hfile.write(mesh, "/mesh")
        hfile.close()
        runModel(hdf_name)
    print 'done 3'



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
    print 'mouseclick ', e.pos().x(), ' ', e.pos().y()
    if len(vpts) == 0:
        first = True
        x = e.pos().x()
        y = e.pos().y()
        x0p, y0p = projCoord(x,y)
        y0 = np.array([x0p, y0p])
        t0, t1, dt = 0, 80, .1
        r = ode(getProfile).set_integrator('zvode', method='bdf')
        r.set_initial_value(y0, t0)
        ox = []
        oy = []
        ind = 0
        while r.successful() and ind < 2:
            '''
            what the fuck was I trying to do here??
            '''
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
                vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', x, y, d2=velocity.vy)
                surfaceInterp = getInterpolators(surface.data, surface.name, x, y)
                vxd = vxInterp([px], [py], grid=False)
                vyd = vyInterp([px], [py], grid=False)
                sur = surfaceInterp([px], [py], grid=False)
                v0 = sqrt(vxd**2 + vyd**2)
                vpts.append(vpt(x, y, v0, velocity.plotWidget)) # in map coordinates x<10018, y< 17946
                x = int(np.floor(x))
                y = int(np.floor(y))
                txt = 'Point ' + str(len(vpts)-1) + \
                ':\n====================\n' + \
                'v: ' + "{:.3f}".format(velocity.data[y][x]) + \
                '\nbed: ' + "{:.3f}".format(bed.data[y][x]) + \
                '\nsurf: ' + "{:.3f}".format(surface.data[y][x]) + \
                '\nSMB: ' + "{:.3f}".format(smb.data[y][x]) + '\n\n'

                textOut.append(txt)
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
                # Interpolating every spot is too much, causes lag
                # vInt = getInterpolators(sqrt(velocity.vx**2 + velocity.vy**2), 'v', x,y)
                # px, py = projCoord(x,y)
                # vLocal = vInt([px],[py], grid=False)
                # print 'mm 30'
                t = 'v: ' + str(velocity.data[y][x]) + '\tbed: ' + str(bed.data[y][x]) + '\tsurf: ' + str(surface.data[y][x]) + '\tSMB: ' + str(smb.data[y][x])
                # t = t + '\t bed: ' + str(bed[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))])  + '\t surface: ' + str(surfaceVals[int(np.ceil(np.abs(y)))][int(np.floor(np.abs(x)))])
                # mouseCoordinates.setText(t)

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
#fixme velocity.imageItem.hoverEvent = mouseMoved

velocity.imageItem.hoverEvent = mouseMoved
bed.imageItem.hoverEvent = mouseMoved
surface.imageItem.hoverEvent = mouseMoved


velocity.imageItem.mouseClickEvent = mouseClick

# proxy = pg.SignalProxy(bp.getPlotItem().scene().sigMouseMoved, rateLimit=60, slot=mouseMovedBP)

#fixme velocity.imageItem.mouseClickEvent = mouseClick # Default map is vel FIXME

# iiContainer.keyboardGrabber()
# iiContainer.keyPressEvent = ky
# iiContainer.keyReleaseEvent = ky


print 'Loaded in: ', time.time() - startTime, ' seconds!'

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