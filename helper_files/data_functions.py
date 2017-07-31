from gui import *
from dataset_objects import *
from math_functions import *
import time
import pickle
import os
# os.chdir('..')
print 'path d: ', os.getcwd()




dr = 150

def getProfile(t,y):
    '''
    Prints a line
    :param t:
    :param y:
    :return:
    '''
    # print 'getPro ', np.real(y[0]), ' ', np.real(y[1])
    mx, my = mapCoord(np.real(y[0]), np.real(y[1]))
    vxInterp, vyInterp = getInterpolators(velocity.vx, "vxvy", np.floor(mx), np.floor(my), d2=velocity.vy)
    return np.array([t * (-vxInterp([y[0]], [y[1]], grid=False)), t * (-vyInterp([y[0]], [y[1]], grid=False))])


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
    Calculates the width of the ice stream at one point, (x1, y1).  (x0, y0) is there
    to give an idea of where the velocity width begins and ends which should be on a
    line which is perpindicular to the line from (x1, y1) to (x0, y0).
    :param x0:
    :param y0:
    :param x1:
    :param y1:
    :param draw:
    :return:
    '''

    # input is in map coordinates
    #
    # This is with interpolation
    #
    theta = np.arctan2(float(y1 - y0), float(x1-x0))
    #Rotation matrix:
    # rotMatrix =  np.matrix([[np.cos(theta), -1*np.sin(theta)],[np.sin(theta), np.cos(theta)]])
    # cos    -sin    x    =    x*cos + y*-sin    = y*-sin
    # sin     cos    y    =    x*sin + y*cos    = y*cos

    velInterp = getInterpolators(velocity.data, 'velocity', x1, y1)
    tx1, ty1 = projCoord(x1,y1)
    v0 = velInterp(tx1,ty1, grid=False)

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
            velInterp = getInterpolators(velocity.data, 'velocity', (x1 + (dr*dis * -np.sin(theta))), (y1 + (dr*dis * np.cos(theta))))
            tx, ty = projCoord(x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta)))  # Line perpindicular to flow
            currentVelocity = velInterp(tx, ty, grid=False)
            if np.abs(currentVelocity - vOld) > dv[i][2]:
                dv[i][0], dv[i][1] = x1 + (dr*dis * -np.sin(theta)), y1 + (dr*dis * np.cos(theta))
                dv[i][2] = np.abs(currentVelocity - vOld)
            if currentVelocity !=0:
                startEndRatio = v0/currentVelocity
            vOld = currentVelocity
        if currentVelocity < 5:
            # plotting line
            endPoints[i][0], endPoints[i][1] = dv[i][0], dv[i][1]#mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
        else:
            endPoints[i][0], endPoints[i][1] = mapCoord(tx, ty)#mapCoord(xa[ir], ya[ir])
    if draw:
        iiContainer.currentWidget().addItem(pg.PlotDataItem([endPoints[0][0], endPoints[1][0]], [endPoints[0][1], endPoints[1][1]], connect='all', pen=whitePlotPen))

        # circle plotting
        d = (0.5)*sqrt((endPoints[0][0]-endPoints[1][0])**2 + (endPoints[0][1]-endPoints[1][1])**2)
        cax, cay = endPoints[1][0] + (d * -np.sin(theta)), endPoints[1][1] + (d * np.cos(theta))
        xc, yc = circArr(cax, cay)
        iiContainer.currentWidget().addItem(pg.PlotDataItem(xc, yc, connect='all', pen=blackPlotPen))
        iiContainer.currentWidget().addItem(pg.PlotDataItem([cax], [cay], pen=blackPlotPen))
    return endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1]


def getInterpolators(d1, choice, x0, y0, x1=-99, y1=-99, d2=None):
    '''
    Determines the local interpolator and returns it.
    Interpolates data d1 (and d2 if necessary)
    If x1,y1 not specified then creates a local interpolator of size 10x10 with
    the point x0,y0 in the middle.

    :param d1:
    :param choice:  Which dataset to process
    :param x0:  IN MAP COORDINATES
    :param y0:  IN MAP COORDINATES
    :param x1:
    :param y1:
    :param d2:
    :return:
    '''
    # vel_x0 = -638000  # first x coordinate
    # vel_x1 = 864550  # last x coordinate
    # vel_y0 = -657600  # first y coordinate
    # vel_y1 = -3349350  # last y coordinate
    # vel_xarray = linspace(vel_x0, vel_x1, 10018, endpoint=True)
    # vel_yarray = linspace(vel_y1, vel_y0, 17946, endpoint=True)

    # Make sure points are in bounds
    p0 = [x0, y0]
    p1 = [x1, y1]
    p0[0], p0[1] = np.floor(p0[0]), np.floor(p0[1])

    if p0[0] < 0:
        p0[0] = 0
    if p0[0] > 10018:
        p0[0] = 10018
    if p0[1] < 0:
        p0[1] = 0
    if p0[1] > 17946:
        p0[1] = 17946
    if p1[0] == -99:
        # if the function is sent a single point then create interpolator around the point
        # in this case a 10x10 interpolator
        p1[0] = p0[0] - 5
        p1[1] = p0[1] - 5
        p0[0] = p0[0] + 5
        p0[1] = p0[1] + 5
    else:
        p1[0], p1[1] = np.floor(p1[0]), np.floor(p1[1])
        if p1[0] < 0:
            p1[0] = 0
        if p1[0] > 10018:
            p1[0] = 10018
        if p1[1] < 0:
            p1[1] = 0
        if p1[1] > 17946:
            p1[1] = 17946


    minSpacing = 10
    # p1, p2 in map coordinates
    projx0, projy0 = projCoord(p0[0], p0[1])
    projx1, projy1 = projCoord(p1[0], p1[1])
    #FIXME Should there be a minimum dx, dy?
    dx = 1 + math.fabs(p1[0] - p0[0])
    dy = 1 + math.fabs(p1[1] - p0[1])

    if p0[1] < p1[1]:
        vel_yarray = linspace(projy1, projy0, int(dy), endpoint=True)
        y0 = p0[1]
        y1 = p1[1] + 1
    else:
        vel_yarray = linspace(projy0, projy1, int(dy), endpoint=True)
        y0 = p1[1]
        y1 = p0[1] + 1
    if p0[0] < p1[0]:
        vel_xarray = linspace(projx0, projx1, int(dx), endpoint=True)
        x0 = p0[0]
        x1 = p1[0] + 1
    else:
        vel_xarray = linspace(projx1, projx0, int(dx), endpoint=True)
        x0 = p1[0]
        x1 = p0[0] + 1
    y0, y1, x0, x1 = int(y0), int(y1), int(x0), int(x1)
    possibleChoices = ['velocity', 'bed', 'surface', 'smb', 'thickness']
    if choice == 'vxvy':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose()), \
               RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d2[y0:y1, x0:x1])).transpose())
    elif choice in possibleChoices:
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    # elif choice is 'bed':
    #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    # elif choice is 'surface':
    #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    # elif choice is 'smb':
    #     #smb and all other RACMO data is 'upside down' compared to the rest of the data
    #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    # elif choice is 'v':
    #     #just velocity, not vx and vy
    #     return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    else:
        print "ERROR: No interpolator selected!  ./helper_files/data_functions.getInterpolators()"



    # RectBivariateSpline(vel_xarray, vel_yarray, (npk.flipud(vy)).transpose()) #FIXME SHOULD IT BE FLIPPED!?!?!?
    #FIXME Make sure flipping the transpose works!!!
    #FIXME changed the interpolator


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
        d += sqrt((vpts[i - 1].x - vpts[i].x) ** 2 + (vpts[i - 1].y - vpts[i].y) ** 2)
    print 'distance is: ', d*150, 'divide', int(d*(150/dr))

    for i in range(1, len(vpts)):
        '''
        This part compares neighbor points to each other. 
        '''

        theta = np.arctan2(float(vpts[i].getY() - vpts[i - 1].getY()), float(vpts[i].getX() - vpts[i - 1].getX()))
        # pvx0, pvy0 = projCoord(vpts[i-1].x, vpts[i-1].y)
        # pvx1, pvy1 = projCoord(vpts[i].x, vpts[i].y)
        # distance = (sqrt((vpts[i].getY() - vpts[i - 1].getY()) ** 2 + (vpts[i].getX() - vpts[i - 1].getX()) ** 2))
        distance = sqrt((vpts[i - 1].x - vpts[i].x) ** 2 + (vpts[i - 1].y - vpts[i].y) ** 2)
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
            # FIXME probably more elegant way to do this
            # transform coordinates into projected coordinates
            tx, ty = projCoord(vpts[i - 1].getX() + t[0, 0], vpts[i - 1].getY() + t[1, 0])
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
        if vWidthCheck == 2 or runModel:
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
        if smbCheck == 2 or runModel:
            # smbInterp = getInterpolators(smb.data, smb.name, mix, miy, x1=mxx, y1=mxy)
            localSMB = smb.interp(px, py, grid=False)
            smbValues.append(localSMB)

        ########################################
        ##   CALCULATE THICKNESS              ##
        ########################################
        # if smbCheck == 2 or runModel:
        # thickInterp   = getInterpolators(thickness.data, thickness.name, mix, miy, x1=mxx, y1=mxy)
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
            thickness.pathData = np.append(thickness.pathData, thickValues[i])

    smb.pathData = smb.pathData*(1.0/1000.0)*(916.7/1000.0) # millimeters -> meters then water-equivalent to ice-equivalent
    # print 'graphx[-1]: ', graphX[len(graphX) - 1]
    dist = 0
    for i in range(len(vpts) - 1):
        # print 'calc distance...'
        xd0, yd0 = projCoord(vpts[i].x, vpts[i].y)
        xd1, yd1 = projCoord(vpts[i + 1].x, vpts[i + 1].y)
        # print xd0, yd0, xd1, yd1
        dist += sqrt(((xd1 - xd0) ** 2 + (yd1 - yd0) ** 2))
    # print 'dist: ', dist
    thickness.distanceData = graphX
    if surfaceCheck.checkState() == 2 or runModel:
        surface.distanceData = graphX
    if bedCheck.checkState() == 2 or runModel:
        bed.distanceData = graphX
    if runModel:
        if velocityCheck.checkState() == 2 or runModel:
            velocity.distanceData = graphX
        if smbCheck.checkState() == 2 or runModel:
            smb.distanceData = graphX
        if vWidthCheck.checkState() == 2 or runModel:
            velocityWidth.distanceData = graphX
        if bedCheck.checkState() == 2 or runModel:
            bed.distanceData = graphX
        if surfaceCheck.checkState() == 2 or runModel:
            surface.distanceData = graphX
        # return linePoints, np.array(graphX)
        # else:
        #     return nbed, nsurf  # , linePoints
    print 'SHAPES'
    print len(thickness.distanceData)
    print len(bed.distanceData)



# def interpolateData_eh(runModel):
#     print 'Starting interpolateData: '
#     if velocityCheck.checkState() == 2 or runModel:
#         # vxInterp, vyInterp = getInterpolators(velocity.vx, velocity.name, mix, miy, x1=mxx, y1=mxy, d2=velocity.vy)
#         _interpolateData(velocity, './data/velInterp.pkl')
#     if surfaceCheck.checkState() == 2 or runModel:
#         _interpolateData(surface, './data/surfaceInterp.pkl')
#     if bedCheck.checkState() == 2 or runModel:
#         _interpolateData(bed, './data/bedInterp.pkl')
#     if smbCheck.checkState() == 2 or runModel:
#         _interpolateData(smb, './data/smbInterp.pkl')
#     _interpolateData(thickness, './data/thicknessInterp.pkl')
#
# def _interpolateData(ds, interpName):
#     '''
#     Calculate the data for bottom plot or to run the model.
#     If botPlotBool, calculate all the data.  Else, calculate just bed/surface.
#     :return:
#     '''
#     start_time = time.time()
#     global dr, bpLegend, dataLen, botPlot
#     botPlot = True
#     arrayValues = []
#     # xValues = []
#     # smbValues = []
#     # surfValues = []
#     # bedValues = []
#     # thickValues = []
#     linePoints = [0]
#     # vwValues = []
#     graphX = []
#
#     interpPickle = open(interpName, 'r')
#     interp = pickle.load(interpPickle)
#     interpPickle.close()
#     print 'Loaded interp: ', time.time() - start_time
#
#     # thickInterp = getInterpolators(thickness.data, thickness.name, mix, miy, x1=mxx, y1=mxy)
#
#
#     # Find a distance ~150m which gets close to dividing the distance between first 2 spots
#
#     for i in range(1, len(vpts)):
#         '''
#         This part compares neighbor points to each other.
#         '''
#
#         theta = np.arctan2(float(vpts[i].getY() - vpts[i - 1].getY()), float(vpts[i].getX() - vpts[i - 1].getX()))
#         distance = sqrt((vpts[i - 1].x - vpts[i].x) ** 2 + (vpts[i - 1].y - vpts[i].y) ** 2)
#         # remainder = distance % dr
#         # Xline needs to be in map coordinates because tx,ty are in map coordinates
#         xline = linspace(0, distance, distance, endpoint=True)  # * const makes it every 150/const meters
#         '''
#         #FIXME NEED TO CHANGE SO IT LINES UP WITH FENICS MESH
#         '''
#
#         linePoints.append(distance * (1 / dr) + linePoints[-1])
#
#         # Rotation matrix:
#         rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])
#         px = []  # px, py are projected coordinates used to get values from the interpolator.  Projected meaning IRL values
#         py = []
#
#         for j in range(len(xline)):
#             # rotate coordinates
#             t = rotMatrix * np.matrix([[xline[j]], [0.0]])
#             # FIXME probably more elegant way to do this
#             # transform coordinates into projected coordinates
#             tx, ty = projCoord(vpts[i - 1].getX() + t[0, 0], vpts[i - 1].getY() + t[1, 0])
#             px.append(tx)
#             py.append(ty)
#             if len(px) > 1:
#                 graphX.append(graphX[len(graphX) - 1] + sqrt((px[-1] - px[-2]) ** 2 + (py[-1] - py[-2]) ** 2))
#                 # print 'dist: ', graphX[-1] + sqrt((px[-1]-px[-2])**2 + (py[-1]-py[-2])**2)
#             elif len(graphX) == 0:
#                 graphX.append(0)
#             else:
#                 graphX.append(graphX[len(graphX) - 1])  # + sqrt((px[-1]) ** 2 + (py[-1]) ** 2))
#                 #     print 'dis3: ', graphX[-1] + sqrt((px[-1]) ** 2 + (py[-1]) ** 2)
#
#         ########################################
#         ##    CALCULATE VALUES                ##
#         ########################################
#         localValues = interp(px, py, grid=False)
#         arrayValues.append(localValues)
#
#         ########################################
#         ##   COMPILE DATA                     ##
#         ########################################
#         ds.pathData = np.array(arrayValues)
#
#         for i in range(1, len(arrayValues)):
#             ds.pathData = np.append(ds.pathData, arrayValues[i])
#
#     if ds.name == 'smb':
#         smb.pathData = smb.pathData * (1.0 / 1000.0) * (916.7 / 1000.0)  # millimeters -> meters then water-equivalent to ice-equivalent
#     # print 'graphx[-1]: ', graphX[len(graphX) - 1]
#     dist = 0
#     for i in range(len(vpts) - 1):
#         # print 'calc distance...'
#         xd0, yd0 = projCoord(vpts[i].x, vpts[i].y)
#         xd1, yd1 = projCoord(vpts[i + 1].x, vpts[i + 1].y)
#         # print xd0, yd0, xd1, yd1
#         dist += sqrt(((xd1 - xd0) ** 2 + (yd1 - yd0) ** 2))
#     # print 'dist: ', dist
#     ds.distanceData = graphX
#     print 'Finished one dataset: ', time.time() - start_time
#
