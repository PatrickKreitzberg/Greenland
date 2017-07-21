import numpy as np
from pylab import sqrt,linspace,array,argmax
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
import math

def projCoord(x,y):
    # returns x,y in the global projected coordinates
    if y >= -3349350 and y <= -657600:
        if x >= -638000 and x <= 864550:
            return x, y
        else:
            print 'ERROR IN PROJCOORD'
            return -1
    else:
        return ((150*x) - 637925), ((-150*y) - 657675)

def mapCoord(x,y):
    # returns x,y in map coordinates which is x <- (0,1018) y<- (0,1746) roughly
    if y >= 0 and y <=17946:
        if x <= 10018:
            return x,y
        else:
            print 'ERROR IN MAPCOORD'
            return -1
    else:
        return ((637925 + x)/150.0), (-(657675 + y)/150.0)


def findSlopes(lines, vlist):

    for i in range(len(vlist)-1):
        m = float(vlist[i+1][1] -vlist[i][1])/float(vlist[i+1][0] -vlist[i][0])
        b = float(vlist[i][1]) - m*(float(vlist[i][0]))
        lines.append([m,b])
    #for last point
    j = len(vlist)-1
    m = float(vlist[j][1] -vlist[j-1][1])/float(vlist[j][0] -vlist[j-1][0])
    b = float(vlist[j][1]) - m*(float(vlist[j][0]))
    lines.append([m, b])
    return lines

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
    print 'vel_xarray, vel_yarray ', vel_xarray, vel_yarray
    print '              d1.shape ',d1[y0:y1, x0:x1].shape
    if choice is 'velocity':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose()), RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d2[y0:y1, x0:x1])).transpose())
    elif choice is 'bed':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    elif choice is 'surface':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    elif choice is 'smb':
        #smb and all other RACMO data is 'upside down' compared to the rest of the data
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())
    elif choice is 'v':
        #just velocity, not vx and vy
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1[y0:y1, x0:x1])).transpose())



    # RectBivariateSpline(vel_xarray, vel_yarray, (npk.flipud(vy)).transpose()) #FIXME SHOULD IT BE FLIPPED!?!?!?
    #FIXME Make sure flipping the transpose works!!!
    #FIXME changed the interpolator

def circArr(x,y):
    r = 1
    t = linspace(0,2*np.pi,50)
    xArr = x + r*np.cos(t)
    yArr = y + r*np.sin(t)
    return xArr, yArr


def curveDistance(x0, y0, cData):
    imin = -1
    found = False
    snapDistance = 30
    minD = 31
    for i in range(len(cData[0])):
        d = sqrt((x0-cData[0][i])**2 + (y0-cData[1][i])**2)
        if not found and d < snapDistance:
            found = True
            minD = d
            imin = i
        elif found and d < minD:
            minD = d
            imin = i
    return imin




