import numpy as np
from pylab import sqrt,linspace,array,argmax
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
import math


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

def getInterpolators(d1, choice, d2=None):
    vel_x0 = -638000  # first x coordinate
    vel_x1 = 864550  # last x coordinate
    vel_y0 = -657600  # first y coordinate
    vel_y1 = -3349350  # last y coordinate
    vel_xarray = linspace(vel_x0, vel_x1, 10018, endpoint=True)
    vel_yarray = linspace(vel_y1, vel_y0, 17946, endpoint=True)
    if choice is 'velocity':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1)).transpose()), RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d2)).transpose()),
    elif choice is 'bed':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1)).transpose())
    elif choice is 'surface':
        return RectBivariateSpline(vel_xarray, vel_yarray, (np.flipud(d1)).transpose())
    elif choice is 'smb':
        return RectBivariateSpline(vel_xarray, vel_yarray, d1.transpose())


    # RectBivariateSpline(vel_xarray, vel_yarray, (npk.flipud(vy)).transpose()) #FIXME SHOULD IT BE FLIPPED!?!?!?
    #FIXME Make sure flipping the transpose works!!!
    #FIXME changed the interpolator

def circArr(x,y):
    r = 1
    t = linspace(0,2*np.pi,50)
    xArr = x + r*np.cos(t)
    yArr = y + r*np.sin(t)
    return xArr, yArr

def projCoord(x,y):
    return ((150*x) - 637925), ((-150*y) - 657675)

def mapCoord(x,y):
    #returns x,y in map coordinates which is x <- (0,1018) y<- (0,1746) roughly
    return ((637925 + x)/150.0), (-(657675 + y)/150.0)

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




