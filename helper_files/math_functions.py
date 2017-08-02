import numpy as np
from pylab import sqrt,linspace,array,argmax
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
import math
from constants import *

def projCoord(x,y):
    # returns x,y in the global projected coordinates
    if y >= proj_y1 and y <= proj_y0:
        if x >= proj_x0 and x <= proj_x1:
            return x, y
        else:
            print 'ERROR IN PROJCOORD'
            return -1
    else:
        return ((150*x) + proj_x0), ((-150*y) + proj_y0)

def mapCoord(x,y):
    # returns x,y in map coordinates which is x <- (0,1018) y<- (0,1746) roughly
    if y >= map['y0'] and y <=map['y1']:
        if x >= map['x0'] and x <= map['x1']:
            return x,y
        else:
            print 'ERROR IN MAPCOORD'
            return -1
    else:
        return ((-proj_x0 + x)/150.0), (-(-proj_y0 + y)/150.0)

def colorToData(x,y):
    # turns colormap data point into data point
    return x*(map['x1']/map['cmap_x1']), x*(map['x1']/10018)

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




