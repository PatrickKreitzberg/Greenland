import numpy as np
from pylab import sqrt,linspace,array,argmax
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline
import math
from constants import *

'''
colorToProj:  color -> projected
dataToProj:   data -> projected
colorCoord:   data, projected -> color
dataToColor:  data -> color
colorToData:  color -> data
dataCoord:    color, projected -> data
'''

def colorToProj(x,y):
    # returns x,y in the global projected coordinates
    # no way to check if in color or map coord
    if y >= map['cmap_proj_y1'] and y <= map['cmap_proj_y0']:
        if x >= map['cmap_proj_x0'] and x <= map['cmap_proj_x1']:
            return x, y
        else:
            print 'ERROR IN PROJCOORD'
            return -1
    else:
        return ((150*x) + float(map['cmap_proj_x0'])), ((-150*y) + float(map['cmap_proj_y0']))

def dataToProj(x,y):
    return dataToColor(colorToProj(x,y))

def colorCoord(x,y):
    # must assume data in data coord or proj coord
    # first return data to color
    # else return proj to map
    if y >= map['y0'] and y <=map['y1']:
        if x >= map['x0'] and x <= map['x1']:
            return dataToColor(x,y)
        else:
            print 'ERROR IN MAPCOORD'
            return -1
    else:
        return ((-float(map['cmap_proj_x0']) + x)/150.0), (-(-float(map['cmap_proj_y0']) + y)/150.0)

def dataToColor(x,y):
    # turns colormap data point into data point
    return x*(float(map['cmap_x1'])/float(map['x1'])), y*(float(map['cmap_y1'])/float(map['y1']))

def colorToData(x,y):
    # turns colormap data point into data point
    return x*(float(map['x1'])/float(map['cmap_x1'])), y*(float(map['y1'])/float(map['cmap_y1']))

def dataCoord(x,y):
    #m ust assume data is in colormap or projected coord
    if y >= map['cmap_y0'] and y <=map['cmap_y1']:
        if x >= map['cmap_x0'] and x <= map['cmap_x1']:
           return colorToData(x,y)
    elif y >= map['cmap_proj_y0'] and y <=map['cmap_proj_y1']:
        if x >= map['cmap_proj_x0'] and x <= map['cmap_proj_x1']:
            return colorToData(colorCoord(x,y))
    else:
        print 'ERROR: dataCoord(x, y) error'
        return -1

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
    snapDistance = 2
    minD = 31
    i = 0
    while i < len(cData[0]): # and not found:
        d = sqrt((x0-cData[0][i])**2 + (y0-cData[1][i])**2)
        if not found and d < snapDistance:
            found = True
            minD = d
            imin = i
        elif found and d < minD:
            minD = d
            imin = i
        i += 1
    return imin




