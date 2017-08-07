import distmesh
from math_functions import *
from pens import *
from gui import *
from dataset_objects import *
import numpy as np
from data_functions import *
from scipy.integrate import ode
from classes.PlotPoint import *
from velocity_functions import *
from classes.StaticPlotter import *
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg

# Local imports.
import distmesh as dm

# Polygon example:
# pv are the vertices
# @dpoly gives the distance function
# @huniform Implements the trivial uniform mesh size function h=1.

def polygon(pv,d, minx,miny, maxx,maxy):  # GOES COUNTER_CLOCKWISE
    """Polygon"""
    # pv0 = np.array([(-0.4,-0.5),(0.4,-0.2),(0.4,-0.7),(1.5,-0.4),(0.9,0.1),
    #                (1.6,0.8),(0.5,0.5),(0.2,1.0),(0.1,0.4),(-0.7,0.7),
    #                (-0.4,-0.5)])
    # print type(pv0)
    # print type(pv0[0])
    fd = lambda p: dm.dpoly(p, pv)
    print 'bleh'
    return dm.distmesh2d(fd, dm.huniform, d/5, (minx,miny, maxx, maxy), pv)

def fstats(p, t):
    print('%d nodes, %d elements, min quality %.2f'
          % (len(p), len(t), dm.simpqual(p,t).min()))


def runPoly():
    cx, cy = 1,1 #map['cmap_x1'], map['cmap_y1']
    pv = [[vpts[0].cx/cx, vpts[0].cy/cy]]
    dmin = 99999
    minx = 99999
    miny = 99999
    maxx = -1
    maxy = -1

    for i in range(1,len(vpts)):
        d = sqrt((vpts[i].cx/cx - vpts[i-1].cx/cx)**2 + (vpts[i].cy/cy - vpts[i-1].cy/cy)**2)
        if d < dmin:
            dmin = d
        if vpts[i].cx < minx:
            minx = vpts[i].cx

        if vpts[i].cx > maxx:
            maxx = vpts[i].cx

        if vpts[i].cy < miny:
            miny = vpts[i].cy

        if vpts[i].cy > maxy:
            maxy = vpts[i].cy

        pv.append([vpts[i].cx/cx, vpts[i].cy/cy])
    print dmin
    print pv
    pause = lambda : None
    # plt.ion()
    np.random.seed(1) # Always the same results
    p, t = polygon(np.array(pv),dmin, minx/cx, miny/cy, maxx/cx, maxy/cy)
    fstats(p, t)
    # pause()
    print('meshout')
    print 'p: The node positions len: ', len(p)
    print p
    print 't: The triangle indices len: ', len(t)
    print t
    return p, t

def meshGui():
    meshWindow = QtGui.QMainWindow(mw)
    cw = QtGui.QWidget()
    meshWindow.setCentralWidget(cw)
    cwLayout = QtGui.QHBoxLayout()
    cw.setLayout(cwLayout)
    meshPW = pg.PlotWidget()
    meshPW.invertY(True)
    meshPW.setAspectLocked(True)
    cwLayout.addWidget(meshPW)
    meshWindow.show()
    p, t = runPoly()
    x, y, c = [], [], []
    for row in t:
        x.append(p[row[0]][0])
        y.append(p[row[0]][1])
        c.append(1)
        x.append(p[row[1]][0])
        y.append(p[row[1]][1])
        c.append(1)
        x.append(p[row[2]][0])
        y.append(p[row[2]][1])
        c.append(1)
        x.append(p[row[0]][0])
        y.append(p[row[0]][1])
        c.append(0)

    meshPW.getPlotItem().plot(x, y, pen=(255, 0, 0), connect=np.array(c))


'''

An essential decision is how to represent the geometry 
(the shape of the region). Our code uses a signed distance 
function d(x, y), negative inside the region. 

*let user pick length*
The edge lengths should be close to the relative size h(x) specified by
the user (the lengths are nearly equal when the user chooses h(x) = 1)




'''