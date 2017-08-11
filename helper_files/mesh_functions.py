import distmesh
from math_functions import *
from pens import *
from gui import *
from dataset_objects import *
from data_functions import *
from classes.PlotPoint import *
from velocity_functions import *
from classes.StaticPlotter import *
import numpy as np
import matplotlib.pyplot as plt
from constants import *
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
from dolfin.cpp.mesh import *
import fenics as fc
# Local imports.
import distmesh as dm
from dolfin.cpp.io import File

# Polygon example:
# pv are the vertices
# OUT:
# p is the node positions   (nx2)
# t is the triangle indices (nx3)
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
    px, py = 1,1 #map['cmap_x1'], map['cmap_y1']
    pv = [[vpts[0].px, vpts[0].py]]
    dmin = 99999
    minx = 99999
    miny = 99999
    maxx = -1
    maxy = -1

    for i in range(1,len(vpts)):
        d = sqrt((vpts[i].px - vpts[i-1].px)**2 + (vpts[i].py - vpts[i-1].py)**2)
        if d < dmin:
            dmin = d
        if vpts[i].px < minx:
            minx = vpts[i].px

        if vpts[i].px > maxx:
            maxx = vpts[i].px

        if vpts[i].py < miny:
            miny = vpts[i].py

        if vpts[i].py > maxy:
            maxy = vpts[i].py

        pv.append([vpts[i].px, vpts[i].py])
    pause = lambda : None
    # plt.ion()
    np.random.seed(1) # Always the same results
    p, t = polygon(np.array(pv),dmin, minx, miny, maxx, maxy)
    fstats(p, t)
    # pause()
    return p, t

def saveMeshAsXML(p, t, fname):
    print 'Writing mesh to ' + str(fname)
    f = open(fname, 'w')
    f.write('<dolfin>\n')
    f.write('\t<mesh celltype=\"triangle\" dim=\"2\">\n')
    f.write('\t\t<vertices size=\"\t'+ str(len(p)) +'\">\n')
    for i in range(len(p)):
        f.write('\t\t\t<vertex index=\"\t' + str(i) +'\" x=\"\t' + str(p[i][0]) +'\t\" y=\"\t' + str(p[i][1]) +'\t\" z=\"0\"/>\n')
    f.write('\t\t</vertices>\n')
    f.write('\t\t<cells size=\"\t' + str(len(t)) + '\">\n')
    for i in range(len(t)):
        f.write('\t\t\t<triangle index=\"\t' + str(i) +'\" v0=\"\t' + str(t[i][0]) +'\" v1=\"\t' + str(t[i][1]) +'\" v2=\"\t' + str(t[i][2]) +'\"/>\n')
    f.write('\t\t</cells>\n')
    f.write('\t</mesh>\n')
    f.write('</dolfin>\n')
    f.close( )
    print 'Mesh file done.'

def writeToHDF5(p, t, fname, meshname):
    mesh = Mesh(meshname)
    hfile = fc.HDF5File(mesh.mpi_comm(), fname, "w")
    V = fc.FunctionSpace(mesh, 'CG', 1)

    thicknessModelData = thickness.interp(mesh.coordinates()[::, 0], mesh.coordinates()[::, 1], grid=False)
    bedModelData       = bed.interp(mesh.coordinates()[::, 0], mesh.coordinates()[::, 1], grid=False)
    surfaceModelData   = surface.interp(mesh.coordinates()[::, 0], mesh.coordinates()[::, 1], grid=False)
    smbModelData       = smb.interp(mesh.coordinates()[::, 0], mesh.coordinates()[::, 1], grid=False)
    velocityModelData  = velocity.interp(mesh.coordinates()[::, 0], mesh.coordinates()[::, 1], grid=False)

    functThickness = fc.Function(V, name="Thickness")
    functBed       = fc.Function(V, name="Bed")
    functSurface   = fc.Function(V, name="Surface")
    functSMB       = fc.Function(V, name='SMB')
    functVelocity  = fc.Function(V, name='Velocity')

    print 'len: ', len(functThickness.vector()[:])
    print 'len: ', len(thicknessModelData)

    functThickness.vector()[:] = thicknessModelData
    functBed.vector()[:]       = bedModelData
    functSurface.vector()[:]   = surfaceModelData
    functSMB.vector()[:]       = smbModelData
    functVelocity.vector()[:]  = velocityModelData

    hfile.write(functThickness, "thickness")
    hfile.write(functBed,       "bed")
    hfile.write(functSurface,   "surface")
    hfile.write(functSMB,       "smb")
    hfile.write(functVelocity,  "velocity")
    hfile.write(mesh,           "mesh")
    hfile.close()
    print 'mesh ', mesh.coordinates()[::, 0], mesh.coordinates()[::, 1]
    print 'velocity', velocityModelData
    print 'thick', thicknessModelData
    print 'bed', bedModelData
    print 'surfae', surfaceModelData
    print 'smb', smbModelData


    paraF = File('paraf.pvd')
    # paraF << mesh
    paraF << functBed
    paraF << functSurface
    paraF << functSMB
    paraF << functVelocity
    paraF << functThickness


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
    saveMeshAsXML(p,t,'./data/mesh2d.xml')
    writeToHDF5(p, t, './data/mesh2d.h5', './data/mesh2d.xml')
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