import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import distmesh
import numpy as np
import matplotlib.pyplot as plt
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
from dolfin.cpp.mesh import *
import fenics as fc
import distmesh as dm
from dolfin.cpp.io import File

# Local imports
from ..math_functions import *
from ..pens import *
from ..gui import *
from ..dataset_objects import *
from ..data_functions import *
from PlotPoint import *
from ..velocity_functions import *
from StaticPlotter import *

class MeshGUI(QtGui.QMainWindow):
    def __init__(self, parent):
        self.parent = parent
        QtGui.QMainWindow.__init__(self, self.parent)
        self.cw = QtGui.QWidget()
        self.setCentralWidget(self.cw)
        self.horLayout = QtGui.QHBoxLayout()
        self.cw.setLayout(self.horLayout)
        self.gui = QtGui.QGridLayout()
        self.pw = pg.PlotWidget()
        self.pw.invertY(True)
        self.pw.setAspectLocked(True)
        self.horLayout.addWidget(self.pw)
        self.horLayout.addWidget(self.gui)

        self.fileLabel = QtGui.QLabel('Save File:')
        self.fileName = QtGui.QLineEdit()

        self.meshLabel = QtGui.QLabel('Save Mesh:')
        self.meshName = QtGui.QLineEdit()

        self.gui.addWidget(self.meshLabel, 0, 0)
        self.gui.addWidget(self.meshName, 0, 1)
        self.gui.addWidget(self.fileLabel, 1, 0)
        self.gui.addWidget(self.fileName, 1, 1)

        self.show()

        self.p, self.t = self.runPoly()
        self.mesh = None


    def graphMesh(self):
        x, y, c = [], [], []
        for row in self.t:
            x.append(self.p[row[0]][0])
            y.append(self.p[row[0]][1])
            c.append(1)
            x.append(self.p[row[1]][0])
            y.append(self.p[row[1]][1])
            c.append(1)
            x.append(self.p[row[2]][0])
            y.append(self.p[row[2]][1])
            c.append(1)
            x.append(self.p[row[0]][0])
            y.append(self.p[row[0]][1])
            c.append(0)
        self.pw.getPlotItem().plot(x, y, pen=(255, 0, 0), connect=np.array(c))


    def polygon(self, pv, d, minx, miny, maxx, maxy):  # GOES COUNTER_CLOCKWISE
        """Polygon"""
        fd = lambda p: dm.dpoly(p, pv)
        print 'bleh'
        return dm.distmesh2d(fd, dm.huniform, d / 5, (minx, miny, maxx, maxy), pv)

    def fstats(self, p, t):
        print('%d nodes, %d elements, min quality %.2f'
              % (len(p), len(t), dm.simpqual(p, t).min()))

    def runPoly(self):
        px, py = 1, 1  # map['cmap_x1'], map['cmap_y1']
        pv = [[vpts[0].px, vpts[0].py]]
        dmin = 99999
        minx = 99999
        miny = 99999
        maxx = -1
        maxy = -1

        for i in range(1, len(vpts)):
            d = sqrt((vpts[i].px - vpts[i - 1].px) ** 2 + (vpts[i].py - vpts[i - 1].py) ** 2)
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
        pause = lambda: None
        # plt.ion()
        np.random.seed(1)  # Always the same results
        p, t = self.polygon(np.array(pv), dmin, minx, miny, maxx, maxy)
        self.fstats(p, t)
        # pause()
        return p, t

    def saveMeshAsXML(self, p, t):
        print 'Writing mesh to ' + str(self.meshName.text())
        f = open(self.meshName.text(), 'w')
        f.write('<dolfin>\n')
        f.write('\t<mesh celltype=\"triangle\" dim=\"2\">\n')
        f.write('\t\t<vertices size=\"\t' + str(len(p)) + '\">\n')
        for i in range(len(p)):
            f.write('\t\t\t<vertex index=\"\t' + str(i) + '\" x=\"\t' + str(p[i][0]) + '\t\" y=\"\t' + str(
                p[i][1]) + '\t\" z=\"0\"/>\n')
        f.write('\t\t</vertices>\n')
        f.write('\t\t<cells size=\"\t' + str(len(t)) + '\">\n')
        for i in range(len(t)):
            f.write('\t\t\t<triangle index=\"\t' + str(i) + '\" v0=\"\t' + str(t[i][0]) + '\" v1=\"\t' + str(
                t[i][1]) + '\" v2=\"\t' + str(t[i][2]) + '\"/>\n')
        f.write('\t\t</cells>\n')
        f.write('\t</mesh>\n')
        f.write('</dolfin>\n')
        f.close()
        print 'Mesh file done.'

    def writeToHDF5(self):
        self.mesh = Mesh(self.meshName.text())
        hfile = fc.HDF5File(self.mesh.mpi_comm(), self.fileName.text(), "w")
        V = fc.FunctionSpace(self.mesh, 'CG', 1)

        thicknessModelData = thickness.interp(self.mesh.coordinates()[::, 0], self.mesh.coordinates()[::, 1], grid=False)
        bedModelData = bed.interp(self.mesh.coordinates()[::, 0], self.mesh.coordinates()[::, 1], grid=False)
        surfaceModelData = surface.interp(self.mesh.coordinates()[::, 0], self.mesh.coordinates()[::, 1], grid=False)
        smbModelData = smb.interp(self.mesh.coordinates()[::, 0], self.mesh.coordinates()[::, 1], grid=False)
        velocityModelData = velocity.interp(self.mesh.coordinates()[::, 0], self.mesh.coordinates()[::, 1], grid=False)

        functThickness = fc.Function(V, name="Thickness")
        functBed = fc.Function(V, name="Bed")
        functSurface = fc.Function(V, name="Surface")
        functSMB = fc.Function(V, name='SMB')
        functVelocity = fc.Function(V, name='Velocity')

        print 'len: ', len(functThickness.vector()[:])
        print 'len: ', len(thicknessModelData)

        functThickness.vector()[:] = thicknessModelData
        functBed.vector()[:] = bedModelData
        functSurface.vector()[:] = surfaceModelData
        functSMB.vector()[:] = smbModelData
        functVelocity.vector()[:] = velocityModelData

        hfile.write(functThickness, "thickness")
        hfile.write(functBed, "bed")
        hfile.write(functSurface, "surface")
        hfile.write(functSMB, "smb")
        hfile.write(functVelocity, "velocity")
        hfile.write(self.mesh, "mesh")
        hfile.close()
        print 'mesh ', self.mesh.coordinates()[::, 0], self.mesh.coordinates()[::, 1]
        print 'velocity', velocityModelData
        print 'thick', thicknessModelData
        print 'bed', bedModelData
        print 'surfae', surfaceModelData
        print 'smb', smbModelData

        paraF = File('paraf.pvd')
        # paraF << self.mesh
        paraF << functBed
        paraF << functSurface
        paraF << functSMB
        paraF << functVelocity
        paraF << functThickness

    def meshGui(self):
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
        p, t = self.runPoly()
        self.saveMeshAsXML(p, t, './data/mesh2d.xml')
        self.writeToHDF5(p, t, './data/mesh2d.h5', './data/mesh2d.xml')

        # x, y, c = [], [], []
        # for row in t:
        #     x.append(p[row[0]][0])
        #     y.append(p[row[0]][1])
        #     c.append(1)
        #     x.append(p[row[1]][0])
        #     y.append(p[row[1]][1])
        #     c.append(1)
        #     x.append(p[row[2]][0])
        #     y.append(p[row[2]][1])
        #     c.append(1)
        #     x.append(p[row[0]][0])
        #     y.append(p[row[0]][1])
        #     c.append(0)
        # meshPW.getPlotItem().plot(x, y, pen=(255, 0, 0), connect=np.array(c))





