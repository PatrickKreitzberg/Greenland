import distmesh
import numpy as np
import matplotlib.pyplot as plt
from pyqtgraph.Qt import QtGui
import pyqtgraph as pg
from dolfin.cpp.mesh import *
import fenics as fc
import distmesh as dm
from dolfin.cpp.io import File
import os

# Local imports
from ..math_functions import *
from ..pens import *
from ..gui import *
from ..dataset_objects import *
from ..data_functions import *
from PlotPoint import *
from ..constants import *

from InterpolateData import interpolateDataClass
from ..velocity_functions import *
from StaticPlotter import *

class MeshGUI(QtGui.QMainWindow):
    def __init__(self, parent):
        self.parent = parent
        QtGui.QMainWindow.__init__(self, self.parent)
        # basic gui setup
        self.cw = QtGui.QWidget()
        self.setCentralWidget(self.cw)
        self.horLayout = QtGui.QHBoxLayout()
        self.cw.setLayout(self.horLayout)
        self.guiWidget = QtGui.QWidget()
        self.gui = QtGui.QGridLayout()
        self.guiWidget.setLayout(self.gui)
        self.pw = pg.PlotWidget()
        self.pw.invertY(True)
        self.pw.setAspectLocked(True)
        self.horLayout.addWidget(self.pw)
        self.horLayout.addWidget(self.guiWidget)

        # right side buttons/textfields


        mn = 999999
        mx = 0
        for i in range(len(vpts)):
            for j in range(len(vpts)):
                if i != j:
                    d = sqrt((vpts[i].cx - vpts[j].cx)**2 + (vpts[i].cy - vpts[j].cy)**2)
                    if d < mn:
                        mn = d
                    if d > mx:
                        mx = d
        minTxt = 'Min side distance: ' + "{:.3f}".format(mn * 150) + 'm'
        maxTxt = 'Max side distance: ' + "{:.3f}".format(mx * 150) + 'm'
        self.lenMinLabel = QtGui.QLabel(minTxt)
        self.lenMaxLabel = QtGui.QLabel(maxTxt)
        self.distLabel   = QtGui.QLabel('Mesh len(m):')
        self.distInput   = QtGui.QLineEdit()
        self.genMeshButt = QtGui.QPushButton('Generate Mesh')
        self.genMeshButt.setEnabled(False)

        self.fileNameLabel   = QtGui.QLabel('Data Dir: ./data/')
        self.fileName = QtGui.QLineEdit()

        self.saveMeshButt = QtGui.QPushButton('Save Mesh')

        mnmxlabelrow = 0
        saverow = mnmxlabelrow + 4

        self.gui.addWidget(self.lenMinLabel, mnmxlabelrow,     0, 1, 2)
        self.gui.addWidget(self.lenMaxLabel, mnmxlabelrow + 1, 0, 1, 2) # row 2
        self.gui.addWidget(self.distLabel,   mnmxlabelrow + 2, 0)
        self.gui.addWidget(self.distInput,   mnmxlabelrow + 2, 1)
        self.gui.addWidget(self.fileNameLabel, saverow+1, 0)
        self.gui.addWidget(self.fileName, saverow+1, 1)
        self.gui.addWidget(self.genMeshButt, saverow + 2, 0, 1, 2)
        self.show()
        self.p, self.t, self.meshLen, self.savedMeshName = None, None, -1, None
        self.distInput.textChanged.connect(self.setMeshLength)
        self.genMeshButt.clicked.connect(self.generateMesh)
        # self.p, self.t = self.runPoly()
        # self.mesh = None
        # self.graphMesh()

    def generateMesh(self):
        if not os.path.exists('./data/' + str(self.fileName.text())):
            os.makedirs('./data/' + str(self.fileName.text()))
        self.hdfFileName = './data/' + str(self.fileName.text()) + '/' + str(self.fileName.text()) + '_mesh.h5'
        self.xmlFileName = './data/' + str(self.fileName.text()) + '/' + str(self.fileName.text()) + '_mesh_data.xml'
        self.paraFileName = './data/' + str(self.fileName.text()) + '/' + str(self.fileName.text()) + '_para.pvd'
        self.p, self.t = self.runPoly()
        self.graphMesh()
        self.saveMeshAsXML()
        self.writeToHDF5()

    def setMeshLength(self):
        if self.distInput is not '':
            self.genMeshButt.setEnabled(True)
            self.genMeshButt.keyboardGrabber()

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


    def polygon(self, pv, h0, minx, miny, maxx, maxy):  # GOES COUNTER_CLOCKWISE
        """Polygon"""
        fd = lambda p: dm.dpoly(p, pv)
        print 'Working on generating the mesh!'
        return dm.distmesh2d(fd, dm.huniform, h0, (minx, miny, maxx, maxy), pv)

    def fstats(self, p, t):
        print('%d nodes, %d elements, min quality %.2f'
              % (len(p), len(t), dm.simpqual(p, t).min()))

    def runPoly(self):
        try:
            h0 = float(self.distInput.text())
            print 'h0 is ', h0
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
            p, t = self.polygon(np.array(pv), h0, minx, miny, maxx, maxy)
            self.fstats(p, t)
            # pause()
            return p, t
        except ValueError:
            print 'Must enter a float for desired mesh length'
            return -1, -1

    def saveMeshAsXML(self):
        print 'Writing mesh'
        f = open(self.xmlFileName, 'w')
        f.write('<dolfin>\n')
        f.write('\t<mesh celltype=\"triangle\" dim=\"2\">\n')
        f.write('\t\t<vertices size=\"\t' + str(len(self.p)) + '\">\n')
        for i in range(len(self.p)):
            f.write('\t\t\t<vertex index=\"\t' + str(i) + '\" x=\"\t' + str(self.p[i][0]) + '\t\" y=\"\t' + str(
                self.p[i][1]) + '\t\" z=\"0\"/>\n')
        f.write('\t\t</vertices>\n')
        f.write('\t\t<cells size=\"\t' + str(len(self.t)) + '\">\n')
        for i in range(len(self.t)):
            f.write('\t\t\t<triangle index=\"\t' + str(i) + '\" v0=\"\t' + str(self.t[i][0]) + '\" v1=\"\t' + str(
                self.t[i][1]) + '\" v2=\"\t' + str(self.t[i][2]) + '\"/>\n')
        f.write('\t\t</cells>\n')
        f.write('\t</mesh>\n')
        f.write('</dolfin>\n')
        f.close()
        print 'Mesh file done, you may now save the data.'
        self.writeToHDF5()


    def writeToHDF5(self):
        print 'self.hdfFileName', self.hdfFileName, type(self.hdfFileName)
        mesh = Mesh(self.xmlFileName)
        hfile = fc.HDF5File(mesh.mpi_comm(), self.hdfFileName, "w")

        V = fc.FunctionSpace(mesh, 'CG', 1)

        thicknessiD = interpolateDataClass(thickness.interp, degree=2)
        bediD       = interpolateDataClass(bed.interp, degree=2)
        surfaceiD   = interpolateDataClass(surface.interp, degree=2)
        smbiD       = interpolateDataClass(smb.interp, degree=2)
        velocityiD  = interpolateDataClass(velocity.interp, degree=2)

        th = project(thicknessiD, V)
        be = project(bediD, V)
        su = project(surfaceiD, V)
        sm = project(smbiD, V)
        ve = project(velocityiD, V)

        hfile.write(th, 'thickness')
        hfile.write(be, 'bed')
        hfile.write(su, 'surface')
        hfile.write(sm, 'smb')
        hfile.write(ve, 'velocity')
        hfile.write(mesh, "mesh")
        hfile.close()

        paraF = File(self.paraFileName)
        paraF << th
        paraF << be
        paraF << su
        paraF << sm
        paraF << ve
        print 'Done with mesh'