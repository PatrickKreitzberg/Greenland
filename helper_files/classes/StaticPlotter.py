from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import pyqtgraph as pg
from pyqtgraph import LegendItem
import PyQt4
from pylab import plot,show,ion,subplots
from scipy import sqrt
from dolfin import project
from scipy.optimize import curve_fit
import numpy as np
import scipy.signal as signal

# LOCAL IMPORTS
from MyLegend import *
from ..dataset_objects import *
from ..peakdetect import *
from ..data_functions import interpolateData
from ..gui import *


def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gaussian(x, A, x0, sig):
    return A*math.exp(-(x-x0)**2/(2.0*sig**2))

def fit(p,x):
    return np.sum([gaussian(x, p[i*3],p[i*3+1],p[i*3+2])
                   for i in xrange(len(p)/3)],axis=0)

class StaticPlot(QtGui.QMainWindow):
    def __init__(self, parent):
        self.parent = parent
        QtGui.QMainWindow.__init__(self, self.parent)
        self.setWindowTitle('Static Plotter')
        self.mainWidget = QtGui.QWidget()
        self.setCentralWidget(self.mainWidget)
        self.mainLayout = QtGui.QHBoxLayout()
        self.mainWidget.setLayout(self.mainLayout)

        # LEFT PANEL
        self.leftPanelWidget = QtGui.QWidget()

        self.leftPanelLayout = QtGui.QVBoxLayout()
        self.leftPanelWidget.setLayout(self.leftPanelLayout)
        self.plt1 = pg.PlotWidget()
        self.plt2 = pg.PlotWidget()
        self.plt3 = pg.PlotWidget()
        self.leftPanelLayout.addWidget(self.plt1)
        self.leftPanelLayout.addWidget(self.plt2)
        self.leftPanelLayout.addWidget(self.plt3)
        self.legend1 = pg.LegendItem(offset=(-1, 1))
        self.legend2 = pg.LegendItem(offset=(-1, 1))
        self.legend3 = pg.LegendItem(offset=(-1, 1))

        # RIGHT PANEL
        self.rightPanelWidget = QtGui.QWidget()
        self.rightPanelLayout = QtGui.QGridLayout()
        self.rightPanelWidget.setLayout(self.rightPanelLayout)

        self.resLabel = QtGui.QLabel('Spatial Resolution(m)')
        self.warningLabel = QtGui.QLabel('Higher resolution takes time')
        self.resLineEdit = QtGui.QLineEdit()
        self.plotButt = QtGui.QPushButton('Plot')
        self.errorLabel = QtGui.QLabel('')
        self.plotButt.clicked.connect(self.run)

        self.allCheck        = QtGui.QCheckBox('Plot All')
        self.velocityCheck   = QtGui.QCheckBox('Velocity')
        self.vWidthCheck     = QtGui.QCheckBox('Velocity Width')
        self.smbCheck        = QtGui.QCheckBox('SMB')
        self.surfaceCheck    = QtGui.QCheckBox('Surface')
        self.bedCheck        = QtGui.QCheckBox('Bed')
        self.thicknessCheck  = QtGui.QCheckBox('Thickness')

        self.checkBoxW = QtGui.QWidget()
        self.checkBLayout = QtGui.QVBoxLayout()
        self.checkBoxW.setLayout(self.checkBLayout)
        self.allCheck.setTristate(False)
        self.velocityCheck.setTristate(False)
        self.vWidthCheck.setTristate(False)
        self.smbCheck.setTristate(False)
        self.surfaceCheck.setTristate(False)
        self.bedCheck.setTristate(False)
        self.thicknessCheck.setTristate(False)
        self.allCheck.setCheckState(2)
        self.velocityCheck.setCheckState(2)
        self.vWidthCheck.setCheckState(2)
        self.smbCheck.setCheckState(2)
        self.surfaceCheck.setCheckState(2)
        self.bedCheck.setCheckState(2)
        self.thicknessCheck.setCheckState(2)
        self.checkBLayout.addWidget(QtGui.QLabel('Plot Checked Data:'))
        self.checkBLayout.addWidget(self.allCheck)
        self.checkBLayout.addWidget(self.velocityCheck)
        self.checkBLayout.addWidget(self.vWidthCheck)
        self.checkBLayout.addWidget(self.smbCheck)
        self.checkBLayout.addWidget(self.surfaceCheck)
        self.checkBLayout.addWidget(self.bedCheck)
        self.checkBLayout.addWidget(self.thicknessCheck)
        self.checkBLayout.setSpacing(0)
        self.allCheck.stateChanged.connect(self.allCheckChange)

        self.rightPanelLayout.addWidget(self.resLabel,     0, 0)
        self.rightPanelLayout.addWidget(self.resLineEdit,  0, 1)
        self.rightPanelLayout.addWidget(self.warningLabel, 1, 0, 1, 2)
        self.rightPanelLayout.addWidget(self.plotButt,     2, 0, 1, 2)
        self.rightPanelLayout.addWidget(self.checkBoxW,    3, 0, 1, 2)
        self.rightPanelLayout.addWidget(self.errorLabel,   3, 0, 1, 2)
        self.rightPanelLayout.setAlignment(QtCore.Qt.AlignTop)

        self.mainLayout.addWidget(self.leftPanelWidget)
        self.mainLayout.addWidget(self.rightPanelWidget)


        # if surfaceCheck.checkState() == 2:
        #     surfPlt  = self.plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=surface.pen)
        #     self.legend1.addItem(surfPlt, 'Surface(m)')
        #
        # if vWidthCheck.checkState() == 2:
        #     vWidthPlt = self.plt1.getPlotItem().plot(velocityWidth.distanceData, velocityWidth.pathData, pen=velocityWidth.pen)
        #     self.legend1.addItem(vWidthPlt, 'Vel. Width(m)')
        #
        # if bedCheck.checkState() == 2:
        #     bedPlt   = self.plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bed.pen)
        #     self.legend1.addItem(bedPlt, 'Bed(m)')
        #
        # if thicknessCheck.checkState() == 2:
        #     thickPlt = self.plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=thickness.pen)
        #     self.legend1.addItem(thickPlt, 'Thickness(m)')
        #
        # if velocityCheck.checkState() == 2:
        #     velocityPlt = self.plt2.getPlotItem().plot(velocity.distanceData, velocity.pathData, pen=velocity.pen)
        #     self.legend2.addItem(velocityPlt, 'Velocity(m/yr)')
        #
        # if smbCheck.checkState() == 2:
        #     smbPlt = self.plt3.getPlotItem().plot(smb.distanceData, smb.pathData, pen=smb.pen)
        #     self.legend3.addItem(smbPlt, 'SMB(m)')

        self.legend1.setParentItem(self.plt1.getPlotItem())
        self.legend2.setParentItem(self.plt2.getPlotItem())
        self.legend3.setParentItem(self.plt3.getPlotItem())
        self.show()

    def allCheckChange(self, e):
        if self.allCheck.checkState() == 2:
            self.velocityCheck.setCheckState(2)
            self.vWidthCheck.setCheckState(2)
            self.smbCheck.setCheckState(2)
            self.surfaceCheck.setCheckState(2)
            self.bedCheck.setCheckState(2)
            self.thicknessCheck.setCheckState(2)
        else:
            self.velocityCheck.setCheckState(0)
            self.vWidthCheck.setCheckState(0)
            self.smbCheck.setCheckState(0)
            self.surfaceCheck.setCheckState(0)
            self.bedCheck.setCheckState(0)
            self.thicknessCheck.setCheckState(0)

    def run(self):
        try:
            dr = float(self.resLineEdit.text())
            self.errorLabel.setText('')
        except ValueError:
            self.errorLabel.setText('ERROR: must enter resolution')
            return -1
        velocity.pathPlotItem.clear()
        surface.pathPlotItem.clear()
        smb.pathPlotItem.clear()
        bed.pathPlotItem.clear()
        velocityWidth.pathPlotItem.clear()
        thickness.pathPlotItem.clear()
        if len(markers) > 0:
            print 'Plotting...'
            # nbed, nsurf, nv, nsmb, nvelWidth, linePoints, graphX = interpolateData(True)
            dataSetsToPopulate = []

            if self.velocityCheck.checkState() == 2:
                dataSetsToPopulate.append('vel')
            if self.smbCheck.checkState() == 2:
                dataSetsToPopulate.append('smb')
            if self.surfaceCheck.checkState() == 2:
                dataSetsToPopulate.append('sur')
            if self.bedCheck.checkState() == 2:
                dataSetsToPopulate.append('bed')
            if self.thicknessCheck.checkState() == 2:
                dataSetsToPopulate.append('thk')
            if self.vWidthCheck.checkState() == 2:
                dataSetsToPopulate.append('wth')

            interpolateData(False, dr, dataSetsToPopulate)
            if self.velocityCheck.checkState() == 2:
                print 'plot velocity!'
                velocityPlt = self.plt2.getPlotItem().plot(velocity.distanceData, velocity.pathData, pen=velocity.pen)
                self.legend2.addItem(velocityPlt, 'Velocity(m/yr)')
                velocity.pathPlotItem.setData(velocity.distanceData          , velocity.pathData)

            if self.smbCheck.checkState() == 2:
                self.smbPlt = self.plt3.getPlotItem().plot(smb.distanceData, smb.pathData, pen=smb.pen)
                self.legend3.addItem(self.smbPlt, 'SMB(m)')
                smb.pathPlotItem.setData(smb.distanceData                    , smb.pathData)

            if self.surfaceCheck.checkState() == 2:
                self.surfPlt  = self.plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=surface.pen)
                self.legend1.addItem(self.surfPlt, 'Surface(m)')
                surface.pathPlotItem.setData(surface.distanceData            , surface.pathData)

            if self.bedCheck.checkState() == 2:
                self.bedPlt = self.plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bed.pen)
                self.legend1.addItem(self.bedPlt, 'Bed(m)')
                bed.pathPlotItem.setData(bed.distanceData                    , bed.pathData)

            if self.thicknessCheck.checkState() == 2:
                self.thickPlt = self.plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=thickness.pen)
                self.legend1.addItem(self.thickPlt, 'Thickness(m)')
                thickness.pathPlotItem.setData(thickness.distanceData, thickness.pathData)

            if self.vWidthCheck.checkState() == 2:
                self.vWidthPlt = self.plt1.getPlotItem().plot(velocityWidth.distanceData, velocityWidth.pathData, pen=velocityWidth.pen)
                self.legend1.addItem(self.vWidthPlt, 'Vel. Width(m)')
                velocityWidth.pathPlotItem.setData(velocityWidth.distanceData, velocityWidth.pathData)

            pg.QtGui.QApplication.processEvents()
            print 'Plotting done.'

# def runStaticPlot():
#     staticPlotWindow = QtGui.QMainWindow(mw)
#     staticPlotWindow.setWindowTitle('Static Plotter')
#     dummyWidget = QtGui.QWidget()
#     staticPlotWindow.setCentralWidget(dummyWidget)
#     layout = QtGui.QVBoxLayout()
#     dummyWidget.setLayout(layout)
#     plt1 = pg.PlotWidget()
#     plt2 = pg.PlotWidget()
#     plt3 = pg.PlotWidget()
#     layout.addWidget(plt1)
#     layout.addWidget(plt2)
#     layout.addWidget(plt3)
#     legend1 = pg.LegendItem(offset=(-1, 1))
#     legend2 = pg.LegendItem(offset=(-1, 1))
#     legend3 = pg.LegendItem(offset=(-1, 1))
#
#     if surfaceCheck.checkState() == 2:
#         surfPlt  = plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=surface.pen)
#         legend1.addItem(surfPlt, 'Surface(m)')
#
#     if vWidthCheck.checkState() == 2:
#         vWidthPlt = plt1.getPlotItem().plot(velocityWidth.distanceData, velocityWidth.pathData, pen=velocityWidth.pen)
#         legend1.addItem(vWidthPlt, 'Vel. Width(m)')
#
#     if bedCheck.checkState() == 2:
#         bedPlt   = plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bed.pen)
#         legend1.addItem(bedPlt, 'Bed(m)')
#
#     if thicknessCheck.checkState() == 2:
#         thickPlt = plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=thickness.pen)
#         legend1.addItem(thickPlt, 'Thickness(m)')
#
#     if velocityCheck.checkState() == 2:
#         velocityPlt = plt2.getPlotItem().plot(velocity.distanceData, velocity.pathData, pen=velocity.pen)
        # dv = velocity.pathData[1:] - velocity.pathData[:-1]
        # dv2 = dv[1:] - dv[:-1]
        #
        # vx = []
        # for i in range(len(dv)):
        #     if math.fabs(dv[i-1]) < 0.01:
        #         vx.append(velocity.distanceData[i])
        # for px in vx:
        #     # vLine = pg.InfiniteLine(angle=90, movable=False, pos=px)
        #     plt1.getPlotItem().addItem(pg.InfiniteLine(angle=90, movable=False, pos=px))

            #
            #   Adds line where slop ~0
            #

        # if len(velocity.pathData)%2 ==0:
        #     wl = min(11, len(velocity.pathData)-1)
        # else:
        #     wl = min(11, len(velocity.pathData))
        # bleh = signal.savgol_filter(velocity.pathData,wl,1)
        # plt2.getPlotItem().plot(velocity.distanceData, bleh, pen=whitePlotPen)
        # dv = bleh[1:] - bleh[:-1]
        # plt3.getPlotItem().plot(velocity.distanceData[:-1], dv, pen=redPlotPen)
        # plt3.getPlotItem().plot(velocity.distanceData, bleh, pen=whitePlotPen)
        # mx = 0
        # mxdv = 0
        # mndv = 0
        # for i in range(len(bleh)-1):
        #     if bleh[i] > bleh[mx]:
        #         mx = i
        #     if dv[i] > dv[mxdv]:
        #         mxdv = i
        #     if dv[i] < dv[mndv]:
        #         mndv = i
        #     if math.fabs(dv[i]) < 0.1:
        #         plt2.addItem(pg.InfiniteLine(angle=90, movable=False, pos=velocity.distanceData[i]))
        # plt2.addItem(pg.InfiniteLine(angle=90, movable=False, pos=velocity.distanceData[mx]))
        # plt2.addItem(pg.InfiniteLine(angle=90, movable=False, pos=velocity.distanceData[mxdv]))
        # plt2.addItem(pg.InfiniteLine(angle=90, movable=False, pos=velocity.distanceData[mndv]))
        #
        #   Line at maximum/minimum dv and maximum v
        #

        # LETS FIND THE PEAKS


        # x = linspace(0,100,500)
        # y = np.sin(x)
        # lah = max((len(velocity.pathData)/velocity.distanceData[-1])*200, int(np.floor(float(model_res_lineEdit.text()))))  # points/meter * 100 meters
        # peaks = peakdetect(velocity.pathData, velocity.distanceData, lookahead=lah)
        # print 'peaks', peaks
        # print 'max peaks', peaks[0]
        # print 'min peaks', peaks[1]
        # for p in peaks[0]:
        #     print 'maxpeak: ', p[1]
        #     plt2.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0]))
        #
        # for p in peaks[1]:
        #     if len(p[0]) > 0:
        #         print 'minpeak: ', p[1]
        #         plt2.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0], pen=redPlotPen))
        #
        #
        #
        # peaksdv = peakdetect(dv, velocity.distanceData[:-1], lookahead=lah)
        #
        # for p in peaksdv:
        #     if len(p[0]) > 0:
        #         print 'peak: ', p
        #         plt3.addItem(pg.InfiniteLine(angle=90, movable=False, pos=p[0]))




    #
    #     legend2.addItem(velocityPlt, 'Velocity(m/yr)')
    #
    #
    # if smbCheck.checkState() == 2:
    #     smbPlt = plt3.getPlotItem().plot(smb.distanceData, smb.pathData, pen=smb.pen)
    #     legend3.addItem(smbPlt, 'SMB(m)')
    #
    # legend1.setParentItem(plt1.getPlotItem())
    #
    # legend2.setParentItem(plt2.getPlotItem())
    #
    # legend3.setParentItem(plt3.getPlotItem())
    # staticPlotWindow.show()

