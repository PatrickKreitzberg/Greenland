from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import pyqtgraph as pg
from pyqtgraph import LegendItem
import PyQt4
from pylab import plot,show,ion,subplots, sqrt
from dolfin import project
from MyLegend import *
from scipy.optimize import curve_fit
import numpy as np
from ..dataset_objects import *
from ..peakdetect import *

from ..gui import *
import scipy.signal as signal

def gauss(x, *p):
    A, mu, sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def gaussian(x, A, x0, sig):
    return A*math.exp(-(x-x0)**2/(2.0*sig**2))

def fit(p,x):
    return np.sum([gaussian(x, p[i*3],p[i*3+1],p[i*3+2])
                   for i in xrange(len(p)/3)],axis=0)

def runStaticPlot():
    staticPlotWindow = QtGui.QMainWindow(mw)
    staticPlotWindow.setWindowTitle('Static Plotter')
    dummyWidget = QtGui.QWidget()
    staticPlotWindow.setCentralWidget(dummyWidget)
    layout = QtGui.QVBoxLayout()
    dummyWidget.setLayout(layout)
    plt1 = pg.PlotWidget()
    plt2 = pg.PlotWidget()
    plt3 = pg.PlotWidget()
    layout.addWidget(plt1)
    layout.addWidget(plt2)
    layout.addWidget(plt3)
    legend1 = pg.LegendItem(offset=(-1, 1))
    legend2 = pg.LegendItem(offset=(-1, 1))
    legend3 = pg.LegendItem(offset=(-1, 1))

    if surfaceCheck.checkState() == 2:
        surfPlt  = plt1.getPlotItem().plot(surface.distanceData, surface.pathData, pen=surface.pen)
        legend1.addItem(surfPlt, 'Surface(m)')

    if vWidthCheck.checkState() == 2:
        vWidthPlt = plt1.getPlotItem().plot(velocityWidth.distanceData, velocityWidth.pathData, pen=velocityWidth.pen)
        legend1.addItem(vWidthPlt, 'Vel. Width(m)')

    if bedCheck.checkState() == 2:
        bedPlt   = plt1.getPlotItem().plot(bed.distanceData, bed.pathData, pen=bed.pen)
        legend1.addItem(bedPlt, 'Bed(m)')

    if thicknessCheck.checkState() == 2:
        thickPlt = plt1.getPlotItem().plot(thickness.distanceData, thickness.pathData, pen=thickness.pen)
        legend1.addItem(thickPlt, 'Thickness(m)')

    if velocityCheck.checkState() == 2:
        velocityPlt = plt2.getPlotItem().plot(velocity.distanceData, velocity.pathData, pen=velocity.pen)
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





        legend2.addItem(velocityPlt, 'Velocity(m/yr)')


    if smbCheck.checkState() == 2:
        smbPlt = plt3.getPlotItem().plot(smb.distanceData, smb.pathData, pen=smb.pen)
        legend3.addItem(smbPlt, 'SMB(m)')

    legend1.setParentItem(plt1.getPlotItem())

    legend2.setParentItem(plt2.getPlotItem())

    legend3.setParentItem(plt3.getPlotItem())
    staticPlotWindow.show()

