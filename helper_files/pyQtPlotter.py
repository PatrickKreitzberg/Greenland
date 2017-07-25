from pyqtgraph.Qt import QtCore, QtGui #, QtWidgets
import pyqtgraph as pg
from pyqtgraph import LegendItem
import PyQt4
from pylab import plot,show,ion,subplots
from dolfin import project
from myLegend import myLegend

class pyqtplotter(object):
    def __init__(self, strs, mesh, plt1, plt2, plt3):
        #####################################################
        ####         CREATE PENS                         ####
        #####################################################
        whitePen = pg.mkPen(color=(255, 255, 255), width=2)
        redPen   = pg.mkPen(color=(100,   0,   0), width=2)
        bluePen  = pg.mkPen(color=(  0,   0, 255), width=2)
        greenPen = pg.mkPen(color=( 76, 153,   0), width=2)
        tealPen  = pg.mkPen(color=(  0, 204, 204), width=2)
        pinkPen  = pg.mkPen(color=(153, 0,  153), width=2)


        self.strs = strs # stresses
        self.mesh = mesh # mesh (2-d)
        self.x = mesh.coordinates().flatten()
        mesh_min = mesh.coordinates().min()
        mesh_max = mesh.coordinates().max()
        S = project(strs.H0 + strs.B)
        Bmin = self.strs.B.vector().min()
        Smax = S.vector().max()
        TD = project(self.strs.tau_d_plot)
        TB = project(self.strs.tau_b_plot)
        TX = project(self.strs.tau_xx_plot)
        TY = project(self.strs.tau_xy_plot)
        TZ = project(self.strs.tau_xz_plot)

        self.plt1 = plt1
        self.plt1.showGrid(x=True, y=True)

        self.plt1.getPlotItem().getViewBox().setRange(xRange=(mesh_min, mesh_max), yRange=(Bmin * 1.05, Smax * 1.2))
        self.ph0   =  self.plt1.plot(self.x,self.strs.B.compute_vertex_values(), pen=bluePen)
        self.ph100 =  self.plt1.plot(self.x,S.compute_vertex_values(), pen=redPen)
        self.ph1   =  self.plt1.plot(self.x,S.compute_vertex_values(), pen=whitePen)
        # self.ax[0].legend([r"$B$",r"$S_o$",r"$S$"])

        # self.legend1 = self.plt1.getPlotItem().addLegend(offset=(-100,50))
        #
        self.legend1 = myLegend(offset=(-50,50))
        self.legend1.setParentItem(self.plt1.graphicsItem())
        self.legend1.addItem(self.ph0,   '<i>B</i>')
        self.legend1.addItem(self.ph100, '<i>S</i><sub>o</sub>')
        self.legend1.addItem(self.ph1,   '<i>S</i>')
        # self.legend1.getViewBox().setBackgroundColor((255,255,255))

        self.plt2 = plt2
        self.plt2.showGrid(x=True, y=True)
        self.plt2.getPlotItem().getViewBox().setRange(xRange=(mesh_min, mesh_max), yRange=(TD.vector().array().min(), TD.vector().array().max()))
        # self.ax[1].legend([r"$\tau_d$",r"$\tau_b$",r"$\tau_{xx}$",r"$\tau_{xy}$",r"$\tau_{xz}$"])
        self.ph2 = self.plt2.plot(self.x,TD.compute_vertex_values(), pen=bluePen)
        self.ph3 = self.plt2.plot(self.x,TB.compute_vertex_values(), pen=greenPen)
        self.ph4 = self.plt2.plot(self.x,TX.compute_vertex_values(), pen=redPen)
        self.ph5 = self.plt2.plot(self.x,TY.compute_vertex_values(), pen=tealPen)
        self.ph6 = self.plt2.plot(self.x,TZ.compute_vertex_values(), pen=pinkPen)

        self.legend2 = myLegend(offset=(-50,50))
        self.legend2.setParentItem(self.plt2.graphicsItem())
        self.legend2.addItem(self.ph2, '&tau;<sub>d</sub>')
        self.legend2.addItem(self.ph3, '&tau;<sub>b</sub>')
        self.legend2.addItem(self.ph4, '&tau;<sub>xx</sub>')
        self.legend2.addItem(self.ph5, '&tau;<sub>xy</sub>')
        self.legend2.addItem(self.ph6, '&tau;<sub>xz</sub>')


        self.plt3 = plt3
        self.plt3.showGrid(x=True, y=True)
        self.plt3.getPlotItem().getViewBox().setRange(xRange=(mesh_min, mesh_max), yRange=(0, 500))
        # self.ax[2].legend([r"$u_s$",r"$u_b$"])
        self.legend3 = myLegend(offset=(-50,50))
        self.legend3.setParentItem(self.plt3.graphicsItem())

        us = project(self.strs.u(0))
        ub = project(self.strs.u(1))

        self.ph7 = self.plt3.plot(self.x, us.compute_vertex_values(), pen=bluePen)
        self.ph8 = self.plt3.plot(self.x, ub.compute_vertex_values(), pen=greenPen)
        self.legend3.addItem(self.ph7, '&mu;<sub>s</sub>')
        self.legend3.addItem(self.ph8, '&mu;<sub>b</sub>')


    def refresh_plot(self):
        BB = self.strs.B.compute_vertex_values()
        HH = self.strs.H0.compute_vertex_values()
        TD = project(self.strs.tau_d_plot)
        TB = project(self.strs.tau_b_plot)
        TX = project(self.strs.tau_xx_plot)
        TY = project(self.strs.tau_xy_plot)
        TZ = project(self.strs.tau_xz_plot)
        us = project(self.strs.u(0))
        ub = project(self.strs.u(1))

        self.ph0.setData(self.x, BB)
        self.ph1.setData(self.x, (BB + HH))

        self.ph2.setData(self.x, TD.compute_vertex_values())
        self.ph3.setData(self.x, TB.compute_vertex_values())
        self.ph4.setData(self.x, TX.compute_vertex_values())
        self.ph5.setData(self.x, TY.compute_vertex_values())
        self.ph6.setData(self.x, TZ.compute_vertex_values())

        self.ph7.setData(self.x, us.compute_vertex_values())
        self.ph8.setData(self.x, ub.compute_vertex_values())



