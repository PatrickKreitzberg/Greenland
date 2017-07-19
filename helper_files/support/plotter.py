from pylab import plot,show,ion,subplots
from dolfin import project

class plotter(object):
    def __init__(self,strs,mesh):
        self.strs = strs
        self.mesh = mesh

        # Initialize plot
        ion()
        self.fig, self.ax = subplots(nrows=3,sharex=True,figsize=(10,12))
        self.x = mesh.coordinates()
        mesh_min = mesh.coordinates().min()
        mesh_max = mesh.coordinates().max()

        S    = project(strs.H0+strs.B)
        Bmin = self.strs.B.vector().min()
        Smax = S.vector().max()

        # ax[0] is first plot
        # ax[1] is second plot etc
        self.ph0,   =  self.ax[0].plot(self.x,self.strs.B.compute_vertex_values(),'b-',linewidth=2)
        self.ph100, =  self.ax[0].plot(self.x,S.compute_vertex_values(),'r-',linewidth=2)
        self.ph1,   =  self.ax[0].plot(self.x,S.compute_vertex_values(),'k-',linewidth=2)
        self.ax[0].legend([r"$B$",r"$S_o$",r"$S$"])

        self.ax[0].set_xlim(mesh_min,mesh_max)

        self.ax[0].set_ylim(Bmin*1.05,Smax*1.2)

        self.ax[0].grid()

        TD = project(self.strs.tau_d_plot)
        TB = project(self.strs.tau_b_plot)
        TX = project(self.strs.tau_xx_plot)
        TY = project(self.strs.tau_xy_plot)
        TZ = project(self.strs.tau_xz_plot)

        self.ph2, = self.ax[1].plot(self.x,TD.compute_vertex_values(),linewidth=2)
        self.ph3, = self.ax[1].plot(self.x,TB.compute_vertex_values(),linewidth=2)
        self.ph4, = self.ax[1].plot(self.x,TX.compute_vertex_values(),linewidth=2)
        self.ph5, = self.ax[1].plot(self.x,TY.compute_vertex_values(),linewidth=2)
        self.ph6, = self.ax[1].plot(self.x,TZ.compute_vertex_values(),linewidth=2)

        self.ax[1].set_xlim(mesh_min,mesh_max)
        self.ax[1].set_ylim(TD.vector().array().min(),TD.vector().array().max())
        self.ax[1].legend([r"$\tau_d$",r"$\tau_b$",r"$\tau_{xx}$",r"$\tau_{xy}$",r"$\tau_{xz}$"])
        self.ax[1].grid()

        us = project(self.strs.u(0))
        ub = project(self.strs.u(1))

        self.ph7, = self.ax[2].plot(self.x,us.compute_vertex_values(),linewidth=2)
        self.ph8, = self.ax[2].plot(self.x,ub.compute_vertex_values(),linewidth=2)

        self.ax[2].set_xlim(mesh_min,mesh_max)
        self.ax[2].set_ylim(0,500)
        self.ax[2].legend([r"$u_s$",r"$u_b$"])
        self.ax[2].grid()

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

        self.ph0.set_ydata(BB)
        self.ph1.set_ydata((BB + HH))

        self.ph2.set_ydata(TD.compute_vertex_values())
        self.ph3.set_ydata(TB.compute_vertex_values())
        self.ph4.set_ydata(TX.compute_vertex_values())
        self.ph5.set_ydata(TY.compute_vertex_values())
        self.ph6.set_ydata(TZ.compute_vertex_values())

        self.ph7.set_ydata(us.compute_vertex_values())
        self.ph8.set_ydata(ub.compute_vertex_values())

        self.fig.canvas.draw()
