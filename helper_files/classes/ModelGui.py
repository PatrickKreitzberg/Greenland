import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import fenics as fc
import h5py
from PyQt4 import QtCore
from time import *
from dolfin import *
from ..gui import *
from ..dataset_objects import *
from ..data_functions import *
from ModelPlotter import *
from helper_files.support.expressions import *
from helper_files.support.fenics_optimizations import *
from helper_files.support.momentum import *
from scipy.interpolate import interp1d

class ModelGUI(QtGui.QMainWindow):
    def __init__(self, parent):
        self.hdf_name = None
        self.pPlt = None
        self.parent = parent
        QtGui.QMainWindow.__init__(self, self.parent)

        self.cw = QtGui.QWidget()
        self.setCentralWidget(self.cw)
        self.horLayout = QtGui.QHBoxLayout()
        self.cw.setLayout(self.horLayout)


        # RIGHT PANEL

        self.rightPanelW   = QtGui.QWidget()
        self.rightPanelLay = QtGui.QGridLayout()
        self.rightPanelW.setLayout(self.rightPanelLay)

            # GUI TIME COMPONENTS
        # self.inputWidget = QtGui.QWidget()
        # self.inputContainer  = QtGui.QGridLayout()
        # self.inputWidget.setLayout(self.inputContainer)
        self.tEndLabel       = QtGui.QLabel('t_end(yr):')
        self.tEndLineEdit    = QtGui.QLineEdit('20000')
        self.tStepLabel      = QtGui.QLabel('t_step(yr):')
        self.tStepLineEdit   = QtGui.QLineEdit('10')
        self.tCurrent        = QtGui.QLabel('Current year: ')
        self.rightPanelLay.setSpacing(4)
        self.sptlResLabel    = QtGui.QLabel('Spatial Res(m)')
        self.sptlResLineEdit = QtGui.QLineEdit('2000')

        self.rightPanelLay.addWidget(self.tEndLabel,       0, 0)
        self.rightPanelLay.addWidget(self.tEndLineEdit,    0, 1)
        self.rightPanelLay.addWidget(self.tStepLabel,      1, 0)
        self.rightPanelLay.addWidget(self.tStepLineEdit,   1, 1)
        self.rightPanelLay.addWidget(self.sptlResLabel,    2, 0)
        self.rightPanelLay.addWidget(self.sptlResLineEdit, 2, 1)
        self.rightPanelLay.addWidget(self.tCurrent,        3, 0, 1, 2)


            # SAVE FILE
        # self.saveFileW = QtGui.QWidget()
        # self.saveFileLay = QtGui.QHBoxLayout()
        # self.saveFileW.setLayout(self.saveFileLay)
        self.saveFileLabel = QtGui.QLabel('Out File')
        self.saveFileLineEdit = QtGui.QLineEdit('./data/outModel.h5')
        self.rightPanelLay.addWidget(self.saveFileLabel,    4, 0)
        self.rightPanelLay.addWidget(self.saveFileLineEdit, 4, 1)
        self.rightPanelLay.setAlignment(QtCore.Qt.AlignTop)
        # self.saveFileLay.addWidget(self.saveFileLabel)
        # self.saveFileLay.addWidget(self.saveFileLineEdit)

            # BUTTONS
        self.runButt = QtGui.QPushButton('Run Model')
        self.pauseButt = QtGui.QPushButton('Pause Model')


            # ADD TO LAYOUT

        # self.rightPanelLay.addWidget(self.inputWidget)
        # self.rightPanelLay.addWidget(self.saveFileW, 4, 0, 1, 2)
        self.rightPanelLay.addWidget(self.runButt,   5, 0, 1, 2)
        self.rightPanelLay.addWidget(self.pauseButt, 6, 0, 1, 2)
        self.rightPanelW.setMaximumWidth(300)
        self.pauseButt.setEnabled(False)

            # SAVE 1D MESH

        self.save1DMeshLabel = QtGui.QLabel('Mesh File:')
        self.save1DMeshLineEdit = QtGui.QLineEdit('./data/out1DMesh.h5')
        self.save1DMeshButton = QtGui.QPushButton('Save 1D Mesh')
        self.save1DMeshButton.clicked.connect(self.write1DMesh)
        self.rightPanelLay.addWidget(self.save1DMeshLabel, 7, 0)
        self.rightPanelLay.addWidget(self.save1DMeshLineEdit, 7, 1)
        self.rightPanelLay.addWidget(self.save1DMeshButton, 8, 0, 1, 2)



        # LEFT SIDE
        self.leftPanelW = QtGui.QWidget()
        self.leftPanelLay = QtGui.QVBoxLayout()  # CENTRAL WIDGET LAYOUT (layout of the central widget)
        self.leftPanelW.setLayout(self.leftPanelLay)
        self.plt1 = pg.PlotWidget()
        self.plt2 = pg.PlotWidget()
        self.plt3 = pg.PlotWidget()
        self.slider = QtGui.QSlider(QtCore.Qt.Horizontal)
        self.slider.setEnabled(False)
        self.slider.sliderChange = self.sliderChange
        self.sliderLabel = QtGui.QLabel('Year:')
        self.leftPanelLay.addWidget(self.plt1)
        self.leftPanelLay.addWidget(self.plt2)
        self.leftPanelLay.addWidget(self.plt3)
        self.leftPanelLay.addWidget(self.slider)
        self.leftPanelLay.addWidget(self.sliderLabel)


        # ADD SIDES TO MAIN WINDOW

        self.horLayout.addWidget(self.leftPanelW)
        self.horLayout.addWidget(self.rightPanelW)

        self.runButt.clicked.connect(self.runModelButt)
        self.pauseButt.clicked.connect(self.pause)
        self.showMaximized()
        self.show()
        self.closeEvent = self.windowClosed

    def write1DMesh(self):
        print "Writing 1D mesh to file ", self.save1DMeshLineEdit.text()
        self.initModel(False)
        print "Done writing 1D mesh to file ", self.save1DMeshLineEdit.text()

    def runModelButt(self):
        self.initModel(True)

    def windowClosed(self, e):
        if self.pPlt:
            self.pPlt.run = False
            self.pPlt.closePlots()

    def sliderChange(self, e):
        key = self.pPlt.outFTimeList[self.slider.value()]
        self.pPlt.ph1.setData(self.pPlt.x, self.pPlt.outF['surface'][key])
        self.sliderLabel.setText('Year: ' + key)

    def pause(self):
        if self.pPlt.run:
            self.pPlt.run = False
            self.pauseButt.setText('Resume')
            self.slider.setEnabled(True)
            self.slider.setTickInterval(10)#int(self.tStepLineEdit.text()))
            self.slider.setRange(0,len(self.pPlt.outFTimeList)-1)

        else:
            self.pPlt.run = True
            self.slider.setEnabled(False)
            self.pauseButt.setText('Pause')
            self.runLoop()

    def initModel(self, run):
        if len(vpts) > 0:
            try:
                self.runButt.setEnabled(False)
                self.dr = float(self.sptlResLineEdit.text())  # dr = 150
                self.pauseButt.setEnabled(True)
                interpolateData(True, self.dr)

                ###########################################
                ###  INTERPOLATE DATA SO EVEN INTERVAL  ###
                ###########################################
                '''
                Interpolating data at an even interval so the data points align with the 
                invterval mesh.
                '''
                self.thickness1dInterp = interp1d(thickness.distanceData, thickness.pathData)
                self.bed1dInterp       = interp1d(bed.distanceData, bed.pathData)
                self.surface1dInterp   = interp1d(surface.distanceData, surface.pathData)
                self.smb1dInterp       = interp1d(smb.distanceData, smb.pathData)
                self.velocity1dInterp  = interp1d(velocity.distanceData, velocity.pathData)

                # N is the number of total data points including the last
                # Data points on interval [0, N*dr] inclusive on both ends

                self.N = int(np.floor(bed.distanceData[-1] / float(self.dr)))  # length of path / resolution


                self.x = np.arange(0, (self.N + 1) * self.dr, self.dr)  # start point, end point, number of segments. END POINT NOT INCLUDED!
                self.mesh = fc.IntervalMesh(self.N, 0, self.dr * self.N)  # number of cells, start point, end point

                self.thicknessModelData = self.thickness1dInterp(self.x)
                self.bedModelData       = self.bed1dInterp(self.x)
                self.surfaceModelData   = self.surface1dInterp(self.x)
                self.smbModelData       = self.smb1dInterp(self.x)
                self.velocityModelData  = self.velocity1dInterp(self.x)

                self.THICKLIMIT = 10.  # Ice is never less than this thick
                self.H = self.surfaceModelData - self.bedModelData
                self.surfaceModelData[self.H <= self.THICKLIMIT] = self.bedModelData[self.H <= self.THICKLIMIT]

                # FIXME the intervalMesh is consistantly 150 between each datapoint this not true for the data being sent
                # FIXME this happens at line self.x = ...


                self.hdf_name = str(self.save1DMeshLineEdit.text())
                self.hfile = fc.HDF5File(self.mesh.mpi_comm(), self.hdf_name, "w")
                self.V = fc.FunctionSpace(self.mesh, "CG", 1)
                self.functThickness = fc.Function(self.V, name="Thickness")
                self.functBed       = fc.Function(self.V, name="Bed")
                self.functSurface   = fc.Function(self.V, name="Surface")
                self.functSMB       = fc.Function(self.V, name='SMB')
                self.functVelocity  = fc.Function(self.V, name='Velocity')

                surface.pathPlotItem.setData(self.x, self.surfaceModelData)
                pg.QtGui.QApplication.processEvents()

                self.functThickness.vector()[:] = self.thicknessModelData
                self.functBed.vector()[:]       = self.bedModelData
                self.functSurface.vector()[:]   = self.surfaceModelData
                self.functSMB.vector()[:]       = self.smbModelData
                self.functVelocity.vector()[:]  = self.velocityModelData

                self.hfile.write(self.functThickness.vector(), "/thickness")
                self.hfile.write(self.functBed.vector(), "/bed")
                self.hfile.write(self.functSurface.vector(), "/surface")
                self.hfile.write(self.functSMB.vector(), "/smb")
                self.hfile.write(self.functVelocity.vector(), "/velocity")
                self.hfile.write(self.mesh, "/mesh")
                self.hfile.close()
                if run:
                    self.runModel()

            except ValueError:
                print 'ERROR: Must have valid spatial resolution.'

    def runModel(self):

        ##########################################################
        ################           FILES         #################
        ##########################################################

        #FIXME Probably dont have to save then open mesh
        self.mesh = Mesh()
        self.in_file  = HDF5File(self.mesh.mpi_comm(), self.hdf_name, "r")  #mesh.mpi_comm() is ussed to read in parallel?


        # out_file = HDF5File(mesh.mpi_comm(),"./output_data/peterman.h5","w")
        # cell_indices: self-explanatory
        # coordinates:  self-explanatory
        # topology:     Shows which nodes are linked together

        ##########################################################
        ################           GUI           #################
        ##########################################################



        ##########################################################
        ################           MESH          #################
        ##########################################################
        self.in_file.read(self.mesh,"/mesh", False)


        # H5FILE Data:
        # bed and surface are datasets shape (378,)
        # mesh is a group with datasets
        # cell_indices  Shape (377,)    Type i8     Array [0 1 2 ... n-1 n]              incremented by 1
        # coordinates   Shape (378,1)   Type f8     Array [[0.] [1000.] ... [ 377000.]]  incremented by 1,000
        # topology      Shape (377,2)   Type i8     Array [[0 1] [1 2] ... [n-1 n]]      incremented by 1, shows which points are attched to each other


        #########################################################
        #################  FUNCTION SPACES  #####################
        #########################################################
        self.E_Q = FiniteElement("CG",self.mesh.ufl_cell(),1)
        self.Q   = FunctionSpace(self.mesh, self.E_Q) # E_Q defined over the mesh
        self.E_V = MixedElement(self.E_Q, self.E_Q, self.E_Q)
        self.V   = FunctionSpace(self.mesh ,self.E_V)


        # For moving data between vector functions and scalar functions
        self.assigner_inv = FunctionAssigner([self.Q,self.Q,self.Q],self.V)
        self.assigner     = FunctionAssigner(self.V, [self.Q, self.Q, self.Q])

        # For solution
        self.U   = Function(self.V)
        self.dU  = TrialFunction(self.V)
        self.Phi = TestFunction(self.V)
        self.u,     self.u2, self.H   = split(self.U) # H will be the thickness at current time
        self.phi, self.phi1, self.xsi = split(self.Phi) # Individual test functions needed for forms

        # Placeholders needed for assignment
        self.un  = Function(self.Q)
        self.u2n = Function(self.Q)


        # Zeros, an initial guess for velocity for when solver fails
        self.zero_sol   = Function(self.Q)
        #zero_sol.assign(Constant(0))

        #########################################################
        #############      FIELD INITIALIZATION   ###############
        #########################################################
        self.S0 = Function(self.Q) # Initial surface elevation
        self.B  = Function(self.Q) # Bed elevation
        self.H0 = Function(self.Q) # Thickness at previous time step
        self.A  = Function(self.Q) # SMB data

        # in_file.read(H0.vector(), "/thickness", True) #NEW
        self.in_file.read(self.S0.vector(), "/surface",   True)
        self.in_file.read(self.B.vector(),  "/bed",       True)
        self.in_file.read(self.A.vector(),  "/smb",       True)
        self.H0.assign(self.S0-self.B)   # Initial thickness  #OLD

        # A generalization of the Crank-Nicolson method, which is theta = .5
        self.Hmid = theta*self.H + (1-theta)*self.H0

        # Define surface elevation
        self.S = self.B + self.Hmid  #OLD
        # S = Max(Hmid+B, rho/rho_w * Hmid)
        # Expressions for now, later move to data from files
        # adot = interpolate(Adot(degree=1),Q)
        # dolfin.fem.interpolation.interpolate: adot is an expression, Q is a functionspace
            # returns an interpolation of a given function into a given finite element space
            # i think this is a vector
        # print 'adot: ', adot.vector()[:][0]
        # Q is a functionspace
        # print 'adot: ', adot  #printed 'f_51
        self.width = interpolate(Width(degree=2), self.Q)

        #############################################################################
        #######################  MOMENTUM CONSERVATION  #############################
        #############################################################################
        # This object stores the stresses
        self.strs = Stresses(self.U, self.Hmid, self.H0, self.H, self.width, self.B, self.S, self.Phi)
        # Conservation of momentum form:
        self.R = -(self.strs.tau_xx + self.strs.tau_xz + self.strs.tau_b + self.strs.tau_d + self.strs.tau_xy)*dx

        #############################################################################
        ########################  MASS CONSERVATION  ################################
        #############################################################################
        self.h = CellSize(self.mesh)
        self.D = self.h*abs(self.U[0])/2.
        self.area = self.Hmid*self.width

        self.mesh_min = self.mesh.coordinates().min()
        self.mesh_max = self.mesh.coordinates().max()

        # Define boundaries
        self.ocean = FacetFunctionSizet(self.mesh,0)
        self.ds = fc.ds(subdomain_data=self.ocean) #THIS DS IS FROM FENICS! border integral

        for f in facets(self.mesh):
            if near(f.midpoint().x(),self.mesh_max):
               self.ocean[f] = 1
            if near(f.midpoint().x(),self.mesh_min):
               self.ocean[f] = 2

        # Directly write the form, with SPUG and area correction,
        self.R += ((self.H-self.H0)/dt*self.xsi - self.xsi.dx(0)*self.U[0]*self.Hmid + self.D*self.xsi.dx(0)*self.Hmid.dx(0) - (self.A - self.U[0]*self.H/self.width*self.width.dx(0))*self.xsi)*dx\
               + self.U[0]*self.area*self.xsi*self.ds(1) - self.U[0]*self.area*self.xsi*self.ds(0)


        #####################################################################
        #########################  SOLVER SETUP   ###########################
        #####################################################################

        # Bounds
        self.l_thick_bound = project(Constant(thklim),self.Q)
        self.u_thick_bound = project(Constant(1e4),self.Q)

        self.l_v_bound = project(-10000.0,self.Q)
        self.u_v_bound = project(10000.0,self.Q)


        self.l_bound = Function(self.V)
        self.u_bound = Function(self.V)

        self.assigner.assign(self.l_bound,[self.l_v_bound]*2+[self.l_thick_bound])
        self.assigner.assign(self.u_bound,[self.u_v_bound]*2+[self.u_thick_bound])

        # This should set the velocity at the divide (left) to zero
        self.dbc0 = DirichletBC(self.V.sub(0),0,lambda x,o:near(x[0],self.mesh_min) and o)
        # Set the velocity on the right terminus to zero
        self.dbc1 = DirichletBC(self.V.sub(0),0,lambda x,o:near(x[0],self.mesh_max) and o)
        # overkill?
        self.dbc2 = DirichletBC(self.V.sub(1),0,lambda x,o:near(x[0],self.mesh_max) and o)
        # set the thickness on the right edge to thklim
        self.dbc3 = DirichletBC(self.V.sub(2),thklim,lambda x,o:near(x[0],self.mesh_max) and o)

        #Define variational solver for the mass-momentum coupled problem
        self.J = derivative(self.R,self.U,self.dU)

        self.coupled_problem = NonlinearVariationalProblem(self.R,self.U,bcs=[self.dbc0,self.dbc1,self.dbc3],J=self.J)

        self.coupled_problem.set_bounds(self.l_bound, self.u_bound)

        self.coupled_solver = NonlinearVariationalSolver(self.coupled_problem)

        # Accquire the optimizations in fenics_optimizations
        set_solver_options(self.coupled_solver)

        ######################################################################
        #######################   TIME LOOP   ################################
        ######################################################################

        # Time interval
        self.t = 0
        # t_end = 20000.
        self.t_end = float(self.tEndLineEdit.text())
        self.dt_float = float(self.tStepLineEdit.text())
        self.pPlt = pyqtplotter(self.strs, self.mesh, self.plt1, self.plt2, self.plt3, self.t, self.dt_float, str(self.saveFileLineEdit.text()))
        self.pPlt.refresh_plot(0)
        pg.QtGui.QApplication.processEvents()
        self.in_file.close()
        self.runLoop()


    def runLoop(self):
        while self.t < self.t_end and self.pPlt.run:
            # time0 = time.time()
            print( "Solving for time: ", self.t)
            self.tCurrent.setText("Current year: " + str(self.t))
            self.coupled_problem = NonlinearVariationalProblem(self.R, self.U, bcs=[self.dbc0, self.dbc1, self.dbc3], J=self.J)
            self.coupled_problem.set_bounds(self.l_bound, self.u_bound)
            self.coupled_solver = NonlinearVariationalSolver(self.coupled_problem)

            # Accquire the optimizations in fenics_optimizations
            set_solver_options(self.coupled_solver)

            try:
                self.coupled_solver.solve(set_solver_options())
            except:
                print ("Exception Triggered!")
                self.coupled_solver.parameters['snes_solver']['error_on_nonconvergence'] = False
                self.assigner.assign(self.U,[self.zero_sol,self.zero_sol,self.H0])
                self.coupled_solver.solve()
                self.coupled_solver.parameters['snes_solver']['error_on_nonconvergence'] = True

            self.assigner_inv.assign([self.un,self.u2n,self.H0],self.U)
            self.t += self.dt_float
            self.pPlt.refresh_plot(self.t)
            pg.QtGui.QApplication.processEvents()



        




