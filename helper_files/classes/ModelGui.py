import pyqtgraph as pg
from pyqtgraph.Qt import QtGui
import fenics as fc
import h5py
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
    def __init__(self, parent, hdf_name):
        self.hdf_name = None
        self.parent = parent
        QtGui.QMainWindow.__init__(self, self.parent)
        self.cw = QtGui.QWidget()
        self.setCentralWidget(self.cw)
        self.horLayout = QtGui.QHBoxLayout()
        self.cw.setLayout(self.horLayout)


        # RIGHT PANEL

        self.rightPanelW = QtGui.QWidget()
        self.rightPanelLay = QtGui.QGridLayout()
        self.rightPanelW.setLayout(self.rightPanelLay)

            # GUI TIME COMPONENTS
        self.timeWidget = QtGui.QWidget()
        self.timeContainer = QtGui.QGridLayout()
        self.timeWidget.setLayout(self.timeContainer)
        self.tEndLabel = QtGui.QLabel('t_end(yr):')
        self.tEndLineEdit = QtGui.QLineEdit('20000')
        self.tStepLabel = QtGui.QLabel('t_step(yr):')
        self.tStepLineEdit = QtGui.QLineEdit('10')
        self.tCurrent = QtGui.QLabel('Current year: ')
        self.timeContainer.addWidget(self.tEndLabel, 1, 0)
        self.timeContainer.addWidget(self.tEndLineEdit, 1, 1)
        self.timeContainer.addWidget(self.tStepLabel, 2, 0)
        self.timeContainer.addWidget(self.tStepLineEdit, 2, 1)
        self.timeContainer.addWidget(self.tCurrent, 3, 0, 1, 2)

            # SAVE FILE
        self.saveFileW = QtGui.QWidget()
        self.saveFileLay = QtGui.QHBoxLayout()
        self.saveFileLabel = QtGui.QLabel('Out File')
        self.saveFileLineEdit = QtGui.QLineEdit('./data/outModel.h5')
        self.saveFileLay.addWidget(self.saveFileLabel)
        self.saveFileLay.addWidget(self.saveFileLineEdit)

            # BUTTONS
        self.runButt = QtGui.QPushButton('Run Model')
        self.pauseButt = QtGui.QPushButton('Pause Model')

            # ADD TO LAYOUT

        self.rightPanelLay.addWidget(self.timeContainer)
        self.rightPanelLay.addWidget(self.saveFileW)
        self.rightPanelLay.addWidget(self.runButt)
        self.rightPanelLay.addWidget(self.pauseButt)
        self.pauseButt.setEnabled(False)

        # LEFT SIDE
        self.leftPanelW = QtGui.QWidget()
        self.leftPanelLay = QtGui.QVBoxLayout()  # CENTRAL WIDGET LAYOUT (layout of the central widget)
        self.leftPanelW.setLayout(self.leftPanelLay)
        self.plt1 = pg.PlotWidget()
        self.plt2 = pg.PlotWidget()
        self.plt3 = pg.PlotWidget()
        self.slider = QtGui.QSlider()
        self.leftPanelLay.addWidget(self.plt1)
        self.leftPanelLay.addWidget(self.plt2)
        self.leftPanelLay.addWidget(self.plt3)
        self.leftPanelLay.addWidget(self.slider)

        # ADD SIDES TO MAIN WINDOW

        self.horLayout.addWidget(self.leftPanelW)
        self.horLayout.addWidget(self.rightPanelW)

    def runModelButt(self):
        dr = float(model_res_lineEdit.text())  # dr = 150
        if len(vpts) > 0:
            interpolateData(True)

            ###########################################
            ###  INTERPOLATE DATA SO EVEN INTERVAL  ###
            ###########################################
            '''
            Interpolating data at an even interval so the data points align with the 
            invterval mesh.
            '''
            thickness1dInterp = interp1d(thickness.distanceData, thickness.pathData)
            bed1dInterp       = interp1d(bed.distanceData, bed.pathData)
            surface1dInterp   = interp1d(surface.distanceData, surface.pathData)
            smb1dInterp       = interp1d(smb.distanceData, smb.pathData)
            velocity1dInterp  = interp1d(velocity.distanceData, velocity.pathData)

            # N is the number of total data points including the last
            # Data points on interval [0, N*dr] inclusive on both ends

            N = int(np.floor(bed.distanceData[-1] / float(dr)))  # length of path / resolution

            x = np.arange(0, (N + 1) * dr, dr)  # start point, end point, number of segments. END POINT NOT INCLUDED!
            print 'N, dr, dr*N', N, dr, dr * N
            mesh = fc.IntervalMesh(N, 0, dr * N)  # number of cells, start point, end point

            thicknessModelData = thickness1dInterp(x)
            bedModelData = bed1dInterp(x)
            surfaceModelData = surface1dInterp(x)
            smbModelData = smb1dInterp(x)
            velocityModelData = velocity1dInterp(x)

            THICKLIMIT = 10.  # Ice is never less than this thick
            H = surfaceModelData - bedModelData
            surfaceModelData[H <= THICKLIMIT] = bedModelData[H <= THICKLIMIT]

            # FIXME the intervalMesh is consistantly 150 between each datapoint this not true for the data being sent
            self.hdf_name = '.data/latest_profile.h5'
            hfile = fc.HDF5File(mesh.mpi_comm(), self.hdf_name, "w")
            V = fc.FunctionSpace(mesh, "CG", 1)

            functThickness = fc.Function(V, name="Thickness")
            functBed = fc.Function(V, name="Bed")
            functSurface = fc.Function(V, name="Surface")
            functSMB = fc.Function(V, name='SMB')
            functVelocity = fc.Function(V, name='Velocity')

            surface.pathPlotItem.setData(x, surfaceModelData)
            pg.QtGui.QApplication.processEvents()

            functThickness.vector()[:] = thicknessModelData
            functBed.vector()[:] = bedModelData
            functSurface.vector()[:] = surfaceModelData
            functSMB.vector()[:] = smbModelData
            functVelocity.vector()[:] = velocityModelData

            hfile.write(functThickness.vector(), "/thickness")
            hfile.write(functBed.vector(), "/bed")
            hfile.write(functSurface.vector(), "/surface")
            hfile.write(functSMB.vector(), "/smb")
            hfile.write(functVelocity.vector(), "/velocity")
            hfile.write(mesh, "/mesh")
            hfile.close()
            self.runModel()

    def runModel(self):

        ##########################################################
        ################           FILES         #################
        ##########################################################

        #FIXME Probably dont have to save then open mesh
        self.mesh = Mesh()
        self.in_file  = HDF5File(mesh.mpi_comm(), self.hdf_name, "r")  #mesh.mpi_comm() is ussed to read in parallel?
        # outF = h5py.File('./data/modelOut')


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
        self.in_file.read(mesh,"/mesh", False)


        # H5FILE Data:
        # bed and surface are datasets shape (378,)
        # mesh is a group with datasets
        # cell_indices  Shape (377,)    Type i8     Array [0 1 2 ... n-1 n]              incremented by 1
        # coordinates   Shape (378,1)   Type f8     Array [[0.] [1000.] ... [ 377000.]]  incremented by 1,000
        # topology      Shape (377,2)   Type i8     Array [[0 1] [1 2] ... [n-1 n]]      incremented by 1, shows which points are attched to each other


        #########################################################
        #################  FUNCTION SPACES  #####################
        #########################################################
        self.E_Q = FiniteElement("CG",mesh.ufl_cell(),1)
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
        self.Hmid = self.theta*self.H + (1-self.theta)*self.H0

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
        self.h = CellSize(mesh)
        self.D = self.h*abs(self.U[0])/2.
        self.area = self.Hmid*self.width

        self.mesh_min = self.mesh.coordinates().min()
        self.mesh_max = self.mesh.coordinates().max()

        # Define boundaries
        self.ocean = FacetFunctionSizet(self.mesh,0)
        self.ds = fc.ds(subdomain_data=self.ocean) #THIS DS IS FROM FENICS! border integral

        for f in facets(self.mesh):
            if near(f.midpoint().x(),self.mesh_max):
               ocean[f] = 1
            if near(f.midpoint().x(),self.mesh_min):
               ocean[f] = 2

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

        self.coupled_problem = NonlinearVariationalProblem(self.R,self.U,bcs=[dbc0,dbc1,dbc3],J=J)

        self.coupled_problem.set_bounds(self.l_bound, self.u_bound)

        self.coupled_solver = NonlinearVariationalSolver(self.coupled_problem)

        # Accquire the optimizations in fenics_optimizations
        self.set_solver_options(self.coupled_solver)

        ######################################################################
        #######################   TIME LOOP   ################################
        ######################################################################

        # Time interval
        t = 0
        # t_end = 20000.
        t_end = float(t_end_lineEdit.text())
        dt_float = float(t_step_lineEdit.text())

        # PyQt gui items
        mw2 = QtGui.QMainWindow(mw)
        mw2.setWindowTitle('PyQt PLOTTER')  # MAIN WINDOW
        cw2 = QtGui.QWidget()  # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
        mw2.setCentralWidget(cw2)
        l = QtGui.QVBoxLayout()  # CENTRAL WIDGET LAYOUT (layout of the central widget)
        cw2.setLayout(l)
        plt1 = pg.PlotWidget()
        plt2 = pg.PlotWidget()
        plt3 = pg.PlotWidget()
        l.addWidget(plt1)
        l.addWidget(plt2)
        l.addWidget(plt3)
        mw2.show()

        pPlt = pyqtplotter(strs, mesh, plt1, plt2, plt3, t, dt_float)
        pPlt.refresh_plot(0)
        mw2.closeEvent = pPlt.closed
        pg.QtGui.QApplication.processEvents()
        print 'mw2..isActive', mw2.isActiveWindow()
        mw2.activateWindow()
        print 'mw2..isActive', mw2.isActiveWindow()

        while t<t_end and pPlt.run:
            # time0 = time.time()
            print( "Solving for time: ",t)
            t_current.setText("Current year: " + str(t))
            coupled_problem = NonlinearVariationalProblem(R, U, bcs=[dbc0, dbc1, dbc3], J=J)
            coupled_problem.set_bounds(l_bound, u_bound)
            coupled_solver = NonlinearVariationalSolver(coupled_problem)

            # Accquire the optimizations in fenics_optimizations
            set_solver_options(coupled_solver)

            try:
                coupled_solver.solve(set_solver_options())
            except:
                print ("Exception Triggered!")
                coupled_solver.parameters['snes_solver']['error_on_nonconvergence'] = False
                assigner.assign(U,[zero_sol,zero_sol,H0])
                coupled_solver.solve()
                coupled_solver.parameters['snes_solver']['error_on_nonconvergence'] = True

            assigner_inv.assign([un,u2n,H0],U)
            t += dt_float
            pPlt.refresh_plot(t)
            pg.QtGui.QApplication.processEvents()
        in_file.close()



        




