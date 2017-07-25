import fenics as fc
from time import *
from dolfin import *

from gui import *
from helper_files.pyQtPlotter import *
from support.expressions import *
from support.fenics_optimizations import *
from support.momentum import *

run = True

def winclose(e):
    global run
    run = False


def runModel(hdf_name):
    global run
    run = True

    ##########################################################
    ################           FILES         #################
    ##########################################################

    mesh = Mesh()
    in_file  = HDF5File(mesh.mpi_comm(), hdf_name, "r")  #mesh.mpi_comm() is ussed to read in parallel?


    # out_file = HDF5File(mesh.mpi_comm(),"./output_data/peterman.h5","w")

    # cell_indices: self-explanatory
    # coordinates:  self-explanatory
    # topology:     Shows which nodes are linked together

    ##########################################################
    ################           GUI           #################
    ##########################################################


    # PyQt gui items
    mw2 = QtGui.QMainWindow()
    mw2.setWindowTitle('PyQt PLOTTER')  # MAIN WINDOW
    mw2.closeEvent = winclose
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

    ##########################################################
    ################           MESH          #################
    ##########################################################
    in_file.read(mesh,"/mesh", False)

    # H5FILE Data:
    # bed and surface are datasets shape (378,)
    # mesh is a group with datasets
    # cell_indices  Shape (377,)    Type i8     Array [0 1 2 ... n-1 n]              incremented by 1
    # coordinates   Shape (378,1)   Type f8     Array [[0.] [1000.] ... [ 377000.]]  incremented by 1,000
    # topology      Shape (377,2)   Type i8     Array [[0 1] [1 2] ... [n-1 n]]      incremented by 1, shows which points are attched to each other


    #########################################################
    #################  FUNCTION SPACES  #####################
    #########################################################
    E_Q = FiniteElement("CG",mesh.ufl_cell(),1)
    Q   = FunctionSpace(mesh,E_Q) # E_Q defined over the mesh
    E_V = MixedElement(E_Q,E_Q,E_Q)
    V   = FunctionSpace(mesh,E_V)


    # For moving data between vector functions and scalar functions
    assigner_inv = FunctionAssigner([Q,Q,Q],V)
    assigner     = FunctionAssigner(V,[Q,Q,Q])

    # For solution
    U   = Function(V)
    dU  = TrialFunction(V)
    Phi = TestFunction(V)
    u,     u2, H   = split(U) # H will be the thickness at current time
    phi, phi1, xsi = split(Phi) # Individual test functions needed for forms

    # Placeholders needed for assignment
    un  = Function(Q)
    u2n = Function(Q)


    # Zeros, an initial guess for velocity for when solver fails
    zero_sol   = Function(Q)
    #zero_sol.assign(Constant(0))

    #########################################################
    #############      FIELD INITIALIZATION   ###############
    #########################################################
    S0 = Function(Q) # Initial surface elevation
    B  = Function(Q) # Bed elevation
    H0 = Function(Q) # Thickness at previous time step
    A  = Function(Q) # SMB data

    in_file.read(S0.vector(), "/surface", True)
    in_file.read(B.vector(),  "/bed",     True)
    in_file.read(A.vector(),  "/smb",     True)
    H0.assign(S0-B)   # Initial thickness

    # A generalization of the Crank-Nicolson method, which is theta = .5
    Hmid = theta*H + (1-theta)*H0

    # Define surface elevation
    S = B + Hmid

    # Expressions for now, later move to data from files
    # adot = interpolate(Adot(degree=1),Q)
    # dolfin.fem.interpolation.interpolate: adot is an expression, Q is a functionspace
        # returns an interpolation of a given function into a given finite element space
        # i think this is a vector
    # print 'adot: ', adot.vector()[:][0]
    # Q is a functionspace
    # print 'adot: ', adot  #printed 'f_51
    width = interpolate(Width(degree=2),Q)

    #############################################################################
    #######################  MOMENTUM CONSERVATION  #############################
    #############################################################################
    # This object stores the stresses
    strs = Stresses(U,Hmid,H0,H,width,B,S,Phi)
    # Conservation of momentum form:
    R = -(strs.tau_xx + strs.tau_xz + strs.tau_b + strs.tau_d + strs.tau_xy)*dx

    #############################################################################
    ########################  MASS CONSERVATION  ################################
    #############################################################################
    h = CellSize(mesh)
    D = h*abs(U[0])/2.
    area = Hmid*width

    mesh_min = mesh.coordinates().min()
    mesh_max = mesh.coordinates().max()

    # Define boundaries
    ocean = FacetFunctionSizet(mesh,0)
    ds = fc.ds(subdomain_data=ocean) #THIS DS IS FROM FENICS! border integral

    for f in facets(mesh):
        if near(f.midpoint().x(),mesh_max):
           ocean[f] = 1
        if near(f.midpoint().x(),mesh_min):
           ocean[f] = 2

    # Directly write the form, with SPUG and area correction,
    R += ((H-H0)/dt*xsi - xsi.dx(0)*U[0]*Hmid + D*xsi.dx(0)*Hmid.dx(0) - (A - U[0]*H/width*width.dx(0))*xsi)*dx\
           + U[0]*area*xsi*ds(1) - U[0]*area*xsi*ds(0)
    print 'smb: ', A.vector()[:]


    #####################################################################
    #########################  SOLVER SETUP   ###########################
    #####################################################################

    # Bounds
    l_thick_bound = project(Constant(thklim),Q)
    u_thick_bound = project(Constant(1e4),Q)

    l_v_bound = project(-10000.0,Q)
    u_v_bound = project(10000.0,Q)


    l_bound = Function(V)
    u_bound = Function(V)

    assigner.assign(l_bound,[l_v_bound]*2+[l_thick_bound])
    assigner.assign(u_bound,[u_v_bound]*2+[u_thick_bound])

    # This should set the velocity at the divide (left) to zero
    dbc0 = DirichletBC(V.sub(0),0,lambda x,o:near(x[0],mesh_min) and o)
    # Set the velocity on the right terminus to zero
    dbc1 = DirichletBC(V.sub(0),0,lambda x,o:near(x[0],mesh_max) and o)
    # overkill?
    dbc2 = DirichletBC(V.sub(1),0,lambda x,o:near(x[0],mesh_max) and o)
    # set the thickness on the right edge to thklim
    dbc3 = DirichletBC(V.sub(2),thklim,lambda x,o:near(x[0],mesh_max) and o)

    #Define variational solver for the mass-momentum coupled problem
    J = derivative(R,U,dU)

    coupled_problem = NonlinearVariationalProblem(R,U,bcs=[dbc0,dbc1,dbc3],J=J)

    coupled_problem.set_bounds(l_bound,u_bound)

    coupled_solver = NonlinearVariationalSolver(coupled_problem)

    # Accquire the optimizations in fenics_optimizations
    set_solver_options(coupled_solver)

    ######################################################################
    #######################   TIME LOOP   ################################
    ######################################################################

    # Time interval
    t = 0
    # t_end = 20000.
    t_end = float(t_end_lineEdit.text())
    dt_float = float(t_step_lineEdit.text())
    pPlt = pyqtplotter(strs, mesh, plt1, plt2, plt3)
    pPlt.refresh_plot()
    pg.QtGui.QApplication.processEvents()

    while t<t_end and run:
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
        t+=dt_float
        pPlt.refresh_plot()
        pg.QtGui.QApplication.processEvents()
    print 'done 1'
    in_file.close()
    # pg.exit()
    print 'done 2'
    # mw2.at



# m = model()
## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
