#!/usr/bin/env python
"""
Balance velocity optimization script
Including facilities to perform an optimization

Modified on Friday Aug 7 2015
@author: jessej
"""
from dolfin import *
from dolfin.cpp.function import *
from dolfin.cpp.common import *
from dolfin.cpp.io import *
from dolfin.cpp.fem import *
from dolfin.cpp.mesh import *
import fenics as fc
import numpy as np
import scipy
import pylab
from smooth import smooth
from scipy.optimize import fmin_l_bfgs_b

#set_log_level(PROGRESS)

# Key paths and files
output_path = "./model_output_fenics_only/totten_output_3H/"
data_file   = "mesh2d.h5"
mesh        = Mesh('mesh2d.xml')

# Smoothing paramters. This many ice thicknesses are averaged over
# to get the directions from surface slope.
kappa_i = 12. # In regions where there are no data
kappa_v = 1.  # Smooth the data directions too

# This is the speed where the directional data is 
# more trustworthy. Low speed often has strange directions...
goodU = 5. 

# Regularization parameter, penalty for gradients in adot
alpha = 1e9

# Bounds on adot, which is our control variable
adot_err = 2.0 # plus or minus this amount

# Function spaces assigned prior to data read.
Q = fc.FunctionSpace(mesh,"CG",1) # First order continuous galerkin
H    = Function(Q) # Thickness
S    = Function(Q) # Surface elevation
B    = Function(Q) # Bed elevation
adot = Function(Q) # Surface mass balance
u    = Function(Q) # x-horizontal component of velocity
v    = Function(Q) # y-horizontal component of velocity

# Import data
f = HDF5File(mesh.mpi_comm(),data_file,'r')
print f.has_dataset('bed')
f.read(B,        'bed')
f.read(S,        'surface')
f.read(H,        'thickness')
f.read(adot,     'smb')
# Read previously saved adot to start from previous simulation's results
# File('./totten_2H_data/adot_current.xml') >> adot
# f.read(u,        'u')
f.read(v,        'velocity')

# Simple data pruning
thklim = 10.
H0 = H.copy()
H0.vector()[H0.vector()<thklim] = thklim 
H.vector()[H.vector()  <thklim] = thklim

# Convenience
dSdx = project(.1,Q)
dSdy = project(.1,Q)
Unorm = v

# Mark facets where velocity is non-zero, 
# Objective function is only computed in these regions
print mesh.coordinates()
print mesh.cells()
class ValidData(SubDomain):
    def __init__(self,u):
        SubDomain.__init__(self)
        self.u = u

    def inside(self,x,on_boundary):
        print 'x,y', x[0], x[1]
        u = self.u(x[0], x[1])
        return u != 0.

subdomains = CellFunction("size_t", mesh)
subdomains.set_all(0)
validdata = ValidData(u)
validdata.mark(subdomains,1)
dx = Measure("dx")[subdomains]

# Smooth the surface and find flow directions based on
# Smoothed surface gradients
rho = 911
g   = 9.81
Nx = smooth(Q, H, -rho * g * H * dSdx, kappa_i)
Ny = smooth(Q, H, -rho * g * H * dSdy, kappa_i)

# Replace smoothed surface slope direction with data where ice is fast enough
# to lower errors:
uhat = project(u / Unorm,Q)
vhat = project(v / Unorm,Q)
# Smooth these directions out too, but over fewer ice thicknesses
uhat = smooth(Q, H, rho * g * H * uhat, kappa_v)
vhat = smooth(Q, H, rho * g * H * vhat, kappa_v)

Nx.vector()[Unorm.vector() > goodU] = uhat.vector()[Unorm.vector() > goodU]
Ny.vector()[Unorm.vector() > goodU] = vhat.vector()[Unorm.vector() > goodU]

# Finally, smooth again to easy the transitions between data and surface slope
Nx = smooth(Q, H, rho * g * H * Nx, kappa_v)
Ny = smooth(Q, H, rho * g * H * Ny, kappa_v)
Nnorm=project(sqrt(Nx**2+Ny**2)     + 1e-16,Q)
dS = as_vector([project(Nx / Nnorm, Q), project(Ny / Nnorm, Q)])

# Boundary conditions
def boundary(x,on_boundary):
  return on_boundary

dbc_U     = DirichletBC(Q,0.,boundary)
dbc_adj   = DirichletBC(Q,0.0,boundary)

# SPUG test functions
h = CellSize(mesh)
U_eff = H0
tau = h / (2.0*U_eff)

# Forward function spaces
Ub  = Function(Q)
dUb = TestFunction(Q)
tUb = TrialFunction(Q)
# Adjoint function spaces
lamda  = Function(Q)
dlamda = TestFunction(Q)
tlamda = TrialFunction(Q)
da     = TestFunction(Q)


# Lagrangian
I   = 0.5*ln(sqrt(Ub**2+1)/sqrt(Unorm**2+1))**2*dx(1) + \
      alpha*dot(grad(adot),grad(adot))*dx + \
      (lamda + tau*div(H0*dS*lamda))*(div(Ub*dS*H0) - adot)*dx

# Unconstrainted lagrangian
I_uc= 0.5*ln(sqrt(Ub**2+1)/sqrt(Unorm**2+1))**2*dx(1) + \
      alpha*dot(grad(adot),grad(adot))*dx 

# Derivitive of objective function
dI  = derivative(ln(sqrt(Ub**2+1)/sqrt(Unorm**2+1))**2,Ub,dUb)

# Forward continuity model
fwd = (dlamda + tau*div(H0*dS*dlamda))*(div(H0*dS*Ub) - adot)*dx
# Adjoint model
adj = dI*dx(1) + (lamda + tau*div(dS*H0*lamda))*div(dS*H0*dUb)*dx
gradient = alpha*dot(grad(adot),grad(da))*dx - (lamda + tau*div(dS*H0*lamda))*da*dx

J_fwd = derivative(fwd,Ub,tUb) 
J_adj = derivative(adj,lamda,tlamda) 

# Setup output
Ufile    = File(output_path+"U.pvd")
adotfile = File(output_path+"adot.pvd")
U        = Function(VectorFunctionSpace(mesh,'CG',1,dim=2))
adot_    = Function(Q)

# The _I_fun computes the value of the objective function from doing a forward
# solve.

def _I_fun(s,*args):
    adot.vector()[:] = s
    solve(fwd == 0, Ub, dbc_U, J = J_fwd,\
          solver_parameters={"newton_solver":{"maximum_iterations": 100,\
              "error_on_nonconvergence": False}},\
          form_compiler_parameters={"optimize": True})
    # Output to files
    Ut = project(as_vector([Ub*dS[0],Ub*dS[1]]),VectorFunctionSpace(mesh,'CG',1,dim=2))
    U.interpolate(Ut)
    Ufile << U
    File(output_path+"Ub_current.xml")<< U

    adot_.assign(adot)
    adotfile << adot_
    File(output_path+"adot_current.xml")<< adot_
    return assemble(I_uc)

# The _J_fun computes the slope of the objective function with respect to the
# control parameter by solving the adjoint system
def _J_fun(s,*args):
    adot.vector()[:] = s
    solve(adj == 0, lamda, dbc_adj, J = J_adj,\
          solver_parameters={"newton_solver":{"maximum_iterations": 100,\
              "error_on_nonconvergence": False}},\
          form_compiler_parameters={"optimize": True})
    g = assemble(gradient) 
    return g.array()

# Set some bounds on adot in the optimization, else it goes way off
adot_bounds = [(aa - adot_err, aa + adot_err) for aa in adot.vector().array()]

# Use BFGS to solve the optimization problem.
# Use maxfun = 1200 or more for really nice results.
aopt,f,d = fmin_l_bfgs_b(_I_fun,adot.vector().array(),fprime=_J_fun,\
        iprint=1,bounds=adot_bounds, maxfun = 1200)

# Copy the optimized balance velocity into the blank spots in the observed velocity
ub = project(Ub * Nx / Nnorm, Q)
u.vector()[u.vector() == 0] = ub.vector()[u.vector() == 0]

vb = project(Ub * Ny / Nnorm, Q)
v.vector()[v.vector() == 0] = vb.vector()[v.vector() == 0]

Uf  = project(as_vector([u,v]),VectorFunctionSpace(mesh,'CG',1,dim=2))
Ubf = project(as_vector([Ub*dS[0],Ub*dS[1]]),VectorFunctionSpace(mesh,'CG',1,dim=2))

File(output_path+"Ufinal.pvd")  << Uf
File(output_path+"Ubfinal.pvd") << Ubf
File(output_path+"Ubfinal.xml") << Ubf
File(output_path+"Ufinal.xml")  << Uf
