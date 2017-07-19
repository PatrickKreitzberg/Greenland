# NOTE: This code, pasted into profilev2, will work...
# But, it'd be nicer if there were a seperate class file!

from physical_constants    import * 
from simulation_parameters import * 
from numpy import array
from dolfin import Constant,Max,sqrt

class VerticalBasis(object):
    def __init__(self,u,coef,dcoef):
        self.u = u
        self.coef = coef
        self.dcoef = dcoef

    def __call__(self,s):
        return sum([u*c(s) for u,c in zip(self.u,self.coef)])

    def ds(self,s):
        return sum([u*c(s) for u,c in zip(self.u,self.dcoef)])

    def dx(self,s,x):
        return sum([u.dx(x)*c(s) for u,c in zip(self.u,self.coef)])

class VerticalIntegrator(object):
    def __init__(self,points,weights):
        self.points = points
        self.weights = weights
    def integral_term(self,f,s,w):
        return w*f(s)
    def intz(self,f):
        return sum([self.integral_term(f,s,w) for s,w in zip(self.points,self.weights)])


# What follows should rightly be wrapped in a class, but it doesn't work (see momentum_class.py)
# For now, just use the following definitions (which is a little sloppy in terms of name spaces)

coef  = [lambda s:1.0, lambda s:1./4.*(5*s**4-1.)] 
dcoef = [lambda s:0.0, lambda s:5*s**3]           # Derivatives of above forms

# This is a quadrature rule for vertical integration
points  = array([0.0,0.4688,0.8302,1.0])
weights = array([0.4876/2.,0.4317,0.2768,0.0476])

# This object does the integration
vi = VerticalIntegrator(points,weights)

u    = VerticalBasis([U[0]  ,  U[1]],coef,dcoef)
phi  = VerticalBasis([Phi[0],Phi[1]],coef,dcoef)

# These last two do not require vertical integration
# Note the impact this has on the test functions associated with them.
def _tau_xy():
    return 2.*Hmid*b/width*((n+2)/(width))**(1./n)*(abs(U[0])+1.0)**(1./n - 1)*U[0]

def _tau_b():
    P_i = rho*g*Hmid                 # Overburden, or hydrostatic pressure of ice
    P_w = Max(-rho_w*g*B,thklim)     # Water pressure is either ocean pressure or zero
    N   = Max((P_i - P_w)/(rho*g),thklim) # Effective pressure is P_i - P_w, rho*g appears in A_s
    normalx = (B.dx(0))/sqrt((B.dx(0)**2 + 1.0))
    normalz = sqrt(1 - normalx**2)
    return mu*A_s*(abs(u(1))+1.0)**(1./n-1)*u(1)*(abs(N)+1.0)**(1./m)*(1-normalx**2)

# Sigma coordinate derivative transforms
def _dsdx(s):
    return 1./Hmid*(S.dx(0) - s*H0.dx(0))  # was H on second term!

def _dsdz(s):
    return -1./Hmid

def _eta_v(s):
    return b/2.*((u.dx(s,0) + u.ds(s)*_dsdx(s))**2 \
                +0.25*((u.ds(s)*_dsdz(s))**2) \
                + eps_reg)**((1.-n)/(2*n))

def _tau_xx(s):
    return (phi.dx(s,0) + phi.ds(s)*_dsdx(s))*Hmid*_eta_v(s)\
    *(4*(u.dx(s,0) + u.ds(s)*_dsdx(s))) + (phi.dx(s,0) \
    + phi.ds(s)*_dsdx(s))*H0*_eta_v(s)*\
    (2*u(s)/width*width.dx(0))       # Final Hmid was H0

def _tau_xz(s):
    return _dsdz(s)**2*phi.ds(s)*Hmid*_eta_v(s)*u.ds(s)
def _tau_d(s):
    return rho*g*Hmid*S.dx(0)*phi(s)

# The stresses ready to be accessed
tau_d       = vi.intz(_tau_d)
tau_xx      = vi.intz(_tau_xx)
tau_xz      = vi.intz(_tau_xz)
tau_xy      = Phi[0]      * _tau_xy()
tau_b       = phi(1) * _tau_b()

# Forms of stresses ammenable to plotting (no test functions)
tau_b_plot  = _tau_b()
tau_d_plot  = rho*g*Hmid*S.dx(0)
tau_xx_plot = vi.intz(lambda s: _dsdx(s)*Hmid*_eta_v(s)*\
                          (4*(u.dx(s,0) + u.ds(s)*_dsdx(s)))\
                           +_dsdx(s)*H0*_eta_v(s)*\
                           (2*u(s)/width*width.dx(0)))
tau_xz_plot = vi.intz(lambda s:_dsdz(s)**2*Hmid*_eta_v(s)*u.ds(s))
tau_xy_plot = _tau_xy()
