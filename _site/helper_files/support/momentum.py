from physical_constants    import * 
from simulation_parameters import * 
from numpy import array
from dolfin import Constant,Max,sqrt
########################################################
#################   SUPPORT FUNCTIONS  #################
########################################################
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

class Stresses(object):
    """ 
    The purpose of this class is to return the vertically integrated stresses.
    This is achived by accessing one of:
        tau_xx - longitudinal
        tau_xy - lateral drag
        tau_xz - vertical shear
        tau_d  - driving stress
        tau_b  - basal drag
    """


    def __init__(self,U,Hmid,H0,H,width,B,S,Phi):
        # These functions represent our prior about the vertical velocity profiles
        # Functional forms through the vertical, element 0) sliding, 1) deformation
        coef  = [lambda s:1.0, lambda s:1./4.*(5*s**4-1.)] 
        dcoef = [lambda s:0.0, lambda s:5*s**3]           # Derivatives of above forms

        # This is a quadrature rule for vertical integration
        points  = array([0.0,0.4688,0.8302,1.0])
        weights = array([0.4876/2.,0.4317,0.2768,0.0476])

        # This object does the integration
        self.vi = VerticalIntegrator(points,weights)

        self.U    = U
        self.Hmid = Hmid
        self.H0   = H0
        self.H    = H
        self.width= width
        self.B    = B
        self.S    = S
        self.u    = VerticalBasis([U[0]  ,  U[1]],coef,dcoef)
        self.phi  = VerticalBasis([Phi[0],Phi[1]],coef,dcoef)

        # The stresses ready to be accessed
        self.tau_d       = self.vi.intz(self._tau_d)
        self.tau_xx      = self.vi.intz(self._tau_xx)
        self.tau_xz      = self.vi.intz(self._tau_xz)
        self.tau_xy      = Phi[0]      * self._tau_xy()
        self.tau_b       = self.phi(1) * self._tau_b()
        # Forms of stresses ammenable to plotting (no test functions)
        self.tau_b_plot  = self._tau_b()
        self.tau_d_plot  = rho*g*self.Hmid*self.S.dx(0)
        self.tau_xx_plot = self.vi.intz(lambda s: self._dsdx(s)*self.Hmid*self._eta_v(s)*\
                                  (4*(self.u.dx(s,0) + self.u.ds(s)*self._dsdx(s)))\
                                   +self._dsdx(s)*self.H0*self._eta_v(s)*\
                                   (2*self.u(s)/self.width*self.width.dx(0)))
        self.tau_xz_plot = self.vi.intz(lambda s:self._dsdz(s)**2*self.Hmid*self._eta_v(s)*self.u.ds(s))
        self.tau_xy_plot = self._tau_xy()

    # These last two do not require vertical integration
    # Note the impact this has on the test functions associated with them.
    def _tau_xy(self):
        return 2.*self.Hmid*b/self.width*((n+2)/(self.width))**(1./n)*(abs(self.U[0])+1.0)**(1./n - 1)*self.U[0]

    def _tau_b(self):
        P_i = rho*g*self.Hmid                 # Overburden, or hydrostatic pressure of ice
        P_w = Max(-rho_w*g*self.B,thklim)     # Water pressure is either ocean pressure or zero
        N   = Max((P_i - P_w)/(rho*g),thklim) # Effective pressure is P_i - P_w, rho*g appears in A_s
        normalx = (self.B.dx(0))/sqrt((self.B.dx(0)**2 + 1.0))
        normalz = sqrt(1 - normalx**2)
        return mu*A_s*(abs(self.u(1))+1.0)**(1./n-1)*self.u(1)*(abs(N)+1.0)**(1./m)*(1-normalx**2)

    # Sigma coordinate derivative transforms
    def _dsdx(self,s):
        return 1./self.Hmid*(self.S.dx(0) - s*self.H0.dx(0))  # was H on second term!

    def _dsdz(self,s):
        return -1./self.Hmid

    def _eta_v(self,s):
        return b/2.*((self.u.dx(s,0) + self.u.ds(s)*self._dsdx(s))**2 \
                    +0.25*((self.u.ds(s)*self._dsdz(s))**2) \
                    + eps_reg)**((1.-n)/(2*n))

    def _tau_xx(self,s):
        return (self.phi.dx(s,0) + self.phi.ds(s)*self._dsdx(s))*self.Hmid*self._eta_v(s)\
        *(4*(self.u.dx(s,0) + self.u.ds(s)*self._dsdx(s))) + (self.phi.dx(s,0) \
        + self.phi.ds(s)*self._dsdx(s))*self.H0*self._eta_v(s)*\
        (2*self.u(s)/self.width*self.width.dx(0))       # Final Hmid was H0

    def _tau_xz(self,s):
        return self._dsdz(s)**2*self.phi.ds(s)*self.Hmid*self._eta_v(s)*self.u.ds(s)
    def _tau_d(self,s):
        return rho*g*self.Hmid*self.S.dx(0)*self.phi(s)
