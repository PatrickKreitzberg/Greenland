from dolfin import TestFunction, TrialFunction, solve, dx, Function
def smooth(Q,H,orig,kappa):
  # Gaussian smoothing:
  # INPUTS:
  # Q is a function space
  # H is the thickness
  # orig is the original field to be smoothed 
  # kappa ice thicknesses smoothed over
  # OUTPUT
  # u the smoothed version of orig
  # NOTES:
  # Due to scaling concerns with diffusivity (kappa*H)**2,
  # the orig field has to be in terms of a driving
  # stress, so it should be written as rho*g*H*orig

  phi        = TestFunction(Q)
  u          = TrialFunction(Q)

  a = + u * phi * dx \
         + (kappa*H)**2 * (phi.dx(0)*u.dx(0) + phi.dx(1)*u.dx(1)) * dx
  L = + orig * phi * dx 

  u = Function(Q)
  solve(a == L, u)
  return u
