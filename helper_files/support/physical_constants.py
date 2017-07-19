from dolfin import Constant
spy = 60**2*24*365          # Seconds per year
rho = 917.                  # Ice density
rho_w = 1029.0              # Seawater density
g = 9.81                    # Gravitational acceleration
n = 3.0                     # Glen's exponent
m = 3.0                     # Sliding law exponent
b = 1e-17**(-1./n)          # Ice hardness
mu = Constant(1.0)          # Basal sliding law constants:
A_s = Constant(rho*g*315.0/500.)
