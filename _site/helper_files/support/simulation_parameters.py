from dolfin import Constant
dt_float = 5.0              # Time step
thklim = 10.0
eps_reg = 1e-5              # Regularization parameter
dt = Constant(dt_float)
theta = Constant(0.5)       # Crank-Nicholson parameter
