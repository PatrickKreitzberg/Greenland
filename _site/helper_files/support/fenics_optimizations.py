from dolfin import parameters, NonlinearVariationalSolver

parameters['allow_extrapolation'] = True


def set_solver_options(nls):
    nls.parameters['nonlinear_solver'] = 'snes'
    nls.parameters['snes_solver']['method'] = 'vinewtonrsls'
    nls.parameters['snes_solver']['relative_tolerance'] = 1e-6
    nls.parameters['snes_solver']['absolute_tolerance'] = 1e-6
    nls.parameters['snes_solver']["preconditioner"]="hypre-amg"
    nls.parameters['snes_solver']['error_on_nonconvergence'] = True
    nls.parameters['snes_solver']['linear_solver'] = 'mumps'
    nls.parameters['snes_solver']['maximum_iterations'] = 100
    nls.parameters['snes_solver']['report'] = True
