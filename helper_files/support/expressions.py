from dolfin import Expression,exp

# For idealized cases these are functions, will later be read from files.
# Surface mass balance, if using a mathematical expression: 

class Adot(Expression):
  def eval(self,values,x):
    '''

    :param values: value 1d array of length 1 with value ~0.69 in one case
    :param x:
    :return:
    '''
    # L = mesh_max-mesh_min
    # values[0] = 2.0 - (x[0]-mesh_min) / L * 6.0

    values[0] = 8.7/(1 + 200*exp(0.05e-3*(x[0]-300.e3))) - 8.0

# Width of glacier, if using mathematical expression
class Width(Expression):
  def eval(self,values,x):
    values[0] = (75000./(1 + 200*exp(0.05e-3*(x[0]-300e3))) + 5000)

