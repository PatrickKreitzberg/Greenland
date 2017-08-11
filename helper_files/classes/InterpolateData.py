from fenics import *

#local imports
# from ..dataset_objects import *

class interpolateData(Expression):
    def __init__(self,interpolator, **kwargs):
        #print 'has element: ', V.has_element()
        self.interpolator = interpolator
    def eval(self, values, x):
        values[0] = self.interpolator(x[0], x[1])