from fenics import *

#local imports
# from ..dataset_objects import *

class interpolateDataClass(Expression):
    def __init__(self,interpolator, **kwargs):
        # Expression.__init__(self,kwargs)
        #print 'has element: ', V.has_element()
        self.interpolator = interpolator
    def eval(self, values, x):
        values[0] = self.interpolator(x[0], x[1])