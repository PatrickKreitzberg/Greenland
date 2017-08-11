from dolfin.cpp.mesh import *
from dolfin.cpp.io import *

import h5py


f = h5py.File('dataCMValues.h5','r')
print f.keys()
f.close()
