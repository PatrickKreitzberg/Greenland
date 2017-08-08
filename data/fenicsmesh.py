from dolfin.cpp.mesh import *
from dolfin.cpp.io import *

import h5py

mesh = Mesh()

me = RectangleMesh(Point(0,0), Point(10,10),10,10)#mesh.mpi_comm(),0.0,0.0,10.0,10.0,5,5)
# print me.cells()
# print me.domains()
# print me.topology()
# print me.ufl_cell()
# print me.coordinates()
mesh = Mesh('example.xml')
print mesh.coordinates()
print mesh.coordinates()[::, 0]
# print mesh.cells()
f = h5py.File('mesh2d.h5','r')
print f.keys()
f.close()


f = HDF5File(mesh.mpi_comm(),'mesh2d.h5','r')
print f.has_dataset('bed')

f.close()