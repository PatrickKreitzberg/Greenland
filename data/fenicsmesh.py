from dolfin.cpp.mesh import *


mesh = Mesh()

me = RectangleMesh(Point(0,0), Point(10,10),10,10)#mesh.mpi_comm(),0.0,0.0,10.0,10.0,5,5)
