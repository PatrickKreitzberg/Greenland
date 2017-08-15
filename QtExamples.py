import numpy as np


# import pyqtgraph.examples
# pyqtgraph.examples.run()

#
import h5py
# f = h5py.File('/home/pat/research/Greenland/data/GreenlandInBedCoord.h5','r')
# print f.keys()
# # f.__delitem__('VY')
# print f['smb'][:]
# print f.keys()
# f.close()

#
# import netCDF4 as ncdf
# bedData = ncdf.Dataset('/home/pat/research/Data/BedMachineGreenland-2017-05-10.nc')
# thickBed = bedData.variables['thickness'][:]
# print np.amin(thickBed[:])
# print np.amax(thickBed[:])
# # print bedData.keys()
# bedData.close()

from scipy import linspace
import numpy as np
path = linspace(0,150,15)

m = np.matrix([[1,1],[0,1]])


x = linspace(-1,1,2, endpoint=True)
y = [0]*2
print m*np.matrix([x, y])
z = m*np.matrix([x, y])
print ''
print z.shape
print z.itemsize
print z.item((z.shape[0]-1,0))
# print x+10