import netCDF4 as ncdf

from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline, griddata, SmoothBivariateSpline
from pylab import *
import numpy as np
import h5py

np.random.seed(7)
'''
clean version of sampledataset_projectvelontobedmachine_regulargridinter.py
since transform doesnt do anything.
'''
##############################
##      GATHER DATA         ##
##############################

bedData = ncdf.Dataset('/home/pat/research/Data/BedMachineGreenland-2017-05-10.nc')
bedX = bedData.variables['x'][:]
bedY = bedData.variables['y'][:]
bed = bedData.variables['thickness'][:]

velocityData = ncdf.Dataset('/home/pat/research/Data/Greenland_ice_speed_v26May2017.nc')
velX = velocityData.variables['x']
velY = velocityData.variables['y']
v = velocityData.variables['VX']
print 'vel shape: ', v.shape
print 'bed shape: ', bed.shape

##############################
##      GATHER BOUNDS       ##
##############################


# VELOCITY
vel_x0 = -638000  # first x coordinate
vel_x1 = 864550  # last x coordinate
spacing = 150
velx_num = 1 + (vel_x1 - vel_x0) / spacing

vel_y0 =  -657600  # first y coordinate
vel_y1 = -3349350  # last y coordinate
vely_num = 1 + (math.fabs(vel_y1 - vel_y0)) / spacing
print velx_num, vely_num

#bed machine
bed_x0 = -637925  # first x
bed_x1 = 864625  # last x

bed_y0 = -657675
bed_y1 = -3349425


bed_xarray = linspace(bed_x0, bed_x1, 10018, endpoint=True) # FIXME should maybe be one less point?? Prob not because +150 on one side, -150 on the other
bed_yarray = linspace(bed_y1, bed_y0,  17946, endpoint=True)

vel_xarray = linspace(vel_x0, vel_x1, 10018, endpoint=True)
vel_yarray = linspace(vel_y1, vel_y0,  17946, endpoint=True)
print vel_xarray[0], vel_xarray[1]
print vel_yarray[0], vel_yarray[1]

# bmx, bmy = meshgrid(bed_xarray, bed_yarray)
vmx, vmy = meshgrid(vel_xarray, vel_yarray)


##############################
##    CREATE INTERPOLATOR   ##
##############################

print 'bedT be4  min, max: ', np.amin(np.flipud(bed.transpose())), np.amax(np.flipud(bed.transpose()))
# # Create spline from source (Mouginot) data and source (transformed) data coordinates
bedIGrid = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(bed.transpose())) #FIXME is this upside-down??

##############################
##     INTERPOLATE DATA     ##
##############################


bedIVals = np.array(bedIGrid(vel_xarray, vel_yarray))
print 'bed be4  min, max: ', np.amin(bed), np.amax(bed)
print 'bedIVals min, max: ', np.amin(bedIVals), np.amax(bedIVals)
print 'shape before transpose(): ', bed.shape
bedIVals = np.flipud(bedIVals)
bedIVals = bedIVals.transpose()

print 'shape after transpose(): ', bedIVals.shape
# bedIvals = np.flipud(bedIVals)

xx = 2516  #int(14-10018//2)# -637925 + 150*(10018/2)
yy = 13673 #int(150+17846//3)# -657600 - 150*(17946/2)

#-637925 + 150*(10018/2)
#-657600 - 150*(17946/2)
print 'Shapes:'
print bed.shape
print bedIVals.shape
print bed[yy][xx], bedIVals[yy][xx]
# print 'bed[-657600 - 150*(17946/2)][-637925 + 150*(10018/2)]: ', bed[yy][xx]
# approx = bedIGrid(bedX[xx],bedY[yy], grid=False)
# print 'Approximation: ', bedIVals[yy][xx]
# bxi = bxi.reshape(bed.shape)
print "bxi shape: ", bedIVals.shape
#
# bedOut = {}
# bedOut['bed'] = bedIVals
# savemat('/home/pat/Research/Data/bed.mat', bedOut)

# outF = h5py.File('/home/pat/research/Greenland/data/AllDataSets.h5','a')
# outF.create_dataset("oldthick", data=bed)
# outF.close()
# Create a new file using defaut properties.

# outF = h5py.File('/home/pat/research/Greenland/data/AllDataSets.h5','a')
# print 'created file'
# #
# # Create a dataset under the Root group.
# #
# dataset = outF.create_dataset("thickness", data=bedIVals)
#
# print 'created dataset, now filling dataset'
# print "Dataset dataspace is", dataset.shape
# print "Dataset Numpy datatype is", dataset.dtype
# print "Dataset name is", dataset.name
# # print "Dataset is a member of the group", dataset.parent
# print "Dataset was created in the file", dataset.file
#
#
# #
# # Close the file before exiting
# #
# outF.close()

print 'negative: ', (bed < 0).sum()
print 'negative: ', (bedIVals < -100).sum()
for i in range(len(bedIVals)):
    for j in range(len(bedIVals[i])):
        if bedIVals[i][j] < 0 and bedIVals[i][j] > -1:
            bedIVals[i][j] = 0
            # print i,j,bedIVals[i][j], '\t\t', bed[i][j]

outF = h5py.File('/home/pat/research/Greenland/data/AllDataSets.h5','a')
outF.__delitem__('thickness')
outF.create_dataset("thickness", data=bedIVals)
outF.close()


velocityData.close()
bedData.close()