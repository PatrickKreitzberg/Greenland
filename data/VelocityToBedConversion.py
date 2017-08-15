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

#bed machine
bed_x0 = -637925  # first x
bed_x1 = 864625  # last x

bed_y0 = -657675
bed_y1 = -3349425


bed_xarray = linspace(bed_x0, bed_x1, 10018, endpoint=True) # FIXME should maybe be one less point?? Prob not because +150 on one side, -150 on the other
bed_yarray = linspace(bed_y1, bed_y0,  17946, endpoint=True)

vel_xarray = linspace(vel_x0, vel_x1, 10018, endpoint=True)
vel_yarray = linspace(vel_y1, vel_y0,  17946, endpoint=True)



# bmx, bmy = meshgrid(bed_xarray, bed_yarray)
# vmx, vmy = meshgrid(vel_xarray, vel_yarray)



# print 'bedT be4  min, max: ', np.amin(np.flipud(bed.transpose())), np.amax(np.flipud(bed.transpose()))
# # Create spline from source (Mouginot) data and source (transformed) data coordinates
# bedIGrid = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(bed.transpose())) #FIXME is this upside-down??



##############################
##     INTERPOLATE DATA     ##
##############################



def id(data, name, subSample):
    velocityGrid = RectBivariateSpline(vel_xarray[::subSample], vel_yarray[::subSample], np.flipud(data.transpose()))
    velIVals = np.array(velocityGrid(bed_xarray[::subSample], bed_yarray[::subSample]))
    velIVals = np.flipud(velIVals)
    velIVals = velIVals.transpose()
    outF = h5py.File('/home/pat/research/Greenland/data/GreenlandInBedCoord.h5','a')
    outF.create_dataset(name, data=velIVals)
    outF.close()

# velocityData = ncdf.Dataset('/home/pat/research0/Data/Greenland_ice_speed_v26May2017.nc')
# id(velocityData.variables['VX'][:][::subSample, ::subSample], 'VX')
# id(velocityData.variables['VY'][:][::subSample, ::subSample], 'VY')
# velocityData.close()
