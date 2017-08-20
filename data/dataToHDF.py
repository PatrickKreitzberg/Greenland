# import netCDF4 as ncdf
# from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline, griddata, SmoothBivariateSpline
# from pylab import *
# import numpy as np
# from RacmoToBed import *
# from VelocityToBedConversion import *
# from helper_files.createCMDatasets import *
# import time
import h5py


f = h5py.File('outModel.h5')

print f['bed']['10.0']

f.close()


'''

Converts all data to bed coordinates and puts all in hdf file


subSample = 6

bedData = ncdf.Dataset('/home/pat/research0/Data/BedMachineGreenland-2017-05-10.nc')
bed = bedData.variables['bed'][:][::subSample, ::subSample]
surface = bedData.variables['surface'][:][::subSample, ::subSample]
thickness = bedData.variables['thickness'][:][::subSample, ::subSample]
# x = bedData.variables['x'][:][::subSample]
# y = bedData.variables['y'][:][::subSample]
bedData.close()

outF = h5py.File('/home/pat/research/Greenland/data/GreenlandInBedCoord.h5', 'a')
outF.create_dataset('bed', data=bed[:])
outF.create_dataset('surface', data=surface[:])
outF.create_dataset('thickness', data=thickness[:])
# outF.create_dataset('x', data=x[:])
# outF.create_dataset('y', data=y[:])
outF.close()

# convert velocity data
velocityData = ncdf.Dataset('/home/pat/research0/Data/Greenland_ice_speed_v26May2017.nc')
id(velocityData.variables['VX'][:][::subSample, ::subSample], 'VX', subSample)
id(velocityData.variables['VY'][:][::subSample, ::subSample], 'VY', subSample)
velocityData.close()

# convert surface mass balance data
smbConv(subSample)
# createCMData()

'''


