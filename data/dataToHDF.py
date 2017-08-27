import netCDF4 as ncdf
# from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline, griddata, SmoothBivariateSpline
# from pylab import *
# import numpy as np
# from RacmoToBed import *
# from VelocityToBedConversion import *
# from helper_files.createCMDatasets import *
# import time
import h5py

f = h5py.File('GreenlandInBedCoord.h5', 'a')
print f.keys()
print f['x'][:].shape
print f['smb'][:].shape
print f['surface'][:].shape
print f['VX'][:].shape
f.close()

def createLowResDataSet():
    res = 13 # * 150 meters
    out = h5py.File('lowRes.h5', 'w')
    inF = h5py.File('large_GreenlandInBedCoord.h5','r')
    out.create_dataset('VX',        data=inF['VX'][:][::res, ::res])
    out.create_dataset('VY',        data=inF['VY'][:][::res, ::res])
    out.create_dataset('bed',       data=inF['bed'][:][::res, ::res])
    out.create_dataset('smb',       data=inF['smb'][:][::res, ::res])
    out.create_dataset('surface',   data=inF['surface'][:][::res, ::res])
    out.create_dataset('thickness', data=inF['thickness'][:][::res, ::res])
    out.create_dataset('x',         data=inF['x'][:][::res])
    out.create_dataset('y',         data=inF['y'][:][::res])
    inF.close()
    out.close()

f = h5py.File('lowRes.h5', 'r')
print f.keys()
print f['x'][:].shape
print f['smb'][:].shape
print f['surface'][:].shape
print f['VX'][:].shape
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


