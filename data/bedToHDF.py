import netCDF4 as ncdf
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline, griddata, SmoothBivariateSpline
from pylab import *
import numpy as np
import h5py


bedData = ncdf.Dataset('/home/pat/research/Data/BedMachineGreenland-2017-05-10.nc')
bed = bedData.variables['bed'][:]
surface = bedData.variables['surface'][:]
thickness = bedData.variables['thickness'][:]
bedData.close()

outF = h5py.File('/home/pat/research/Greenland/data/GreenlandInBedCoord.h5', 'a')
outF.create_dataset('bed', data=bed[:])
outF.create_dataset('surface', data=surface[:])
outF.create_dataset('thickness', data=thickness[:])
outF.close()