import netCDF4 as ncdf
from scipy.interpolate import RegularGridInterpolator, RectBivariateSpline, griddata, SmoothBivariateSpline
from pylab import *
import numpy as np
import h5py

''''beep boop
    NOTES:
    RACMO DATA
    lon = x            range:   (-639456 , 855544)
    lat = y         range:  (-3355096 ,-656096)
    Spacing is 1,000 meters
    
    VELOCITY DATA
    x                range:   (-638000, 864550)
    y                range:   (-657600, -3349350)
    Spacing is 150 meters
    
'''
#
# # VELOCITY
# vel_x0 = -638000  # first x coordinate
# vel_x1 = 864550  # last x coordinate
# vSpacing = 150
# velx_num = 1 + (vel_x1 - vel_x0) / vSpacing
#
# vel_y0 =  -657600  # first y coordinate
# vel_y1 = -3349350  # last y coordinate
# vely_num = 1 + (math.fabs(vel_y1 - vel_y0)) / vSpacing
#

# RACMO

r_x0 = -639456
r_y0 = -3355096

r_x1 = 855544
r_y1 = -656096

rSpacing = 1000

rx_num = 1 +(r_x1 - r_x0)/rSpacing
ry_num = 1 + (math.fabs(r_y1 - r_y0))/rSpacing

# Bed machine

bed_x0 = -637925  # first x
bed_x1 = 864625  # last x

bed_y0 = -657675
bed_y1 = -3349425

bed_xarray = linspace(bed_x0, bed_x1, 10018, endpoint=True) # FIXME should maybe be one less point?? Prob not because +150 on one side, -150 on the other
bed_yarray = linspace(bed_y1, bed_y0,  17946, endpoint=True)


r_xarray = linspace(r_x0, r_x1, rx_num, endpoint=True)
r_yarray = linspace(r_y0, r_y1, int(ry_num), endpoint=True)

'''
    PROBLEM
    Arrays are (58,x,y) can't get all into one array!!!!!!!

'''

# subSample = 10

def maskConv():
    maskData = ncdf.Dataset('/home/pat/research/Data/Icemask_Topo_Iceclasses_lon_lat_average_1km.nc')
    mask = maskData.variables['Icemask'][:]
    maskIGrid = RectBivariateSpline(r_xarray, r_yarray, mask.transpose())
    maskData.close()


def smbConv(subSample):
    '''
    ONLY DOES THE LATEST YEAR!
    :return:
    '''
    outF = h5py.File('/home/pat/research/Greenland/data/GreenlandInBedCoord.h5','a')
    if 'smb' in outF.keys():
        outF.__delitem__('smb')
    smbData = ncdf.Dataset('/home/pat/research0/Data/SMB_rec_corr_v1.0.1958-2015.BN_1958_2013_1km.YY.nc')
    smb = smbData.variables['SMB_rec'][:][57][::subSample, ::subSample]
    smbData.close()

    smbIGrid = RectBivariateSpline(r_xarray[::subSample], r_yarray[::subSample], smb.transpose())
    smbIVals = np.array(smbIGrid(bed_xarray[::subSample], bed_yarray[::subSample]))  #memoryError

    smbIVals = smbIVals.transpose()
    smbIVals = np.flipud(smbIVals)
    print 'shape: ', smbIVals.shape

    outF.create_dataset("smb", data=smbIVals)
    outF.close()




def t2mConv():
    t2mData = ncdf.Dataset('/home/pat/research/Data/t2m.1958-2015.BN_1958_2013_1km.YY.nc')
    t2m = t2mData.variables['t2m'][:][57]
    # print t2m.shape
    t2mIGrid = RectBivariateSpline(r_xarray, r_yarray, t2m.transpose())
    t2mIVals = np.array(t2mIGrid(bed_xarray, bed_yarray))
    # t2mIVals = np.flipud(t2mIVals)
    t2mIVals = t2mIVals.transpose()
    outF = h5py.File('/home/pat/research/Greenland/data/GreenlandInBedCoord.h5','a')
    outF.__delitem__('t2m')
    outF.create_dataset("t2m", data=t2mIVals)
    outF.close()
    t2mData.close()

# t2mConv()

def mergeSMB():
    fl = h5py.File('/home/pat/research/Data/AllDataSets.h5', 'r')
    fl.create_dataset('SMB_rec', data=fl['SMB_rec_0'][:])
    for i in range(1,58):
        print 'loop ' + str(i)
        smb = [fl['SMB_rec'][:]]
        fl.__delitem__('SMB_rec')
        nm = 'SMB_rec_' + str(i)
        smb.append(fl[nm][:])
        fl.create_dataset('SMB_rec', data=np.array(smb))
        del smb
    fl.close()
    print 'done!'


