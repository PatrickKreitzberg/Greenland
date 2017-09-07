import netCDF4 as ncdf
from colorMaps import *
from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline
import h5py
from constants import *
from math_functions import *

def createCMData():
    subSample = 1

    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    vx = allData['VX'][:][::subSample, ::subSample]
    vy = allData['VY'][:][::subSample, ::subSample]
    v = sqrt(vx**2 + vy**2)
    bed = allData['bed'][:][::subSample, ::subSample]
    smb = allData['smb'][:][::subSample, ::subSample]
    thickness = allData['thickness'][:][::subSample, ::subSample]
    surface = allData['surface'][:][::subSample, ::subSample]
    allData.close()


    velVals = getCM('velocity').mapToByte(v)
    velVals = velVals.astype(np.uint8)

    bedVals = getCM('bed').mapToByte(bed)
    bedVals = bedVals.astype(np.uint8)

    surfaceVals = getCM('surface').mapToByte(surface)
    surfaceVals = surfaceVals.astype(np.uint8)

    smbVals = getCM('smb').mapToByte(smb)
    smbVals = smbVals.astype(np.uint8)

    thicknessVals = getCM('thickness').mapToByte(thickness)
    thicknessVals = thicknessVals.astype(np.uint8)


    outF = h5py.File('../data/dataCMValues.h5','w') #719
    outF.create_dataset('velocity', data=velVals)
    outF.create_dataset('bed', data=bedVals)
    outF.create_dataset('surface', data=surfaceVals)
    outF.create_dataset('smb', data=smbVals)
    outF.create_dataset('thickness', data=thicknessVals)
    outF.close()




def createSMB():
    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    smb = allData['smb'][:]
    print 'smb shape  ', allData['smb'][:].shape
    print 'thick shape', allData['thickness'][:].shape
    allData.close()

    smbVals = getCM('smb').mapToByte(smb)
    smbVals = smbVals.astype(np.uint8)

    outF = h5py.File('../data/dataCMValues.h5','a') #719

    outF.__delitem__('smb')

    outF.create_dataset('smb', data=smbVals)
    outF.close()

def createVel():
    subSample = 1
    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    vx = allData['VX'][:]#[::subSample, ::subSample]
    vy = allData['VY'][:]#::subSample, ::subSample]
    v = sqrt(vx ** 2 + vy ** 2)
    allData.close()
    velVals = getCM('velocity').mapToByte(v)
    velVals = velVals.astype(np.uint8)

    outF = h5py.File('../data/dataCMValues.h5', 'a')  # 719
    print outF.keys()
    outF.__delitem__('velocity')
    outF.create_dataset('velocity', data=velVals)
    outF.close()



def testVel():
    # BedMachineGreenland-2017-05-10.nc
    # bm = ncdf.Dataset('/home/pat/research0/Data/BedMachineGreenland-2017-05-10.nc', 'r')
    # print 'x:', bm.variables['x'][:][0], bm.variables['x'][:][-1]
    # print 'y:', bm.variables['y'][:][0], bm.variables['y'][:][-1]
    # bm.close()

    subSample = 1

    bed_x0 = -637925  # first x
    bed_x1 = 864625  # last x

    bed_y0 = -657675
    bed_y1 = -3349425

    bed_x00 = 0  # first x
    bed_x10 = 10018  # last x

    bed_y00 = 17946
    bed_y10 = 0

    bed_xarray = linspace(bed_x0, bed_x1, 10018, endpoint=True)
    bed_yarray = linspace(bed_y1, bed_y0, 17946, endpoint=True)

    bed_xarray0 = linspace(bed_x00, bed_x10, 10018, endpoint=True)
    bed_yarray0 = linspace(bed_y10, bed_y00, 17946, endpoint=True)

    subSample=1

    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    vx = allData['VX'][:]#[::subSample, ::subSample]
    vy = allData['VY'][:]#[::subSample, ::subSample]
    v = sqrt(vx**2 + vy**2)
    print 'x shape: ', allData.keys()
    allData.close()

    # interp = RectBivariateSpline(bed_xarray[::subSample], bed_yarray[::subSample], v[::subSample, ::subSample].transpose())
    interp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(v).transpose())
    interp0= RectBivariateSpline(bed_xarray0, bed_yarray0, np.flipud(v.transpose()))
    vi = interp(bed_xarray, bed_yarray, grid=True)

    x8, y8 = 5000, 10000
    x9, y9 = 2567, 14944
    px, py = projCoord(x9, y9)
    bleh, py2 = projCoord(0, 17946 - y9)
    lx, ly = projCoord(map_x1, map_y1)
    print v[y9][x9], interp0([x9], [y9], grid=False)
    print v[y9][x9], interp( [px], [py], grid=False)
    print v[y9][x9], interp( [px], [py2], grid=False)
    print 'shape: ', vi.shape

    # velVals = getCM('velocity').mapToByte(vi)
    # velVals = velVals.astype(np.uint8)

    # outF = h5py.File('../data/dataCMValues.h5','a') #719
    # outF.__delitem__('velocity')
    #
    # outF.create_dataset('velocity', data=velVals)
    # outF.close()

# testVel()


def createSurface():
    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    surface = allData['surface'][:]

    surfaceVals = getCM('surface').mapToByte(surface)
    surfaceVals = surfaceVals.astype(np.uint8)
    allData.close()

    outF = h5py.File('../data/dataCMValues.h5', 'a')  # 719
    outF.__delitem__('surface')

    outF.create_dataset('surface', data=surfaceVals)
    outF.close()

def createThickness():
    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    d = allData['thickness'][:]
    allData.close()

    thickVals = getCM('thickness').mapToByte(d)
    thickVals = thickVals.astype(np.uint8)

    outF = h5py.File('../data/dataCMValues.h5', 'a')  # 719
    outF.__delitem__('thickness')

    outF.create_dataset('thickness', data=thickVals)
    outF.close()

def createOldThickness():
    bedData = ncdf.Dataset('/home/pat/research/Data/BedMachineGreenland-2017-05-10.nc')
    bed = bedData.variables['thickness'][:]
    bedData.close()
    thickVals = getCM('oldthick').mapToByte(bed)
    thickVals = thickVals.astype(np.uint8)
    #

    outF = h5py.File('../data/dataCMValues.h5', 'a')  # 719
    outF.__delitem__('oldthick')
    outF.create_dataset('oldthick', data=thickVals)
    outF.close()

