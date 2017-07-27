import netCDF4 as ncdf
from cm import *
from pylab import sqrt
import h5py
def createCMData():
    allData = h5py.File('../data/GreenlandInBedCoord.h5', 'r')
    vx = allData['VX'][:]
    vy = allData['VY'][:]
    v = sqrt(vx**2 + vy**2)
    bed = allData['bed'][:]
    smb = allData['smb'][:]
    thickness = allData['thickness'][:]
    surface = allData['surface'][:]
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

