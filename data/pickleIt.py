import pickle
import h5py
import numpy as np
from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline

import time
bed_x0 = -637925  # first x
bed_x1 = 864625  # last x

bed_y0 = -657675
bed_y1 = -3349425

bed_xarray = linspace(bed_x0, bed_x1, 10018, endpoint=True) # FIXME should maybe be one less point?? Prob not because +150 on one side, -150 on the other
bed_yarray = linspace(bed_y1, bed_y0, 17946, endpoint=True)


def velPickle():
    data = h5py.File('GreenlandInBedCoord.h5','r')
    vx = data['VX'][:]
    vy = data['VY'][:]
    data.close()

    v = sqrt(vx**2 + vy**2)
    velInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(v.transpose()))

    velOut = open('velInterp.pkl', 'wb')
    pickle.dump(velInterp, velOut)
    velOut.close()


def pickleItMeth(dname, fname):
    data = h5py.File('GreenlandInBedCoord.h5', 'r')
    d = data[dname][:]
    data.close()
    i = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(d.transpose()))
    out = open(fname, 'wb')
    pickle.dump(i, out)
    out.close()

# data = h5py.File('GreenlandInBedCoord.h5', 'r')
# print data.keys() #[u'VX', u'VY', u'bed', u'smb', u'surface', u't2m', u'thickness']
# data.close()

# pickleItMeth('bed', 'bedInterp.pkl')
# pickleItMeth('smb', 'smbInterp.pkl')
# pickleItMeth('surface', 'surfaceInterp.pkl')
# pickleItMeth('thickness', 'thicknessInterp.pkl')


def testLoad():
    st = time.time()
    print 'starting: '
    pickFile = open('./bedInterp.pkl', 'r')
    interp = pickle.load(pickFile)
    pickFile.close()
    print 'finished: ', time.time() - st

def testInterpTime():
    st0 = time.time()
    print 'starting: '
    data = h5py.File('GreenlandInBedCoord.h5', 'r')
    d = data['VX'][:]
    data.close()
    i = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(d.transpose()))
    print 'finished: ', time.time() - st0


def testSS(i):
    st0 = time.time()
    print 'starting: '
    data = h5py.File('GreenlandInBedCoord.h5', 'r')
    d = data['VX'][:]
    data.close()
    i = RectBivariateSpline(bed_xarray[::i], bed_yarray[::i], np.flipud(d[::i, ::i].transpose()))
    print 'finished: ', time.time() - st0

# testLoad()
# testInterpTime() # 11.53 seconds
# testSS(4)

a = [[1,2,3,4],[5,6,7,8],[9,10,11,12],[13,14,15,16]]
na = np.array(a)
print na[::1, ::1]


