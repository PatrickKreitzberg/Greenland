import pickle
import h5py
import numpy as np
from pylab import sqrt, linspace
from scipy.interpolate import RectBivariateSpline

bed_x0 = -637925  # first x
bed_x1 = 864625  # last x

bed_y0 = -657675
bed_y1 = -3349425

bed_xarray = linspace(bed_x0, bed_x1, 10018, endpoint=True) # FIXME should maybe be one less point?? Prob not because +150 on one side, -150 on the other
bed_yarray = linspace(bed_y1, bed_y0, 17946, endpoint=True)

data = h5py.File('GreenlandInBedCoord.h5','r')
vx = data['VX'][:]
vy = data['VY'][:]
data.close()

v = sqrt(vx**2 + vy**2)
velInterp = RectBivariateSpline(bed_xarray, bed_yarray, np.flipud(v.transpose()))

velOut = open('velInterp.pkl', 'wb')
pickle.dump(velInterp, velOut)
velOut.close()