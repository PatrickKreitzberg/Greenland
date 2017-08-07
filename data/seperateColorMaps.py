import h5py


def createAll():
    datF = h5py.File('dataCMValues.h5', 'r')

    velF = h5py.File('velocityCM.h5','w')
    velF.create_dataset('colormap',data=datF['velocity'][:])
    velF.close()

    bedF = h5py.File('bedCM.h5','w')
    bedF.create_dataset('colormap',data=datF['bed'][:])
    bedF.close()

    smbF = h5py.File('smbCM.h5','w')
    smbF.create_dataset('colormap',data=datF['smb'][:])
    smbF.close()

    surfaceF = h5py.File('surfaceCM.h5','w')
    surfaceF.create_dataset('colormap',data=datF['surface'][:])
    surfaceF.close()

    thicknessF = h5py.File('thicknessCM.h5','w')
    thicknessF.create_dataset('colormap',data=datF['thickness'][:])
    thicknessF.close()

    datF.close()

# createAll()


