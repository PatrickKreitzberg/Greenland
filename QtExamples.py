

# import pyqtgraph.examples
# pyqtgraph.examples.run()


import h5py
from multiprocessing import *
import time
dataFile = h5py.File('./data/GreenlandInBedCoord.h5', 'r')
def worker(i, num, out_q, dataFile):
    """thread worker function"""

    time.sleep(num)
    outdict = {}
    # print dataFile.keys()
    fd = dataFile['VX'][:]
    print fd
    # print 'Worker:', i, 'seconds', num
    outdict[i] = num
    out_q.put(outdict)

    return

out_q = Queue()

if __name__ == '__main__':
    jobs = []
    nms = [2,1,1,1,3]
    for i in range(5):
        p = Process(target=worker, args=(i,nms[i], out_q, dataFile))
        jobs.append(p)
        p.start()
    resultdict = {}
    for j in jobs:
        resultdict.update(out_q.get())

    for j in jobs:
        j.join()

    print resultdict
    print resultdict[0]

dataFile.close()

