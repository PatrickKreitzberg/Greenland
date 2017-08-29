

import pyqtgraph.examples
pyqtgraph.examples.run()

import threading
#
#   THREADING EXAMPLE
#


import threading
import time

exitFlag = 0

class myThread (threading.Thread):
   def __init__(self, threadID, name, counter):
      threading.Thread.__init__(self)
      self.threadID = threadID
      self.name = name
      self.counter = counter
   def run(self):
      print "Starting " + self.name
      print_time(self.name, self.counter, 5)
      print "Exiting " + self.name

def print_time(threadName, counter, delay):
   while counter:
      if exitFlag:
         threadName.exit()
      time.sleep(delay)
      print "%s: %s" % (threadName, time.ctime(time.time()))
      counter -= 1

# Create new threads
thread1 = myThread(1, "Thread-1", 1)
thread2 = myThread(2, "Thread-2", 2)

# Start new Threads
thread1.start()
thread2.start()

print "Exiting Main Thread"


#
#   THREADING EXAMPLE TWO
#

import time

def myfunc(i):
    print "sleeping 5 sec from thread %d" % i
    time.sleep(5)
    print "finished sleeping from thread %d" % i

for i in range(10):
    t = threading.Thread(target=myfunc, args=(i,))
    t.start()


'''

#
#   MULTIPROCESSING EXAMPLE
#
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

'''

