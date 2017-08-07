import numpy as np
import matplotlib.pyplot as plt

# Local imports.
import distmesh as dm

# Polygon example:
# pv are the vertices
# @dpoly gives the distance function
# @huniform Implements the trivial uniform mesh size function h=1.

def polygon():  # GOES COUNTER_CLOCKWISE
    """Polygon"""
    pv = np.array([(-0.4,-0.5),(0.4,-0.2),(0.4,-0.7),(1.5,-0.4),(0.9,0.1),
                   (1.6,0.8),(0.5,0.5),(0.2,1.0),(0.1,0.4),(-0.7,0.7),
                   (-0.4,-0.5)])
    fd = lambda p: dm.dpoly(p, pv)
    return dm.distmesh2d(fd, dm.huniform, 0.1, (-1,-1, 2,1), pv)



pause = lambda : None

plt.ion()
np.random.seed(1) # Always the same results

def fstats(p, t):
    print('%d nodes, %d elements, min quality %.2f'
          % (len(p), len(t), dm.simpqual(p,t).min()))

p, t = polygon()
fstats(p, t)
pause(); print('')

print p
print t