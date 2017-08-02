import numpy as np
import pyqtgraph as pg
from helper_files.classes.Colorbar import *
from b import geta
from constants import *

def getCM(dataSet):
    if dataSet == 'velocity':
        minimum, maximum = 0.0, 12950.7
        c = [[[255, 255, 255], 0.01],  # zero is white
             [[0, 76, 153], 1],  # dark blue 1 - >10
             [[0, 0, 200], 10],
             [[0, 255, 255], 100],
             [[0, 153, 76], 1000],
             [[255, 255, 0], maximum],
             [[200, 100, 0], 0],
             [[255, 0, 0], 0],
             [[153, 0, 0], 0],
             [[120, 0, 0], 0]]

        logmax = np.log10(maximum - 0.01)
        div = logmax / (len(c) - 2)
        linearEnd = 0.01
        # for i in range(1,int(linearEnd)):
        #     pos.append(i/linearEnd)
        pos = []
        pos.append(0)
        pos.append(0.01)
        for i in range(1, len(c) - 2):  # int(linearEnd)):
            pos.append(np.power(10, i * div))
        pos.append(8000)
        return pg.ColorMap(pos, [row[0] for row in c], mode='byte')

    elif dataSet == 'thickness' or dataSet == 'oldthick':
        thickMin = 0.0
        thickMax = 3428
        colors = [[255,0,255],
                  [255,255,255],
                  [6,6,94],
                  [0, 213, 253],
                  [255,253,3],
                  [138,0,0],
                  [110, 0, 0]]
        pos = [-0.00001, 0.0, 0.00001, 1000.0, 2000.0, 3000.0, 3500.0]
        return pg.ColorMap(pos, colors, mode='byte')

    elif dataSet == 'bed':
        bedMin, bedMax = -5052.39510208, 3675.78314839
        lcArr = geta()
        div = float(np.abs(bedMin)) / float(len(lcArr) - 1)
        linPos = []
        for i in range(len(lcArr)):
            linPos.append(bedMin + i * div)

        ######### ELEVATION COLORMAP
        lcArr1 = [
            [0, 153, 0],
            [102, 255, 102],
            [255, 255, 0],
            [255, 0, 0],
            [153, 0, 0],
            [255, 255, 255]
        ]
        linPos1 = [0, 3678/5, (2*3678)/5, (3*3678)/5,(4*3678)/5,(5*3678)/5]
        bedColor = lcArr + lcArr1
        bedPos = linPos + linPos1
        return pg.ColorMap(bedPos, bedColor, mode='byte')

    elif dataSet == 'surface':
        colors = [
            [255,   0, 255],
            [255, 255, 255],
            [0,   132,  53],
            [51,  204,   0],
            [244, 240, 113],
            [244, 189,  69],
            [178,   0,   0],
            [255, 255, 255],
        ]
        mx = 3677
        spos = [
            -324,
            0,
            0.0001,
            .125 * mx,
            .25 * mx,
            .5 * mx,
            .75 * mx,
            mx,
        ]
        return pg.ColorMap(spos, colors, mode='byte')

    elif dataSet == 'smb':
        colorsSMB = [[255,   0,   0],
                  [255, 255, 255],
                  [0,     0, 255]]
        posSMB = [-11000, 0, 6061]
        return pg.ColorMap(posSMB, colorsSMB)


def getColorBar(dataSet, cm):
    colorbarHeight = 200
    colorbarWidth = 20
    if dataSet == 'velocity':
        cbItem = LogColorBar(cm, colorbarWidth, colorbarHeight,
                             label='Velocity(m/yr)',
                             tick_labels=['0', '10', '100', '1,000', '8,000'],
                             ticks=[0, 10, 100, 1000, 8000])
        return cbItem

    elif dataSet == 'bed':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='Bed Ele.(m)',
                          tick_labels=['-5,000','-4,000','-3,000','-2,000','-1,000','0', '1,000', '2,000', '3,000', '3,500'],
                          ticks=[-5000,-4000,-3000,-2000,-1000,0, 1000, 2000, 3000, 3500])
        return cbItem
    elif dataSet == 'surface':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='Surface Ele.(m)',
                          tick_labels=['0', '500', '1,000', '1,500', '2,000', '2,500', '3,000', '3,500'],
                          ticks=[0, 500, 1000, 1500, 2000, 2500, 3000, 3500])
        return cbItem

    elif dataSet == 'smb':
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='SMB(m)',
                          tick_labels=['-11', '0',  '6'],
                          ticks=[-11000, 0, 6000],
                          name='smb')
        return cbItem

    elif dataSet == 'thickness' or dataSet == 'oldthick':
        #thickness range:  0 3428
        cbItem = ColorBar(cm, colorbarWidth, colorbarHeight,
                          label='Thickness(m)',
                          tick_labels=['0', '500', '1000', '1500', '2000', '2500', '3000', '3500'],
                          ticks=[0, 500, 1000, 1500, 2000, 2500, 3000, 3500])
        return cbItem