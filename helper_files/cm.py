import pyqtgraph as pg
from colorbar import *

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

    elif dataSet == 'bed':
        bedMin, bedMax = -5052.39510208, 3675.78314839
        lcArr = []
        linCMFile = open('./helper_files/bathCPT.txt', 'r')
        for line in linCMFile:
            vals = line.split()
            v0 = []
            v0.append(int(vals[1]))
            v0.append(int(vals[2]))
            v0.append(int(vals[3]))
            lcArr.append(v0)

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
        linPos1 = [0, 500, 1000, 1500, 2000, 2700]
        bedColor = lcArr + lcArr1
        bedPos = linPos + linPos1
        return pg.ColorMap(bedPos, bedColor, mode='byte')

    elif dataSet == 'surface':
        colors = [
            [255, 0, 255],
            [0, 132, 53],
            [51, 204, 0],
            [244, 240, 113],
            [244, 189, 69],
            [153, 100, 43],
            [255, 255, 255],
        ]
        mx = 3677
        spos = [
            -324,
            0,
            .125 * mx,
            .25 * mx,
            .5 * mx,
            .75 * mx,
            mx,
        ]
        return pg.ColorMap(spos, colors, mode='byte')


def getColorBar(dataSet, cm):
    if dataSet == 'velocity':
        cbItem = LogColorBar(cm, 10, 200, label='Velocity (m/yr)', tick_labels=['0', '10', '100', '1000', '8000'],
                             ticks=[0, 10, 100, 1000, 8000])
        cbItem.scale(50, 50)
        cbItem.translate(400.0, 90.0)  # may need to be in the regular script
        return cbItem

    elif dataSet == 'bed':
        cbItem = ColorBar(cm, 10, 200, label='Bed Ele. (m)',
                          tick_labels=['-1000', '-500', '0', '500', '1000', '1500', '2000', '2700'],
                          ticks=[-1000, -500, 0, 500, 1000, 1500, 2000, 2700])  # , [0., 0.5, 1.0])
        cbItem.scale(50, 50)
        cbItem.translate(400.0, 90.0)
        return cbItem
    elif dataSet == 'surface':
        cbItem = ColorBar(cm, 10, 200, label='Surface Ele. (m)',
                          tick_labels=['-323', '0', '460', '920', '1839', '2758', '3677'],
                          ticks=[-323, 0, 460, 920, 1839, 2758, 3677])  # , [0., 0.5, 1.0])
        cbItem.scale(50, 50)
        cbItem.translate(400.0, 90.0)
        return cbItem


########################  PROCESS DATA  #################################
# def getVelocityCM():
#     minimum, maximum =  0.0, 12950.7
#     c = [[[255,255,255],0.01],    # zero is white
#         [[0,76,153], 1],        # dark blue 1 - >10
#         [[0,0,200],10],
#         [[0,255,255],100],
#         [[0,153,76],1000],
#         [[255,255,0], maximum],
#         [[200,100,0],0],
#         [[255, 0, 0], 0],
#         [[153, 0, 0], 0],
#         [[120, 0, 0], 0]]
#
#     logmax = np.log10(maximum - 0.01)
#     div = logmax/(len(c)-2)
#     linearEnd = 0.01
#     # for i in range(1,int(linearEnd)):
#     #     pos.append(i/linearEnd)
#     pos = []
#     pos.append(0)
#     pos.append(0.01)
#     for i in range(1,len(c)-2):#int(linearEnd)):
#         pos.append(np.power(10,i*div))
#     pos.append(8000)
#     return pg.ColorMap(pos,[row[0] for row in c],mode='byte')

# def getVelocityColorbar():
#     cm = getVelocityCM()
#     cbItem = LogColorBar(cm, 10, 200, label='Velocity (m/yr)', tick_labels=['0', '10', '100', '1000', '8000'],
#                          ticks=[0, 10, 100, 1000, 8000])
#     cbItem.scale(50, 50)
#     cbItem.translate(400.0, 90.0) # may need to be in the regular script
#     return cbItem

# def getBedCM():
#     bedMin, bedMax = -5052.39510208,    3675.78314839
#     lcArr = []
#     linCMFile = open('./bathCPT.txt', 'r')
#     for line in linCMFile:
#         vals = line.split()
#         v0 = []
#         v0.append(int(vals[1]))
#         v0.append(int(vals[2]))
#         v0.append(int(vals[3]))
#         lcArr.append(v0)
#
#     div = float(np.abs(bedMin)) / float(len(lcArr) - 1)
#     linPos = []
#     for i in range(len(lcArr)):
#         linPos.append(bedMin + i * div)
#
#     ######### ELEVATION COLORMAP
#     lcArr1 = [
#      [0,153,0],
#      [102,255,102],
#      [255,255,0],
#      [255,0,0],
#      [153,0,0],
#      [255,255,255]
#     ]
#     linPos1 = [0,500,1000,1500,2000, 2700]
#     bedColor = lcArr + lcArr1
#     bedPos = linPos + linPos1
#     return pg.ColorMap(bedPos,bedColor,mode='byte')

# def getBedColorbar():
#     cm = getBedCM()
#     cbItem = ColorBar(cm, 10, 200, label='Bed Ele. (m)', tick_labels=['-1000', '-500', '0', '500', '1000', '1500', '2000', '2700'],
#                          ticks=[-1000, -500, 0, 500, 1000, 1500, 2000, 2700])  # , [0., 0.5, 1.0])
#     cbItem.scale(50, 50)
#     cbItem.translate(400.0, 90.0)
#     return cbItem


# def getSurfaceCM():
#     colors =[
#         [255,0,255],
#         [0,132,53],
#         [51,204,  0],
#         [244,240,113],
#         [244,189, 69],
#         [153,100, 43],
#         [255,255,255],
#     ]
#     mx = 3677
#     pos = [
#         -324,
#         0,
#         .125*mx,
#         .25*mx,
#         .5*mx,
#         .75*mx,
#         mx,
#     ]
#     return pg.ColorMap(pos, colors, mode='byte')

# def getSurfaceColorbar():
#     cm = getSurfaceCM()
#     cbItem = ColorBar(cm, 10, 200, label='Surface Ele. (m)',
#                       tick_labels=['-323', '0', '460', '920', '1839', '2758', '3677'],
#                       ticks=[-323, 0, 460, 920, 1839, 2758, 3677])  # , [0., 0.5, 1.0])
#     cbItem.scale(50, 50)
#     cbItem.translate(400.0, 90.0)
#     return cbItem


