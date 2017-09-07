import numpy as np
from scipy.integrate import ode


# LOCAL IMPORTS
from math_functions import *
from pens import *
from gui import *
from dataset_objects import *
from data_functions import *
from classes.Marker import *
from velocity_functions import *
from classes.StaticPlotter import *

def centerVelocityStream(x, y):
    # x,y in color coordinates
    x0p, y0p = colorToProj(x, y)
    y0 = np.array([x0p, y0p])
    t0, t1, dt = 0, 80, .1
    r = ode(getProfile).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    ox = []
    oy = []
    ind = 0
    while r.successful() and ind < 2:
        '''
        what the fuck was I trying to do here??
        I think, integrate then look perpindicular for the velocity width margins
        '''
        ai = r.integrate(r.t + dt)
        ind += 1
        xi, yi = colorCoord(np.real(ai[0]), np.real(ai[1]))
        ox.append(xi)
        oy.append(yi)
    # iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=blackPlotPen))

    endPoints = [[0, 0], [0, 0]] # either side of the velocity stream
    endPoints[0][0], endPoints[0][1], endPoints[1][0], endPoints[1][1] = calcVelWidth(x, y, ox[-1], oy[-1], True)
    d = (0.5) * np.sqrt((endPoints[0][0] - endPoints[1][0]) ** 2 + (endPoints[0][1] - endPoints[1][1]) ** 2)
    theta = np.arctan2(float(y - oy[-1]), float(x - ox[-1]))
    x, y = endPoints[0][0] + (d * -np.sin(theta)), endPoints[0][1] + (d * np.cos(theta))
    return x, y

def mouseClick(e):
    global vptSel, vptCur#, integrateLine
    first = False
    if len(vpts) == 0:
        first = True
        cx = e.pos().x()
        cy = e.pos().y()
        if autoCorrectVpt.checkState() == 2:
            cx, cy = centerVelocityStream(cx, cy)
    if not vptSel:
        # If there is not a marker selected already
        for i in range(len(vpts)):
            if vpts[i].checkClicked(e.pos()):
                # See if you clicked on an already existing point
                if keysPress['shift']:
                    x0p, y0p = colorToProj(vpts[i].cx, vpts[i].cy)
                    y0 = np.array([x0p, y0p])
                    t0, t1, dt = 0, 80, .1
                    r = ode(getProfile).set_integrator('zvode', method='bdf')
                    r.set_initial_value(y0, t0)
                    ox = [vpts[i].cx]
                    oy = [vpts[i].cy]
                    while r.successful() and r.t < t1:
                        ai = r.integrate(r.t + dt)
                        xi, yi = colorCoord(ai[0], ai[1])
                        if np.sqrt((xi-ox[-1])**2 + (yi-oy[-1])**2):
                            ox.append(np.real(xi))
                            oy.append(np.real(yi))
                    print 'ox, oy', len(ox), len(oy)
                    print ox
                    print oy
                    intLines.append(pg.PlotDataItem(ox, oy, pen=whitePlotPen))
                    iiContainer.currentWidget().addItem(intLines[-1])
                elif keysPress['ctrl']:
                    # if you ctrl+click a maker it is deleted
                    if i > 0:
                        if i + 1 < len(vpts):
                            # connect line from previous node to next point
                            vpts[i - 1].lines[1] = pg.PlotDataItem([vpts[i - 1].cx, vpts[i + 1].cx], [vpts[i - 1].cy, vpts[i + 1].cy], connect='all', pen=skinnyBlackPlotPen)
                            vpts[i + 1].lines[0] = vpts[i - 1].lines[1]
                            iiContainer.currentWidget().addItem(vpts[i - 1].lines[1])
                            pg.QtGui.QApplication.processEvents()
                        else:  # delete line because there is no next point
                            iiContainer.currentWidget().removeItem(vpts[i - 1].lines[1])
                            vpts[i - 1].lines[1] = None
                    elif i == 0 and len(vpts) > 1:  # point is the first so delete line from first to second
                        iiContainer.currentWidget().removeItem(vpts[1].lines[0])
                        vpts[1].lines[0] = None
                    for k in range(i, len(vpts) - 1):
                        vpts[k] = vpts[k + 1]
                    del vpts[-1]
                    if len(vpts) < 2:
                        modelButton.setEnabled(False)
                        cProfButton.setEnabled(False)
                        meshButton.setEnabled(False)
                    print 'Number of markers is ', len(vpts)
                else:
                    # A marker is clicked on while no buttons are being held - it is selected and can be moved now
                    vptSel = True
                    vptCur = vpts[i]
                break # Exit loop if a marker has been clicked on

        if keysPress['ctrl']:
            for m in range(len(intLines)):
                cData = intLines[m].curve.getData()
                imin = curveDistance(e.pos().x(), e.pos().y(), cData)
                found = False
                if imin != -1:
                    # snap the cross to the line
                    iiContainer.currentWidget().removeItem(intLines[m])
                    del (intLines[m])
                    found = True
                if found:
                    break
        elif keysPress['alt']:
            #FIXME finish making
            for m in range(len(intLines)):
                cData = intLines[m].curve.getData()
                imin = curveDistance(e.pos().x(), e.pos().y(), cData)
                if imin != -1:
                    # If you clicked within one pixel of the line
                    # FIXME need to make this line into the path data
                    break
                    found = True


        if not vptSel and not keysPress['ctrl'] and not keysPress['shift']:
            # no you did not click on a point
            if not first:
                cx = e.pos().x() # in color coordinates
                cy = e.pos().y()
                if autoCorrectVpt.checkState() == 2:
                    cx, cy = centerVelocityStream(cx, cy)

            px, py = colorToProj(cx, cy) # color map to projected
            v0 = velocity.interp([px], [py], grid=False)
            dx, dy = colorToData(cx, cy)
            vpts.append(Marker(cx, cy, dx, dy, v0, iiContainer.currentWidget())) # in map coordinates x<10018, y< 17946

            x = int(np.floor(dx))
            y = int(np.floor(dy))
            txt = 'Point ' + str(len(vpts)-1) + \
            ':\n=================\n' + \
            'x: ' + str(cx) + '\n' +\
            'y: ' + str(cy) + '\n' +\
            'v: ' +      "{:.3f}".format(velocity.data[y][x]) + \
            '\nbed: ' +  "{:.3f}".format(bed.data[y][x]) + \
            '\nsurf: ' + "{:.3f}".format(surface.data[y][x]) + \
            '\nSMB: ' +  "{:.3f}".format(smb.data[y][x]*(1.0/1000.0)*(916.7/1000.0)) + \
            '\nSMB: ' +  "{:.3f}".format(smb.data[y][x]) + '\n\n'

            textOut.append(txt)
            iiContainer.currentWidget().addItem(vpts[-1].getCross()[0])
            iiContainer.currentWidget().addItem(vpts[-1].getCross()[1])
            # vpts[-1].setIntLine(calcProf(None))
            if len(vpts) > 1:
                if not modelButton.isEnabled():
                    modelButton.setEnabled(True)
                    cProfButton.setEnabled(True)
                    meshButton.setEnabled(True)
                xa = [vpts[-2].cx, vpts[-1].cx]
                ya = [vpts[-2].cy, vpts[-1].cy]
                vpts[-1].setLine(pg.PlotDataItem(xa, ya, connect='all', pen=skinnyBlackPlotPen), 0)
                vpts[-2].setLine(vpts[-1].lines[0], 1)
                iiContainer.currentWidget().addItem(vpts[-1].lines[0])  # ,pen=plotPen)
    else:
        #FIXME Why delete vptCur ?? I remember it caused a bug
        del vptCur
        vptSel = False


#
# def calcProf(e):
#     global integrateLine
#     x0p, y0p = colorToProj(vpts[-1].cx, vpts[-1].cy)
#     y0 = np.array([x0p, y0p])
#     t0, t1, dt = 0, 80, .1
#     r = ode(getProfile).set_integrator('zvode', method='bdf')
#     r.set_initial_value(y0, t0)
#     print "Printing integration points"
#     ox = [vpts[-1].cx]
#     oy = [vpts[-1].cy]
#     while r.successful() and r.t < t1:
#         # print(r.t+dt, r.integrate(r.t+dt))
#         ai = r.integrate(r.t+dt)
#         xi, yi = colorCoord(ai[0], ai[1])
#         # print 'xi, iy: ', xi, yi
#         ox.append(np.real(xi))
#         oy.append(np.real(yi))
#
#     # print 'ox, oy: ', ox, oy
#     integrateLine = pg.PlotDataItem(ox, oy, pen=whitePlotPen)
#     iiContainer.currentWidget().addItem(integrateLine)


def changeMap(index):
    '''
    Called when data-set drop down menu is changed.
    :param index:
    :return:
    '''

    # print iiContainer.currentWidget()
    vr = iiContainer.currentWidget().getPlotItem().getViewBox().viewRange()
    maps = [velocity, bed, surface, smb, thickness]
    global currentMap
    if index != currentMap:
        oldMap = currentMap
        currentMap = index
        # if not colormaps[maps[index].name]:
        #     print 'not colormap'
        #     maps[index].createColorMap()
        #     iiContainer.addWidget(maps[index].plotWidget)
        #     maps[index].imageItem.hoverEvent = mouseMoved

        iiContainer.setCurrentWidget(maps[index].plotWidget)
        maps[index].imageItem.hoverEvent = mouseMoved
        maps[index].imageItem.mouseClickEvent = mouseClick
        maps[index].plotWidget.getPlotItem().getViewBox().setRange(xRange=vr[0], yRange=vr[1], padding=0.0)
        for ln in intLines:
            maps[oldMap].plotWidget.removeItem(ln)
            maps[currentMap].plotWidget.addItem(ln)
        for pt in vpts:
            pt.plotWidget = maps[currentMap]
            maps[oldMap].plotWidget.removeItem(pt.cross[0])
            maps[oldMap].plotWidget.removeItem(pt.cross[1])
            maps[currentMap].plotWidget.addItem(pt.cross[0])
            maps[currentMap].plotWidget.addItem(pt.cross[1])
            if pt.lines[0]:
                maps[oldMap].plotWidget.removeItem(pt.lines[0])
                maps[currentMap].plotWidget.addItem(pt.lines[0])
            if pt.intLine:
                maps[oldMap].plotWidget.removeItem(pt.intLine)
                maps[currentMap].plotWidget.addItem(pt.intLine)



def calcBP():
    '''
    Calculate the data for bottom plot then populate the plot.
    :return:
    '''
    foo = StaticPlot(mw)
    # t0 = time.time()
    # # Empty the graph
    # velocity.pathPlotItem.clear()
    # surface.pathPlotItem.clear()
    # smb.pathPlotItem.clear()
    # bed.pathPlotItem.clear()
    # velocityWidth.pathPlotItem.clear()
    # thickness.pathPlotItem.clear()
    # if len(vpts) > 0:
    #     print 'Plotting...'
    #     # nbed, nsurf, nv, nsmb, nvelWidth, linePoints, graphX = interpolateData(True)
    #     interpolateData(False, 150)
    #     if velocityCheck.checkState() == 2:
    #         velocity.pathPlotItem.setData(velocity.distanceData          , velocity.pathData)
    #     if smbCheck.checkState() == 2:
    #         smb.pathPlotItem.setData(smb.distanceData                    , smb.pathData)
    #     if vWidthCheck.checkState() == 2:
    #         velocityWidth.pathPlotItem.setData(velocityWidth.distanceData, velocityWidth.pathData)
    #     if surfaceCheck.checkState() == 2:
    #         surface.pathPlotItem.setData(surface.distanceData            , surface.pathData)
    #     if bedCheck.checkState() == 2:
    #         bed.pathPlotItem.setData(bed.distanceData                    , bed.pathData)
    #     if thicknessCheck.checkState() == 2:
    #         thickness.pathPlotItem.setData(thickness.distanceData, thickness.pathData)
    #     if vWidthCheck.checkState() == 2:
    #         velocityWidth.pathPlotItem.setData(velocityWidth.distanceData, velocityWidth.pathData)
    #     pg.QtGui.QApplication.processEvents()
    #     runStaticPlot()

def mouseMoved(e):
    global vptSel, vptCur, integrateLine, currentMap
    if e.isExit() is False:
        if vptSel and len(intLines) != 0 : #integrateLine is not None:  # and integrateLine is not None:
            #
            # Checks to see if you are moving around a point and if it is near a curve
            # which it can be snapped to
            #
            i = 0
            while i < len(intLines):
                cData = intLines[i].curve.getData()
                imin = curveDistance(e.pos().x(), e.pos().y(), cData)
                if imin != -1:
                    # snap the cross to the line
                    vptCur.cx = cData[0][imin]
                    vptCur.cy = cData[1][imin]
                    vptCur.updateCross()
                    i = len(intLines) + 5
                else:
                    # just move the cross
                    vptCur.cx = e.pos().x()
                    vptCur.cy = e.pos().y()
                    vptCur.updateCross()
                if vptCur.lines[0] is not None:
                    vptCur.lines[0].setData([vptCur.lines[0].getData()[0][0], vptCur.cx],
                                            [vptCur.lines[0].getData()[1][0], vptCur.cy])  # previous line
                if vptCur.lines[1] is not None:
                    vptCur.lines[1].setData([vptCur.cx, vptCur.lines[1].getData()[0][1]],
                                            [vptCur.cy, vptCur.lines[1].getData()[1][1]])  # next line
                i += 1
        elif vptSel:
            vptCur.cx = e.pos().x()
            vptCur.cy = e.pos().y()
            vptCur.updateCross()
            if vptCur.lines[0] is not None:
                vptCur.lines[0].setData([vptCur.lines[0].getData()[0][0], vptCur.cx],
                                        [vptCur.lines[0].getData()[1][0], vptCur.cy])  # previous line
            if vptCur.lines[1] is not None:
                vptCur.lines[1].setData([vptCur.cx, vptCur.lines[1].getData()[0][1]],
                                        [vptCur.cy, vptCur.lines[1].getData()[1][1]])  # next

        x = int(np.floor(e.pos().x()))
        y = int(np.floor(e.pos().y()))
        if np.abs(x) <= map['cmap_x1'] and np.abs(y) <= map['cmap_y1']:
            mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y))# + '\n' + 'th: ' + str(thickness.data[y][x]))# + '\n' + 'oth: ' + str(oldthick.data[y][x]))



# def mouseMovedBP(evt):
#     global bpLegend, dataLen
#     print 'mouse moved bottom plot'
#     if botPlot:
#         print 'is not not not none'
#         pos = evt[0]  ## using signal proxy turns original arguments into a tuple
#         if bp.getPlotItem().sceneBoundingRect().contains(pos):
#             mousePoint = bp.getPlotItem().vb.mapSceneToView(pos)
#             index = int(mousePoint.x())
#             if index > 0 and index < len(velocity.pathData):
#                 velocity.legendItem.setText('Velocity' + str(velocity.pathData[index]))
#                 smb.legendItem.setText('SMB ' + str(smb.pathData[index]))
#                 velocityWidth.legendItem.setText('Velocity Width ' + str(velocityWidth.pathData[index]))
#                 surface.legendItem.setTxt('Surface ele ' + str(surface.pathData[index]))
#     # else:
#         # print 'It is NONE'


def arrows():
    ''''
    Plot velocity vectors to see the velocity field and its borders
    '''
    print 'arrows'
    print velocity.imageItem.getViewBox().viewRange() # [y0, y1] [x0, x1]
    rngx, rngy = velocity.imageItem.getViewBox().viewRange()
    print rngx
    print rngy
    print 'Starting arrows'
    for i in range(int(rngx[0]),int(rngx[1])):
        for j in range(int(rngy[0]), int(rngy[1])):
            vdir = [velocity.vx[j][i], -velocity.vy[j][i]]
            vmag = np.sqrt(velocity.vx[j][i]**2 +  velocity.vy[j][i]**2)
            # theta = np.arctan2(vdir[1], vdir[0])
            iiContainer.currentWidget().addItem(pg.PlotDataItem([(i+0.5),(i+0.5 + vdir[0]/(1.5*vmag))], [(j+0.5),(j+0.5 + vdir[1]/(1.5*vmag))], pen=whitePlotPen))
            iiContainer.currentWidget().addItem(pg.PlotDataItem([i+0.5], [j+0.5]), pen=(255,255,255), symbolBrush=(255,0,0), symbolPen='w')
    print 'finished arrows'


def calcProf(e):
    '''
    Calculates the line that shows velocity flow inwards (negative direction).
    :param e:
    :return:
    '''
    global integrateLine
    x0p, y0p = colorToProj(vpts[-1].cx, vpts[-1].cy)
    y0 = np.array([x0p, y0p])
    t0, t1, dt = 0, 80, .1
    r = ode(getProfile).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    ox = [vpts[-1].cx]
    oy = [vpts[-1].cy]
    while r.successful() and r.t < t1:
        # print(r.t+dt, r.integrate(r.t+dt))
        ai = r.integrate(r.t+dt)
        xi, yi = colorCoord(ai[0], ai[1])
        # print 'xi, iy: ', xi, yi
        ox.append(np.real(xi))
        oy.append(np.real(yi))
    integrateLine = pg.PlotDataItem(ox, oy, pen=whitePlotPen)
    iiContainer.currentWidget().addItem(integrateLine)



def regionIntLine(e):
    '''
    DEPRECATED -> WAS USED FOR TESTING
    Calculates and prints the integrated velocity path for several paths in a velocity stream.
    :param e:
    :return:
    '''
    xp0 = vpts[-1].cx
    yp0 = vpts[-1].cy
    xp1 = vpts[-2].cx
    yp1 = vpts[-2].cy
    xril, yril, xrir, yrir = calcVelWidth(xp0, yp0, xp1, yp1, False)  # for vpts[-1]

    theta = np.arctan2((yril - yrir), (xril - xrir))
    d = np.sqrt((yril - yrir) ** 2 + (xril - xrir) ** 2)
    dr = 10
    l = linspace(-d / 2, d / 2, dr, endpoint=True)

    rotMatrix = np.matrix([[np.cos(theta), -1 * np.sin(theta)], [np.sin(theta), np.cos(theta)]])

    lines = []

    for i in range(len(l)):
        '''
            Moves along velocity width line.
            Good spot for parallel programming.
        '''
        trot = rotMatrix * np.matrix([[l[i]], [0.0]])
        x0p, y0p = projCoord((xp1 + trot[0, 0]), (yp1 + trot[1, 0]))
        lines.append(intLine(x0p, y0p))
        iiContainer.currentWidget().addItem(lines[-1])

def ky(e):
    # 16777248 is shift
    # 16777249 is left ctrl
    #6 is press, 7 is release
    if e.type() == 6:
        if e.key() == 16777248:
            keysPress['shift'] = True
        elif e.key() == 16777249:
            keysPress['ctrl'] = True
        elif e.key() == 16777251:
            keysPress['alt'] = True
    else:
        keysPress['ctrl'] = False
        keysPress['shift'] = False
        keysPress['alt'] = False
