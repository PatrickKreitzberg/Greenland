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
    if len(markers) == 0:
        first = True
        cx = e.pos().x()
        cy = e.pos().y()
        if autoCorrectVpt.checkState() == 2:
            cx, cy = centerVelocityStream(cx, cy)
    if not vptSel:
        # If there is not a marker selected already
        keyClicked = False
        for i in range(len(markers)):
            if markers[i].checkClicked(e.pos()):
                keyClicked = True
                # See if you clicked on an already existing point
                if keysPress['shift']:
                    x0p, y0p = colorToProj(markers[i].cx, markers[i].cy)
                    y0 = np.array([x0p, y0p])
                    t0, t1, dt = 0, 80, .1
                    r = ode(getProfile).set_integrator('zvode', method='bdf')
                    r.set_initial_value(y0, t0)
                    ox = [markers[i].cx]
                    oy = [markers[i].cy]
                    try:
                        segLength = float(intResInput.text())
                        while r.successful() and r.t < t1:
                            ai = r.integrate(r.t + dt)
                            xi, yi = colorCoord(ai[0], ai[1])
                            if np.sqrt((xi-ox[-1])**2 + (yi-oy[-1])**2) > segLength/150:
                                ox.append(np.real(xi))
                                oy.append(np.real(yi))
                        print 'ox, oy', len(ox), len(oy)
                        print ox
                        print oy

                        intLines.append(pg.PlotDataItem(ox, oy, pen=whitePlotPen))
                        iiContainer.currentWidget().addItem(intLines[-1])
                    except ValueError:
                        textOut.append(('\nMust enter valid number for integration line resolution!'))

                elif keysPress['ctrl']:
                    # if you ctrl+click a maker it is deleted
                    if i > 0:
                        if i + 1 < len(markers):
                            # connect line from previous node to next point
                            markers[i - 1].lines[1] = pg.PlotDataItem([markers[i - 1].cx, markers[i + 1].cx], [markers[i - 1].cy, markers[i + 1].cy], connect='all', pen=skinnyBlackPlotPen)
                            markers[i + 1].lines[0] = markers[i - 1].lines[1]
                            iiContainer.currentWidget().addItem(markers[i - 1].lines[1])
                            pg.QtGui.QApplication.processEvents()
                        else:  # delete line because there is no next point
                            iiContainer.currentWidget().removeItem(markers[i - 1].lines[1])
                            markers[i - 1].lines[1] = None
                    elif i == 0 and len(markers) > 1:  # point is the first so delete line from first to second
                        iiContainer.currentWidget().removeItem(markers[1].lines[0])
                        markers[1].lines[0] = None
                    for k in range(i, len(markers) - 1):
                        markers[k] = markers[k + 1]
                    del markers[-1]
                    if len(markers) < 2:
                        modelButton.setEnabled(False)
                        cProfButton.setEnabled(False)
                        meshButton.setEnabled(False)
                    print 'Number of markers is ', len(markers)
                else:
                    # A marker is clicked on while no buttons are being held - it is selected and can be moved now
                    vptSel = True
                    vptCur = markers[i]
                break # Exit loop if a marker has been clicked on
        for ln in intLines:
            if sqrt((e.pos().x() - ln.curve.getData()[0][-1])**2 + (e.pos().y() - ln.curve.getData()[1][-1])**2) < 2:
                print 'Clicked end of line'
                globalConstants['moveLine'] = True

        print 'out of forrest'
        # if keysPress['alt']:
        #     print 'Entering with alt pressed'

        if not keyClicked and keysPress['ctrl']:
            print 'control clicked'
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

        elif not keyClicked and keysPress['shift']: # elif not keyClicked and keysPress['alt']:
            print 'Alt clicked'
            #FIXME finish making
            for m in range(len(intLines)):
                cData = intLines[m].curve.getData()
                imin = curveDistance(e.pos().x(), e.pos().y(), cData)
                if imin != -1:
                    # If you clicked within one pixel of the line
                    # FIXME need to make this line into the path data
                    print 'shift clicked line'

                    globalConstants['isPathIntLine'] = True
                    globalConstants['lineData'] = intLines[m].curve.getData()
                    globalConstants['lineIndex'] = m
                    globalConstants['minIndex'] = len(globalConstants['lineData']) - 1
                    globalConstants['minValue'] = 99#globalConstants['lineData'][-1]
                    intLines[m].setPen(purplePlotPen)

                    # Delete markers
                    # del markers[:]
                    # textOut.setText('')
                    # modelButton.setEnabled(False)
                    # cProfButton.setEnabled(False)
                    # velocity.pathPlotItem.clear()
                    # surface.pathPlotItem.clear()
                    # smb.pathPlotItem.clear()
                    # bed.pathPlotItem.clear()
                    # thickness.pathPlotItem.clear()
                    #
                    # # Set markers to int line
                    # for x,y in zip(cData[0], cData[1]):
                    #     px, py = colorToProj(x, y)
                    #     dx, dy = colorToData(x,y)
                    #     markers.append(Marker(x, y, dx, dy, velocity.interp([px], [py], grid=False), iiContainer.currentWidget()))
                    #     iiContainer.currentWidget().addItem(markers[-1].getCross()[0])
                    #     iiContainer.currentWidget().addItem(markers[-1].getCross()[1])
                    break

        elif not keyClicked and not keysPress['shift']:
            # no you did not click on a point
            if not first:
                cx = e.pos().x() # in color coordinates
                cy = e.pos().y()
                if autoCorrectVpt.checkState() == 2:
                    cx, cy = centerVelocityStream(cx, cy)

            px, py = colorToProj(cx, cy) # color map to projected
            v0 = velocity.interp([px], [py], grid=False)
            dx, dy = colorToData(cx, cy)
            markers.append(Marker(cx, cy, dx, dy, v0, iiContainer.currentWidget())) # in map coordinates x<10018, y< 17946

            x = int(np.floor(dx))
            y = int(np.floor(dy))
            txt = 'Point ' + str(len(markers) - 1) + \
            ':\n=================\n' + \
            'x: ' + str(cx) + '\n' +\
            'y: ' + str(cy) + '\n' +\
            'v: ' +      "{:.3f}".format(velocity.data[y][x]) + \
            '\nbed: ' +  "{:.3f}".format(bed.data[y][x]) + \
            '\nsurf: ' + "{:.3f}".format(surface.data[y][x]) + \
            '\nSMB: ' +  "{:.3f}".format(smb.data[y][x]*(1.0/1000.0)*(916.7/1000.0)) + \
            '\nSMB: ' +  "{:.3f}".format(smb.data[y][x]) + '\n\n'

            textOut.append(txt)
            iiContainer.currentWidget().addItem(markers[-1].getCross()[0])
            iiContainer.currentWidget().addItem(markers[-1].getCross()[1])
            # vpts[-1].setIntLine(calcProf(None))
            if len(markers) > 1:
                if not modelButton.isEnabled():
                    modelButton.setEnabled(True)
                    cProfButton.setEnabled(True)
                    meshButton.setEnabled(True)
                xa = [markers[-2].cx, markers[-1].cx]
                ya = [markers[-2].cy, markers[-1].cy]
                markers[-1].setLine(pg.PlotDataItem(xa, ya, connect='all', pen=skinnyBlackPlotPen), 0)
                markers[-2].setLine(markers[-1].lines[0], 1)
                iiContainer.currentWidget().addItem(markers[-1].lines[0])  # ,pen=plotPen)
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
        for pt in markers:
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



def openStaticPlotter():
    '''
    Calculate the data for bottom plot then populate the plot.
    :return:
    '''
    foo = StaticPlot(mw)


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

        elif globalConstants['moveLine']:
            mi = len(globalConstants['lineData']) - 1
            mv = globalConstants['minValue'] = globalConstants['lineData'][-1]
            for i in range(len(globalConstants['lineData'][0])):
                # globalConstants['minIndex'] = len(globalConstants['lineData']) - 1
                # globalConstants['minValue'] = globalConstants['lineData'][-1]

                if sqrt((globalConstants['lineData'][0][globalConstants['minIndex']] - e.pos().x())**2 +
                        (globalConstants['lineData'][1][globalConstants['minIndex']] - e.pos().y())**2) < mv:
                    mi = i
                    mv = sqrt((globalConstants['lineData'][0][globalConstants['minIndex']] - e.pos().x())**2 +
                         (globalConstants['lineData'][1][globalConstants['minIndex']] - e.pos().y())**2)
                intLines[globalConstants['lineIndex']].curve.setData(globalConstants['lineData'][0][:mi], globalConstants['lineData'][0][:mi])


        x = int(np.floor(e.pos().x()))
        y = int(np.floor(e.pos().y()))
        if np.abs(x) <= map['cmap_x1'] and np.abs(y) <= map['cmap_y1']:
            mouseCoordinates.setText('x: ' + str(x) + '\ty: ' + str(y))# + '\n' + 'th: ' + str(thickness.data[y][x]))# + '\n' + 'oth: ' + str(oldthick.data[y][x]))




def calcProf(e):
    '''
    Calculates the line that shows velocity flow inwards (negative direction).
    :param e:
    :return:
    '''
    global integrateLine
    print 'calcProf'
    x0p, y0p = colorToProj(markers[-1].cx, markers[-1].cy)
    y0 = np.array([x0p, y0p])
    t0, t1, dt = 0, 80, .1
    r = ode(getProfile).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    ox = [markers[-1].cx]
    oy = [markers[-1].cy]
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
    xp0 = markers[-1].cx
    yp0 = markers[-1].cy
    xp1 = markers[-2].cx
    yp1 = markers[-2].cy
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
            print 'pressed alt'
    else:
        keysPress['ctrl'] = False
        keysPress['shift'] = False
        keysPress['alt'] = False
