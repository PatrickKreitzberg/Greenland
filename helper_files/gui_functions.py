from mathFunctions import *
from pens import *
from gui import *

#########################################
####          CURRENTLY NOT IN USE   ####
#########################################



def calcProf(e):
    global integrateLine
    x0p, y0p = projCoord(vpts[-1].getX(), vpts[-1].getY())
    y0 = np.array([x0p, y0p])
    t0, t1, dt = 0, 80, .1
    r = ode(getProfile).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    print "Printing integration points"
    ox = [vpts[-1].getX()]
    oy = [vpts[-1].getY()]
    while r.successful() and r.t < t1:
        # print(r.t+dt, r.integrate(r.t+dt))
        ai = r.integrate(r.t+dt)
        xi, yi = mapCoord(ai[0], ai[1])
        # print 'xi, iy: ', xi, yi
        ox.append(np.real(xi))
        oy.append(np.real(yi))

    print 'ox, oy: ', ox, oy
    integrateLine = pg.PlotDataItem(ox,oy,pen=plotPen2)
    iiContainer.currentWidget().addItem(integrateLine)