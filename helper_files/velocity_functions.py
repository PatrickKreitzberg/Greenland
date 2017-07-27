from data_functions import *
from scipy.integrate import ode
from classes.pltPoints import *

def _intMesh0(t,y):
    vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', [math.floor(y[0]) - 1, math.floor(y[1]) - 1],
                                          p1=[math.ceil(y[0]) + 1, math.ceil(y[1]) + 1], d2=velocity.vy)
    return np.array([t * (vxInterp([y[0]], [y[1]], grid=False)), t * (vyInterp([y[0]], [y[1]], grid=False))])

def intLine(x, y):
    '''
    Calculates and returns an integrated velocity flow path from the point clicked.  Calculates both inward
    and outward.
    :param e:
    :return:
    '''

    #
    # Removing the requirement to find the center of velocity stream to integrate from
    # This will be turned on if the check function works
    #

    x0p, y0p = projCoord(x, y)
    print 'x,y ', x0p, '  ----  ', y0p
    y0 = np.array([x0p, y0p])
    print y0
    t0, t1, dt = 0, 80, .1

    r = ode(getProfile).set_integrator('zvode', method='bdf')
    r.set_initial_value(y0, t0)
    print "Printing integration points"

    ox = []
    oy = []
    ivx, ivy = 5,5
    while r.successful() and sqrt(ivx**2 + ivy**2) > 5:
        ai = r.integrate(r.t + dt)
        xi, yi = mapCoord(ai[0], ai[1])
        vxInterp, vyInterp = getInterpolators(velocity.vx, 'velocity', ai[0], ai[1], d2=velocity.vy)
        ivx = vxInterp([ai[0]], [ai[1]], grid=False)
        ivy = vyInterp([ai[0]], [ai[1]], grid=False)
        # print 'xi, iy: ', xi, yi
        ox.append(np.real(xi))
        oy.append(np.real(yi))

    r2 = ode(_intMesh0).set_integrator('zvode', method='bdf')
    r2.set_initial_value(y0, 0)
    print "Printing integration points"
    # ox2 = []
    # oy2 = []
    ivx, ivy = 5, 5
    while r2.successful() and sqrt(ivx**2 + ivy**2) > 5:
        ai2 = r2.integrate(r2.t + dt)
        xi, yi = mapCoord(ai2[0], ai2[1])
        ivx = vxInterp([ai2[0]], [ai2[1]], grid=False)
        ivy = vyInterp([ai2[0]], [ai2[1]], grid=False)
        # print 'xi, iy: ', xi, yi
        ox.append(np.real(xi))
        oy.append(np.real(yi))
    return pg.PlotDataItem(ox, oy, pen=whitePlotPen)
    # iiContainer.currentWidget().addItem(pg.PlotDataItem(ox, oy, pen=plotPen2))
    # iiContainer.currentWidget().addItem(pg.PlotDataItem(ox2, oy2, pen=plotPen2))
