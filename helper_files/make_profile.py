import numpy as np
import matplotlib.pyplot as plt
from pylab import imshow,log10,flipud,sqrt,show,linspace,array,argmax
from scipy.interpolate import interp2d
from scipy.io import loadmat
import fenics as fc

#plt.matplotlib.rcParams['backend'] = 'TkAgg' 
# Import and prepare data:
b = loadmat('data/bamber.mat')
r = loadmat('data/rignot_proj.mat')
# extents
#extents =  (b['map_western_edge'][0][0],b['map_eastern_edge'][0][0],\
#            b['map_southern_edge'][0][0],b['map_northern_edge'][0][0])

RESOLUTION = 1000 # Desired resolution in meter.
THICKLIMIT = 10.  # Ice is never less than this thick

extents = (-987578,793599,-3274259,-746617)
xint = linspace(extents[0],extents[1],b['nx'][0][0])
yint = linspace(extents[2],extents[3],b['ny'][0][0])

Binterp = interp2d(xint,yint,b['B'])
Sinterp = interp2d(xint,yint,b['S'])

vxinterp = interp2d(xint,yint,r['vx'])
vyinterp = interp2d(xint,yint,r['vy'])

def dxdt(x):
  return array([-vxinterp(x[0],x[1]),-vyinterp(x[0],x[1])])

def get_margins(x,y):
  vx = vxinterp(x,y)
  vy = vyinterp(x,y)
  n = (-vx,vy) / sqrt(vx**2 + vy**2)
  d = 50e3
  nsamp = 500
  line = array([x+d*linspace(-n[0] , n[0], nsamp),y+d*linspace(-n[1],n[1],nsamp)])
  vxd = vxinterp(line[0],line[1])
  vyd = vyinterp(line[0],line[1])
  v = sqrt(vxd**2+vyd**2)
  dv = v[1:] - v[:-1]
  il = argmax(dv[:nsamp/2])
  ir = argmax(dv[nsamp/2:])
  print il,ir

  return line[0,il], line[1,il], line[0,nsamp/2+ir],line[1,nsamp/2+ir]

# full res speed
v = sqrt(r['vx']**2 + r['vy']**2)
v = flipud(v)
v[v<=.001] = .001

line_points_x = []
line_points_y = []
l_margin_line_x =[]
l_margin_line_y =[]
r_margin_line_x =[]
r_margin_line_y =[]

# Two subplots, the axes array is 1-d
f, axarr = plt.subplots(1, 2, sharex=True, sharey=True)

axarr[0].imshow(log10(v),extent = extents,interpolation='nearest')
axarr[0].set_title("Use tool to zoom on this image.")

axarr[1].imshow(log10(v),extent = extents,interpolation='nearest')
axarr[1].set_title("Mouse click to create centerline points.")
line = axarr[1].plot(line_points_x,line_points_y,'wo-',lw=5,ms=8)[0]
margin_line_l = axarr[1].plot(line_points_x,line_points_y,'yo-')[0]
margin_line_r = axarr[1].plot(line_points_x,line_points_y,'yo-')[0]

def smooth(S,wl):
    window_len = wl
    s=np.r_[S[wl-1:0:-1],S,S[-1:-wl:-1]]
    w=np.ones(wl,'d')

    y= np.convolve(w/w.sum(),s,mode='valid')
    return y[(int(wl/2)-1):-(int(wl/2)+1)]


def onButtonPress(event):
    global line_points_x,line_points_y
    if event.inaxes is not None:
        ax = event.inaxes
        if ax == axarr[1] and event.button == 1:
            # print('button=%d, x=%d, y=%d, xdata=%f, ydata=%f' \
            #     %(event.button, event.x, event.y, event.xdata, event.ydata))
            line_points_x.append(event.xdata)
            line_points_y.append(event.ydata)
            line.set_xdata(line_points_x)
            line.set_ydata(line_points_y)
            f.canvas.draw()

            # This was supposed to get the margins and plot them, but it isn't working.
            """
            xml,yml,xmr,ymr = get_margins(event.xdata,event.ydata)
            print event.xdata,event.ydata,xml,yml,xmr,ymr
            l_margin_line_x.append(xml)
            l_margin_line_y.append(yml)
            r_margin_line_x.append(xmr)
            r_margin_line_y.append(ymr)
            margin_line_l.set_xdata(l_margin_line_x)
            margin_line_l.set_ydata(l_margin_line_y)
            margin_line_r.set_xdata(r_margin_line_x)
            margin_line_r.set_ydata(r_margin_line_y)
            """

            f.canvas.draw()
        elif ax == axarr[1] and event.button == 3:
            # Collect data
            B = []
            S = []
            xf = []
            yf = []
            """ What a fucking abortion, trying to keep the spacing of interpolated points even
              while moving through the individual points in the profile. 
              This works within a meter or so, meaning that the difference between RESOLUTION and what is actually
              done is at most a meter.
            """
            for i in range(1,len(line_points_x)):
                dx = line_points_x[i]-line_points_x[i-1]
                dy = line_points_y[i]-line_points_y[i-1]
                d = sqrt(dx**2+dy**2)
                xhat = dx / d
                yhat = dy / d
                total_distance = 0.
                if i == 1:
                    p = 0
                    remainder = 0
                else:
                    p = 1
                while (p * RESOLUTION + (RESOLUTION - remainder)) < d:
                    if p <= 1:
                        xf.append(line_points_x[i-1] + xhat * p * (RESOLUTION - remainder))
                        yf.append(line_points_y[i-1] + yhat * p * (RESOLUTION - remainder))
                    else:
                        xf.append(xf[-1] + xhat * RESOLUTION)
                        yf.append(yf[-1] + yhat * RESOLUTION)
                    p+=1
                remainder = d - ((p-2) * RESOLUTION + (RESOLUTION - remainder))
            for x,y in zip(xf,yf):
                B.append(Binterp(x,y))
                S.append(Sinterp(x,y))
            print x,y,B[-1],S[-1]
        B = array(B)[:,0]
        S = array(S)[:,0]

        # Smooth the surface
        #S = smooth(S,15)
        #B = smooth(B,15)
        # Plot data
        ff,ax = plt.subplots(1)
        ax.plot(B,lw=4)
        ax.plot(S,lw=4)
        plt.show()
        # Write data
        N = len(B)
        mesh = fc.IntervalMesh(N-1,0,RESOLUTION*(N-1))
        V = fc.FunctionSpace(mesh,"CG",1)
        hfile = fc.HDF5File(mesh.mpi_comm(),"latest_profile.h5","w")
        Bed     = fc.Function(V,name="Bed")
        Surface = fc.Function(V, name="Surface")
        H = S-B
        S[H<=THICKLIMIT] = B[H<=THICKLIMIT]
        Bed.vector()[:]     =  B
        Surface.vector()[:] =  S
        hfile.write(Bed.vector(),"/bed")
        hfile.write(Surface.vector(),"/surface")
        hfile.write(mesh,"/mesh")
        hfile.close()
              
f.canvas.mpl_connect('button_press_event', onButtonPress)
plt.show()
