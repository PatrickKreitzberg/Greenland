import pyqtgraph as pg
from ..pens import *
from ..math_functions import *

class vpt:
    def __init__(self, cx, cy, dx, dy, velocity, plotWidget):
        self.plotWidget = plotWidget
        self.cx, self.cy = cx, cy  # color coordinates
        self.dx, self.dy = dx, dy  # data coordinates
        self.px, self.py = colorToProj(self.cx, self.cy)
        self.v = velocity
        self.pen = blackPlotPen
        self.pen.setWidth(2)
        c = 3
        xV0 = [self.cx - c, self.cx + c]
        xV1 = [self.cx - c, self.cx + c]
        yV0 = [self.cy - c, self.cy + c]
        yV1 = [self.cy + c, self.cy - c]
        self.cross = [pg.PlotDataItem(xV0, yV0, connect='all', pen=self.pen), pg.PlotDataItem(xV1, yV1, connect='all', pen=self.pen)]
        self.lines = [None] * 2
        self.intLine = None

    def __del__(self):
        if self.lines[0]:
            self.plotWidget.removeItem(self.lines[0])
        if self.lines[1]:
            self.plotWidget.removeItem(self.lines[1])
        self.plotWidget.removeItem(self.cross[0])
        self.plotWidget.removeItem(self.cross[1])
        if self.intLine:
            self.plotWidget.removeItem(self.intLine)

    def updateCross(self):
        c = 3
        xV0 = [self.cx - c, self.cx + c]
        xV1 = [self.cx - c, self.cx + c]
        yV0 = [self.cy - c, self.cy + c]
        yV1 = [self.cy + c, self.cy - c]
        self.dx, self.dy = colorToData(self.cx, self.cy)
        self.cross[0].setData([self.cx - c, self.cx + c], [self.cy - c, self.cy + c], connect='all', pen=self.pen)
        self.cross[1].setData([self.cx - c, self.cx + c], [self.cy + c, self.cy - c], connect='all', pen=self.pen)
        self.cross[0].updateItems()
        self.cross[1].updateItems()

    def setPlotWidget(self, pw):
        self.plotWidget = pw

    def checkClicked(self, pos):
        # if self.cross[0].curve.mouseShape().contains(pos) or self.cross[1].curve.mouseShape().contains(pos):
        c = 3
        if self.cx - c <= pos.x() <= self.cx + c and self.cy - c <= pos.y() <= self.cy + c:
            return True
        else:
            return False

    def setIntLine(self, ln):
        self.intLine = ln

    def getIntLine(self):
        return self.intLine

    def setLine(self, line, i):#, index):
        # 0 connects to the previous Marker, 1 connects to the second
        self.lines[i] = line

    def getLine(self): #, index):
        return self.lines # s[index]

    def getCross(self):
        return self.cross[0], self.cross[1]

    def setV(self, v):
        self.v = v

    def getV(self):
        return self.v