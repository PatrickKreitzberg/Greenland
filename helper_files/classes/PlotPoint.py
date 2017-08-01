import pyqtgraph as pg

class vpt:

    def __init__(self,x, y, velocity, velW):
        self.velW = velW
        self.x = x
        self.y = y
        self.v = velocity
        self.pen = pg.mkPen(color=(255, 255, 255))
        self.pen.setWidth(2)
        c = 3
        xV0 = [x - c, x + c]
        xV1 = [x - c, x + c]
        yV0 = [y - c, y + c]
        yV1 = [y + c, y - c]
        self.cross = [pg.PlotDataItem(xV0, yV0, connect='all', pen=self.pen), pg.PlotDataItem(xV1, yV1, connect='all', pen=self.pen)]
        self.line = None #[None]*2
        self.intLine = None

    def __del__(self):
        self.velW.removeItem(self.line)
        self.velW.removeItem(self.cross[0])
        self.velW.removeItem(self.cross[1])
        self.velW.removeItem(self.intLine)

    def updateCross(self, x, y):
        self.x = x
        self.y = y
        c = 3
        xV0 = [self.x - c, self.x + c]
        xV1 = [self.x - c, self.x + c]
        yV0 = [self.y - c, self.y + c]
        yV1 = [self.y + c, self.y - c]
        self.cross[0].setData(xV0, yV0, connect='all', pen=self.pen)
        self.cross[1].setData(xV1, yV1, connect='all', pen=self.pen)
        self.cross[0].updateItems()
        self.cross[1].updateItems()

    def checkClicked(self, pos):
        if self.cross[0].curve.mouseShape().contains(pos) or self.cross[1].curve.mouseShape().contains(pos):
            return True
        else:
            return False


    def setIntLine(self, ln):
        self.intLine = ln

    def getIntLine(self):
        return self.intLine

    def setLine(self, line):#, index):
        self.line = line
        #0 connects to the previous vpt, 1 connects to the second
        #self.lines[index] = line

    def getLine(self): #, index):
        return self.line # s[index]

    def getCross(self):
        return self.cross[0], self.cross[1]

    def setV(self, v):
        self.v = v

    def getV(self):
        return self.v

    def setPos(self,x,y):
        self.x = x
        self.y = y

    def getPos(self):
        return self.x, self.y

    def getX(self):
        return self.x

    def getY(self):
        return self.y


class bs():
    def __init__(self,x,y):
        self.x = x
        self.y = y
