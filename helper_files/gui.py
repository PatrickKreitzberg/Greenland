import pyqtgraph as pg
from pyqtgraph.Qt import QtGui



#####################################################
####           BUILD GUI                         ####
#####################################################

print 'Building GUI'
#Main window
app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.setWindowTitle('GREENLAND')    # MAIN WINDOW
mw.setMinimumHeight(1000)
mw.setMinimumWidth(1200)
cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
mw.setCentralWidget(cw)
# l = QtGui.QGridLayout()            # CENTRAL WIDGET LAYOUT (layout of the central widget)
mainLayout = QtGui.QHBoxLayout()
cw.setLayout(mainLayout)
iiContainer = QtGui.QStackedWidget()    # STACKED WIDGET (inside the layout)
# l.setRowStretch(0,2)

#####################################################
####      SIDE WIDGET WITH BUTTONS               ####
#####################################################

bbW = QtGui.QWidget()
buttonBox = QtGui.QVBoxLayout()
bbW.setLayout(buttonBox)

mapList = QtGui.QComboBox()
showBedButton = QtGui.QPushButton('Show Bed Data')
showVelButton = QtGui.QPushButton('Show Velocity Data')
maps = ['Velocity', 'Bed', 'Surface', 'Smb']
mapList.addItems(maps)

clearButton      = QtGui.QPushButton('Clear Points')
calcWidthButton  = QtGui.QPushButton('Calculate Velocity Width')
intButton        = QtGui.QPushButton('Integrate')
cProfButton      = QtGui.QPushButton('Plot Path') #FIXME should automatically get profile
cRegionButton    = QtGui.QPushButton('Region')
cVelArrowsButton = QtGui.QPushButton('Arrows')
modelButton      = QtGui.QPushButton('Run Model')
textOut = QtGui.QTextBrowser()
maxWidth = 150

mapList.setMaximumWidth(maxWidth)
clearButton.setMaximumWidth(maxWidth)
clearButton.setMaximumWidth(maxWidth)
calcWidthButton.setMaximumWidth(maxWidth)
intButton.setMaximumWidth(maxWidth)
cProfButton.setMaximumWidth(maxWidth)
cRegionButton.setMaximumWidth(maxWidth)
cVelArrowsButton.setMaximumWidth(maxWidth)
modelButton.setMaximumWidth(maxWidth)
textOut.setMaximumWidth(maxWidth)

buttonBox.addWidget(mapList)
buttonBox.addWidget(clearButton)
buttonBox.addWidget(calcWidthButton)
buttonBox.addWidget(intButton)
buttonBox.addWidget(cProfButton)
buttonBox.addWidget(cRegionButton)
buttonBox.addWidget(cVelArrowsButton)
buttonBox.addWidget(modelButton)
buttonBox.addWidget(textOut)

print 'Done'


#####################################################
####         CREATE BOTTOM PLOT                  ####
#####################################################
bp = pg.PlotWidget()
bpLegend = bp.getPlotItem().addLegend()

iiContainer.setMinimumHeight(mw.height()*(2/3))
lsw = QtGui.QWidget()
leftSide = QtGui.QVBoxLayout()
lsw.setLayout(leftSide)
leftSide.addWidget(iiContainer,2)
leftSide.addWidget(bp,1)

mainLayout.addWidget(lsw)
mainLayout.addWidget(bbW)
bbW.setMaximumWidth(maxWidth+12)

mw.show()
