import pyqtgraph as pg
from pyqtgraph.Qt import QtGui

'''
This is to build all of the gui components for the program.  This means the application and
its main window along with the plots and buttons.

The purpose of this file is to reduce the clutter of the main file.
'''


#####################################################
####           BUILD GUI                         ####
#####################################################

print 'Building GUI'
#Main window
app = QtGui.QApplication([])
mw = QtGui.QMainWindow()
mw.setWindowTitle('GREENLAND')
mw.setMinimumHeight(1000)
mw.setMinimumWidth(1200)
cw = QtGui.QWidget()            # GENERIC WIDGET AS CENTRAL WIDGET (inside main window)
mw.setCentralWidget(cw)
mainLayout = QtGui.QHBoxLayout()
cw.setLayout(mainLayout)
iiContainer = QtGui.QStackedWidget()    # STACKED WIDGET (inside the layout)

#####################################################
####      SIDE WIDGET WITH BUTTONS               ####
#####################################################

bbW = QtGui.QWidget()
buttonBox = QtGui.QVBoxLayout()
bbW.setLayout(buttonBox)

mapList = QtGui.QComboBox()
maps = ['Velocity', 'Bed', 'Surface', 'Smb']
mapList.addItems(maps)

autoCorrectVpt   = QtGui.QCheckBox('Auto-correct vpt pos.')
autoCorrectVpt.setTristate(False)
autoCorrectVpt.setCheckState(2)
clearButton      = QtGui.QPushButton('Clear Points')
calcWidthButton  = QtGui.QPushButton('Calculate Velocity Width')
intButton        = QtGui.QPushButton('Integrate')
cProfButton      = QtGui.QPushButton('Plot Path') #FIXME should automatically get profile
cRegionButton    = QtGui.QPushButton('Region')
cVelArrowsButton = QtGui.QPushButton('Arrows')
modelButton      = QtGui.QPushButton('Run Model')



textOut = QtGui.QTextBrowser()

maxWidth = 150 #Maximum width for all the buttons on the right side panel

mapList.setMaximumWidth(maxWidth)
autoCorrectVpt.setMaximumWidth(maxWidth)
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
buttonBox.addWidget(autoCorrectVpt)
buttonBox.addWidget(clearButton)
buttonBox.addWidget(calcWidthButton)
buttonBox.addWidget(intButton)
buttonBox.addWidget(cProfButton)
buttonBox.addWidget(cRegionButton)
buttonBox.addWidget(cVelArrowsButton)
buttonBox.addWidget(modelButton)


time_widget    = QtGui.QWidget()
time_container = QtGui.QGridLayout()
time_widget.setLayout(time_container)

t_end_label    = QtGui.QLabel('t_end:')
t_end_lineEdit = QtGui.QLineEdit('20000')

t_step_label    = QtGui.QLabel('t_step:')
t_step_lineEdit = QtGui.QLineEdit('10')

time_container.addWidget(t_end_label, 0, 0)
time_container.addWidget(t_end_lineEdit, 0, 1)
time_container.addWidget(t_step_label, 1, 0)
time_container.addWidget(t_step_lineEdit, 1, 1)

buttonBox.addWidget(time_widget)
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


#####################################################
####         MISCELLANEOUS                       ####
#####################################################
currentMap = 0  # selects which data map to show [velocity, bed, surface]
vptSel = False
mouseCoordinates = QtGui.QLabel('x:\ty:')
shift = False