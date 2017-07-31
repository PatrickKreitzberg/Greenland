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

buttonBoxWidget = QtGui.QWidget()
buttonBox = QtGui.QVBoxLayout()
buttonBoxWidget.setLayout(buttonBox)

mapList = QtGui.QComboBox()
maps = ['Velocity', 'Bed', 'Surface', 'Smb', 'Thickness']#, 'OldThickness']
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
hiResButton      = QtGui.QPushButton('High-res interpolators')
mouseCoordinates = QtGui.QLabel('x:\ty:')
textOut = QtGui.QTextBrowser()

maxWidth = 300

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
hiResButton.setMaximumWidth(maxWidth)
textOut.setMaximumWidth(maxWidth)


buttonBox.addWidget(mapList)
buttonBox.addWidget(autoCorrectVpt)
buttonBox.addWidget(clearButton)
buttonBox.addWidget(calcWidthButton)
buttonBox.addWidget(intButton)
buttonBox.addWidget(cProfButton)
# buttonBox.addWidget(cRegionButton)
# buttonBox.addWidget(cVelArrowsButton)
buttonBox.addWidget(modelButton)
buttonBox.addWidget(hiResButton)
buttonBox.addWidget(mouseCoordinates)


# LABELS AND LINE EDITS
time_widget    = QtGui.QWidget()

time_container = QtGui.QGridLayout()
time_widget.setLayout(time_container)

model_res_label    = QtGui.QLabel('Sptl res:')
model_res_lineEdit = QtGui.QLineEdit('150')

t_end_label    = QtGui.QLabel('t_end:')
t_end_lineEdit = QtGui.QLineEdit('20000')

t_step_label    = QtGui.QLabel('t_step:')
t_step_lineEdit = QtGui.QLineEdit('10')
t_current = QtGui.QLabel('Current year: ')

time_container.addWidget(model_res_label,    0, 0)
time_container.addWidget(model_res_lineEdit, 0, 1)
time_container.addWidget(t_end_label,        1, 0)
time_container.addWidget(t_end_lineEdit,     1, 1)
time_container.addWidget(t_step_label,       2, 0)
time_container.addWidget(t_step_lineEdit,    2, 1)
time_container.addWidget(t_current,          3, 0, 1, 2)

buttonBox.addWidget(time_widget)
buttonBox.addWidget(textOut)




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
dataCheckContainer = QtGui.QWidget()
dCClayout = QtGui.QHBoxLayout()


dataCheckContainer.setLayout(dCClayout)

velocityCheck   = QtGui.QCheckBox('Velocity')
vWidthCheck     = QtGui.QCheckBox('Velocity Width')
smbCheck        = QtGui.QCheckBox('SMB')
surfaceCheck    = QtGui.QCheckBox('Surface')
bedCheck        = QtGui.QCheckBox('Bed')


# velocityCheck.setMaximumWidth(80)
# vWidthCheck.setMaximumWidth(80)
# smbCheck.setMaximumWidth(80)
# surfaceCheck.setMaximumWidth(80)
# bedCheck.setMaximumWidth(80)


checkBoxW = QtGui.QWidget()
checkBLayout = QtGui.QVBoxLayout()
checkBoxW.setLayout(checkBLayout)

velocityCheck.setTristate(False)
vWidthCheck.setTristate(False)
smbCheck.setTristate(False)
surfaceCheck.setTristate(False)
bedCheck.setTristate(False)

velocityCheck.setCheckState(2)
vWidthCheck.setCheckState(2)
smbCheck.setCheckState(2)
surfaceCheck.setCheckState(2)
bedCheck.setCheckState(2)

checkBLayout.addWidget(QtGui.QLabel('Plot Checked Data:'))
checkBLayout.addWidget(velocityCheck)
checkBLayout.addWidget(vWidthCheck)
checkBLayout.addWidget(smbCheck)
checkBLayout.addWidget(surfaceCheck)
checkBLayout.addWidget(bedCheck)
checkBLayout.setSpacing(0)

buttonBox.addWidget(checkBoxW)
dataCheckContainer.setMaximumWidth(500)


mainLayout.addWidget(lsw)
mainLayout.addWidget(buttonBoxWidget)
buttonBoxWidget.setMaximumWidth(maxWidth + 12)

mw.show()


#####################################################
####         MISCELLANEOUS                       ####
#####################################################
currentMap = 0  # selects which data map to show [velocity, bed, surface]
vptSel = False
shift = False

print 'Done building gui'