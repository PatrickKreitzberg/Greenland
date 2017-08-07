from pyqtgraph.Qt import QtGui


class Instructions:
    def __init__(self, parent):
        self.mw = QtGui.QMainWindow(parent)
        self.mw.setWindowTitle('Instructions')
        self.widget = QtGui.QWidget()
        self.mw.setCentralWidget(self.widget)
        self.lay = QtGui.QVBoxLayout()
        self.txt = QtGui.QTextBrowser()
        self.widget.setLayout(self.lay)
        self.lay.addWidget(self.txt)
        self.txt.setText('click on map to set a point\n'
                         'shift+click on a point to show integrated velocity line.\n'
                         'ctrl+click integrated line to delete it.\n')
        self.mw.show()
