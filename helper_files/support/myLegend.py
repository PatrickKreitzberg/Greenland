from pyqtgraph import LegendItem, PlotDataItem, ScatterPlotItem, LabelItem
from pyqtgraph.Qt import QtCore, QtGui
from pyqtgraph.graphicsItems.ScatterPlotItem import drawSymbol
from pyqtgraph.graphicsItems.LegendItem import ItemSample
import pyqtgraph.functions as fn

class myLegend(LegendItem):
    def __init__(self, size=None, offset=None):
        LegendItem.__init__(self, size, offset)

    def paint(self, p, *args):
        p.setPen(fn.mkPen(0,0,0)) # outline
        p.setBrush(fn.mkBrush(255,255,255))   # background
        p.drawRect(self.boundingRect())

    def addItem(self, item, name):
        """
        Add a new entry to the legend.

        ==============  ========================================================
        **Arguments:**
        item            A PlotDataItem from which the line and point style
                        of the item will be determined or an instance of
                        ItemSample (or a subclass), allowing the item display
                        to be customized.
        title           The title to display for this item. Simple HTML allowed.
        ==============  ========================================================
        """
        label = LabelItem(name, color=(0,0,0))
        if isinstance(item, ItemSample):
            sample = item
        else:
            sample = ItemSample(item)
        row = self.layout.rowCount()
        self.items.append((sample, label))
        self.layout.addItem(sample, row, 0)
        self.layout.addItem(label, row, 1)
        self.updateSize()
