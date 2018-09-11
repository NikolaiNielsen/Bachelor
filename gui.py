import sys
# from matplotlib.figure import Figure
from cmp import *

from PyQt5 import QtWidgets
# from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class PlotWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QHBoxLayout(self._main)

        static_fig, static_ax = Lattice(returns=True, plots=False)
        static_canvas = FigureCanvas(static_fig)
        layout.addWidget(static_canvas)
        self.addToolBar(NavigationToolbar(static_canvas, self))
        static_ax.mouse_init()


class OptionsWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.title = "Options window"
        self.left = 10
        self.top = 10
        self.width = 320
        self.height = 200
        self.lattices = ['simple cubic', 'fcc', 'bcc']
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.initUI()

    def initUI(self):

        button = QtWidgets.QPushButton("Show plot", self)
        button.setToolTip('this is an example button')
        button.move(0, 0)
        button.clicked.connect(self.on_click)

        self.latticeChoice = QtWidgets.QLabel('Simple cubic', self)
        self.latticeChoice.move(0, 50)

        latticeChooser = QtWidgets.QComboBox(self)
        latticeChooser.addItem('hello')
        latticeChooser.addItem('There')
        latticeChooser.move(0, 100)
        latticeChooser.activated[str].connect(self.change_lattice)

        self.dialog = PlotWindow()

    def on_click(self):
        self.dialog.show()

    def change_lattice(self, text):
        print(text)


def main():
    qapp = QtWidgets.QApplication(sys.argv)
    app = OptionsWindow()
    app.show()
    sys.exit(qapp.exec_())


if __name__ == "__main__":
    main()
