import sys
# from matplotlib.figure import Figure
from cmp import *

from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class ApplicationWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QVBoxLayout(self._main)

        static_fig, static_ax = Lattice(returns=True, plots=False)

        Lattice(fig=static_fig, ax=static_ax, plots=False)

        static_canvas = FigureCanvas(static_fig)
        layout.addWidget(static_canvas)
        self.addToolBar(NavigationToolbar(static_canvas, self))
        static_ax.mouse_init()

        options = OptionsWindow()
        layout.addWidget(options)


class OptionsWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        button = QtWidgets.QPushButton("PyQt5 Button", self)
        button.setToolTip('this is an example button')
        button.move(1000, 70)

        button.clicked.connect(self.on_click)

    # @pyqtSlot()
    def on_click(self):
        print('PyQt5 button click')


if __name__ == "__main__":
    Lattice()
    qapp = QtWidgets.QApplication(sys.argv)
    app = ApplicationWindow()
    app.show()
    qapp.exec_()
