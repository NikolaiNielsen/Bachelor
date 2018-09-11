import sys
# from matplotlib.figure import Figure
from cmp import *

from PyQt5 import QtWidgets
from PyQt5.QtCore import pyqtSlot

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class ApplicationWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        # layout = QtWidgets.QHBoxLayout(self._main)

        static_fig, static_ax = Lattice(returns=True, plots=False)
        static_canvas = FigureCanvas(static_fig)
        # static_fic = plt.figure(figsize=(5, 3))
        # static_canvas = FigureCanvas(static_fic)
        # layout.addWidget(static_canvas)
        self.addToolBar(NavigationToolbar(static_canvas, self))


        # self._static_ax = static_canvas.figure.subplots()
        # t = np.linspace(0, 10, 501)
        # self._static_ax.plot(t, np.tan(t), ".")
        # static_ax.mouse_init()

        # options = OptionsWindow()
        # layout.addWidget(options)
        # self.show()


class OptionsWindow(QtWidgets.QWidget):
    def __init__(self):
        super().__init__()
        self.title = "Options window"
        self.left = 10
        self.top = 10
        self.width = 320
        self.height = 200
        self.initUI()

    def initUI(self):
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)

        button = QtWidgets.QPushButton("PyQt5 Button", self)
        button.setToolTip('this is an example button')
        button.move(100, 70)
        button.clicked.connect(self.on_click)
        self.dialog = ApplicationWindow()

    def on_click(self):
        self.dialog.show()


def main():
    qapp = QtWidgets.QApplication(sys.argv)
    app = OptionsWindow()
    app.show()
    sys.exit(qapp.exec_())


if __name__ == "__main__":
    main()
