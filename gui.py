import sys
from cmp import *

from PyQt5 import QtWidgets

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class plot_window(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QtWidgets.QWidget()
        self.setCentralWidget(self._main)
        layout = QtWidgets.QHBoxLayout(self._main)

        self.static_fig, self.static_ax = Lattice(returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        layout.addWidget(self.static_canvas)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def update_lattice(self, lattice_name):
        self.static_ax.clear()
        self.static_fig, self.static_ax = Lattice(lattice_name=lattice_name,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False)


class options_window(QtWidgets.QWidget):
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
        self.default_config = d
        self.lattice_config = {'a1': d[0],
                               'a2': d[1],
                               'a3': d[2],
                               'basis': d[3],
                               'colors': d[4],
                               'sizes': d[5],
                               'lim_type': d[6],
                               'grid_type': None,
                               'max_': d[8]}
        self.init_ui()

    def init_ui(self):

        self.button_show = QtWidgets.QPushButton("Show plot", self)
        self.button_show.setToolTip('this is an example button')
        self.button_show.move(0, 0)
        self.button_show.clicked.connect(self.show_plot)

        self.lattice_chooser = QtWidgets.QComboBox(self)
        self.lattice_chooser.addItems(self.lattices)
        self.lattice_chooser.move(0, 100)
        self.lattice_chooser.activated[str].connect(self.change_lattice)

        self.plot = plot_window()
        self.plot.show()

    def show_plot(self):
        self.plot.show()

    def change_lattice(self, text):
        self.plot.update_lattice(text)


def main():
    qapp = QtWidgets.QApplication(sys.argv)
    app = options_window()
    app.show()
    sys.exit(qapp.exec_())


if __name__ == "__main__":
    main()
