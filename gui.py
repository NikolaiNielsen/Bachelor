import sys
from cmp import *

from PyQt5 import QtWidgets as QW, QtGui as QG

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class plot_window(QW.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QW.QWidget()
        self.setCentralWidget(self._main)
        layout = QW.QHBoxLayout(self._main)

        self.static_fig, self.static_ax = Lattice(returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        layout.addWidget(self.static_canvas)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def update_lattice(self, a1, a2, a3, basis):
        self.static_ax.clear()
        self.static_fig, self.static_ax = Lattice(a1=a1, a2=a2, a3=a3,
                                                  basis=basis,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False)


class full_window(QW.QWidget):
    def __init__(self):
        super().__init__()
        self._main = QT.QWidget()
        self.setCentralWidget(self._main)
        main_layout = QW.QHBoxLayout(self._main)

        self.title = "Options window"
        self.left = 10
        self.top = 10
        self.width = 320
        self.height = 200
        self.lattices = ['simple cubic', 'bcc', 'fcc', 'base centred cubic',
                         'tetragonal', 'tetragonal body centred',
                         'tetragonal face centred', 'orthorhombic',
                         'orthorhombic body centred',
                         'orthorhombic face centred',
                         'orthorhombic base centred',
                         'simple monoclinic',
                         'base centred monoclinic',
                         'hexagonal',
                         'triclinic',
                         'rhombohedral',
                         'diamond',
                         'wurtzite',
                         'zincblende']
        self.setWindowTitle(self.title)
        # self.setGeometry(self.left, self.top, self.width, self.height)
        self.default_config = d
        self.lattice_config = {'a1': d[0],
                               'a2': d[1],
                               'a3': d[2],
                               'basis': d[3],
                               'colors': d[4],
                               'sizes': d[5],
                               'lim_type': d[6],
                               'grid_type': None,
                               'max_': d[8],
                               'a': 1,
                               'b': 1,
                               'c': 1,
                               'theta': 80,
                               'beta': 90,
                               'gamma': 90,
                               'lattice': 'simple cubic'}
        self.padding = 0
        self.create_options()
        main_layout.addWidget(self.layout_box)

        self.static_fig, self.static_ax = Lattice(returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        main_layout.addWidget(self.static_canvas)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def create_options(self):

        # Create the "show plot" button
        self.button_show = QW.QPushButton("Show plot", self)
        self.button_show.setToolTip('this is an example button')
        self.button_show.clicked.connect(self.update_plot)

        # Create the lattice chooser dropdown
        self.lattice_chooser = QW.QComboBox(self)
        self.lattice_chooser.addItems(self.lattices)
        self.lattice_chooser.activated[str].connect(self.update_lattice_name)

        # Create parameter fields
        self.text_a = QW.QLabel('a', self)
        self.le_a = QW.QLineEdit()
        self.le_a.setValidator(QG.QIntValidator())
        self.le_a.setMaxLength(2)
        self.le_a.textChanged.connect(
            lambda: self.update_config('a', self.le_a.text()))

        self.text_b = QW.QLabel('b', self)
        self.le_b = QW.QLineEdit()
        self.le_b.setValidator(QG.QIntValidator())
        self.le_b.setMaxLength(2)
        self.le_b.textChanged.connect(
            lambda: self.update_config('b', self.le_b.text()))

        self.text_c = QW.QLabel('c', self)
        self.le_c = QW.QLineEdit()
        self.le_c.setValidator(QG.QIntValidator())
        self.le_c.setMaxLength(2)
        self.le_c.textChanged.connect(
            lambda: self.update_config('c', self.le_c.text()))

        self.text_theta = QW.QLabel('theta, degrees', self)
        self.le_theta = QW.QLineEdit()
        self.le_theta.setValidator(QG.QIntValidator())
        self.le_theta.setMaxLength(3)
        self.le_theta.textChanged.connect(
            lambda: self.update_config('theta', self.le_theta.text()))

        self.text_beta = QW.QLabel('beta, degrees', self)
        self.le_beta = QW.QLineEdit()
        self.le_beta.setValidator(QG.QIntValidator())
        self.le_beta.setMaxLength(3)
        self.le_beta.textChanged.connect(
            lambda: self.update_config('beta', self.le_beta.text()))

        self.text_gamma = QW.QLabel('gamma, degrees', self)
        self.le_gamma = QW.QLineEdit()
        self.le_gamma.setValidator(QG.QIntValidator())
        self.le_gamma.setMaxLength(3)
        self.le_gamma.textChanged.connect(
            lambda: self.update_config('gamma', self.le_gamma.text()))

        # Setup stuff in a layout
        self.layout_box = QW.QFormLayout()
        self.layout_box.addRow(self.button_show, self.lattice_chooser)
        self.layout_box.addRow(self.text_a, self.le_a)
        self.layout_box.addRow(self.text_b, self.le_b)
        self.layout_box.addRow(self.text_c, self.le_c)
        self.layout_box.addRow(self.text_theta, self.le_theta)
        self.layout_box.addRow(self.text_beta, self.le_beta)
        self.layout_box.addRow(self.text_gamma, self.le_gamma)
        # self.setLayout(self.layout_box)

    def update_plot(self):
        a = self.lattice_config['a']
        b = self.lattice_config['b']
        c = self.lattice_config['c']
        theta = self.lattice_config['theta'] * np.pi / 180
        beta = self.lattice_config['beta'] * np.pi / 180
        gamma = self.lattice_config['gamma'] * np.pi / 180
        name = self.lattice_config['lattice']
        (a1, a2, a3), basis, _ = lattices.chooser(lattice_name=name,
                                                  a=a, b=b, c=c,
                                                  theta=theta,
                                                  beta=beta,
                                                  gamma=gamma)
        self.lattice_config.update(dict(zip(('a1', 'a2', 'a3', 'basis'),
                                            (a1, a2, a3, basis))))
        self.update_lattice(a1=a1, a2=a2, a3=a3, basis=basis)

    def update_lattice_name(self, text):
        self.lattice_config['lattice'] = text

    def update_config(self, param, text):
        try:
            self.lattice_config[param] = float(text)
        except ValueError:
            pass

    def update_lattice(self, a1, a2, a3, basis):
        self.static_ax.clear()
        self.static_fig, self.static_ax = Lattice(a1=a1, a2=a2, a3=a3,
                                                  basis=basis,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False)


def main():
    qapp = QW.QApplication(sys.argv)
    app = full_window()
    app.show()
    sys.exit(qapp.exec_())


def not_main():
    qapp = QW.QApplication(sys.argv)
    app = options_window()
    app.show()
    return qapp, app


if __name__ == "__main__":
    main()
