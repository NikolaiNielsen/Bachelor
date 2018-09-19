import sys
from cmp import *

from PyQt5 import QtWidgets as QW, QtGui as QG

from matplotlib.backends.backend_qt5agg import (
    FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class full_window(QW.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QW.QWidget()
        self.setCentralWidget(self._main)
        main_layout = QW.QHBoxLayout(self._main)
        # A shortcut to close the app.
        self.closer = QW.QShortcut(QG.QKeySequence('Ctrl+Q'), self, self.quit)
        self.title = "Crystal Structure"
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
                               'b': 1.2,
                               'c': 1.5,
                               'theta': 80,
                               'beta': 90,
                               'gamma': 90,
                               'lattice': 'simple cubic'}
        self.basis = np.zeros((5, 3))
        self.create_options()
        main_layout.addLayout(self.layout_options)

        self.params = [self.le_a, self.le_b, self.le_c,
                       self.le_theta, self.le_beta, self.le_gamma]
        self.needed_params = {
            'simple cubic': [self.le_a],
            'bcc': [self.le_a],
            'fcc': [self.le_a],
            'base centred cubic': [self.le_a],
            'tetragonal': [self.le_a, self.le_b],
            'tetragonal body centred': [self.le_a, self.le_b],
            'tetragonal face centred': [self.le_a, self.le_b],
            'orthorhombic': [self.le_a, self.le_b, self.le_c],
            'orthorhombic body centred': [self.le_a, self.le_b, self.le_c],
            'orthorhombic face centred': [self.le_a, self.le_b, self.le_c],
            'orthorhombic base centred': [self.le_a, self.le_b, self.le_c],
            'simple monoclinic': [self.le_a, self.le_b, self.le_c,
                                  self.le_theta],
            'base centred monoclinic': [self.le_a, self.le_b, self.le_c,
                                        self.le_theta],
            'hexagonal': [self.le_a],
            'triclinic': [self.le_a, self.le_b, self.le_c,
                          self.le_theta, self.le_beta, self.le_gamma],
            'rhombohedral': [self.le_a],
            'diamond': [self.le_a],
            'wurtzite': [self.le_a],
            'zincblende': [self.le_a, self.le_b]
        }
        for le in self.params:
            le.setEnabled(False)
        for le in self.needed_params[self.lattice_config['lattice']]:
            le.setEnabled(True)
        self.static_fig, self.static_ax = Lattice(returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        main_layout.addWidget(self.static_canvas)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def create_options(self):
        self.layout_options = QW.QVBoxLayout()
        # Create the "show plot" button
        self.button_show = QW.QPushButton("Update plot", self)
        self.button_show.setToolTip('this is an example button')
        self.button_show.clicked.connect(self.update_lattice)

        # Create the lattice chooser dropdown
        self.lattice_chooser = QW.QComboBox(self)
        self.lattice_chooser.addItems(self.lattices)
        self.lattice_chooser.activated[str].connect(self.update_lattice_name)

        # Create parameter fields
        self.text_a = QW.QLabel('a', self)
        self.le_a = QW.QLineEdit()
        self.le_a.setValidator(QG.QDoubleValidator(decimals=2))
        # self.le_a.setMaxLength(2)
        self.le_a.setText(str(self.lattice_config['a']))
        self.le_a.returnPressed.connect(
            lambda: self.update_config_parameter('a', self.le_a.text()))

        self.text_b = QW.QLabel('b', self)
        self.le_b = QW.QLineEdit()
        self.le_b.setValidator(QG.QDoubleValidator(decimals=2))
        # self.le_b.setMaxLength(2)
        self.le_b.setText(str(self.lattice_config['b']))
        self.le_b.returnPressed.connect(
            lambda: self.update_config_parameter('b', self.le_b.text()))

        self.text_c = QW.QLabel('c', self)
        self.le_c = QW.QLineEdit()
        self.le_c.setValidator(QG.QDoubleValidator(decimals=2))
        # self.le_c.setMaxLength(2)
        self.le_c.setText(str(self.lattice_config['c']))
        self.le_c.returnPressed.connect(
            lambda: self.update_config_parameter('c', self.le_c.text()))

        self.text_theta = QW.QLabel('theta, degrees', self)
        self.le_theta = QW.QLineEdit()
        self.le_theta.setValidator(QG.QDoubleValidator(decimals=2))
        # self.le_theta.setMaxLength(3)
        self.le_theta.setText(str(self.lattice_config['theta']))
        self.le_theta.returnPressed.connect(
            lambda: self.update_config_parameter(
                'theta', self.le_theta.text()))

        self.text_beta = QW.QLabel('beta, degrees', self)
        self.le_beta = QW.QLineEdit()
        self.le_beta.setValidator(QG.QDoubleValidator(decimals=2))
        # self.le_beta.setMaxLength(3)
        self.le_beta.setText(str(self.lattice_config['beta']))
        self.le_beta.returnPressed.connect(
            lambda: self.update_config_parameter(
                'beta', self.le_beta.text()))

        self.text_gamma = QW.QLabel('gamma, degrees', self)
        self.le_gamma = QW.QLineEdit()
        self.le_gamma.setValidator(QG.QDoubleValidator(decimals=2))
        # self.le_gamma.setMaxLength(3)
        self.le_gamma.setText(str(self.lattice_config['gamma']))
        self.le_gamma.returnPressed.connect(
            lambda: self.update_config_parameter(
                'gamma', self.le_gamma.text()))

        # Setup stuff in a layout
        self.layout_parameters = QW.QFormLayout()
        self.layout_parameters.addRow(self.button_show, self.lattice_chooser)
        self.layout_parameters.addRow(self.text_a, self.le_a)
        self.layout_parameters.addRow(self.text_b, self.le_b)
        self.layout_parameters.addRow(self.text_c, self.le_c)
        self.layout_parameters.addRow(self.text_theta, self.le_theta)
        self.layout_parameters.addRow(self.text_beta, self.le_beta)
        self.layout_parameters.addRow(self.text_gamma, self.le_gamma)

        self.layout_options.addLayout(self.layout_parameters)
        self.create_basis()
        self.layout_options.addLayout(self.layout_basis)

    def create_basis(self):
        # Basis-stuff
        self.basis1_x = QW.QLineEdit()
        self.basis1_x.editingFinished.connect(
            lambda: self.update_basis_val(0, 0, self.basis1_x.text()))

        self.basis1_y = QW.QLineEdit()
        self.basis1_y.editingFinished.connect(
            lambda: self.update_basis_val(0, 1, self.basis1_y.text()))

        self.basis1_z = QW.QLineEdit()
        self.basis1_z.editingFinished.connect(
            lambda: self.update_basis_val(0, 2, self.basis1_z.text()))

        self.basis1_enable = QW.QCheckBox()
        self.basis1_enable.setChecked(False)
        self.basis1_enable.stateChanged.connect(
            lambda: self.hide_basis_widgets(0))

        # Basis 2
        self.basis2_x = QW.QLineEdit()
        self.basis2_x.editingFinished.connect(
            lambda: self.update_basis_val(1, 0, self.basis2_x.text()))

        self.basis2_y = QW.QLineEdit()
        self.basis2_y.editingFinished.connect(
            lambda: self.update_basis_val(1, 1, self.basis2_y.text()))

        self.basis2_z = QW.QLineEdit()
        self.basis2_z.editingFinished.connect(
            lambda: self.update_basis_val(1, 2, self.basis2_z.text()))

        self.basis2_enable = QW.QCheckBox()
        self.basis2_enable.setChecked(False)
        self.basis2_enable.stateChanged.connect(
            lambda: self.hide_basis_widgets(1))

        # Basis 3
        self.basis3_x = QW.QLineEdit()
        self.basis3_x.editingFinished.connect(
            lambda: self.update_basis_val(2, 0, self.basis3_x.text()))

        self.basis3_y = QW.QLineEdit()
        self.basis3_y.editingFinished.connect(
            lambda: self.update_basis_val(2, 1, self.basis3_y.text()))

        self.basis3_z = QW.QLineEdit()
        self.basis3_z.editingFinished.connect(
            lambda: self.update_basis_val(2, 2, self.basis3_z.text()))

        self.basis3_enable = QW.QCheckBox()
        self.basis3_enable.setChecked(False)
        self.basis3_enable.stateChanged.connect(
            lambda: self.hide_basis_widgets(2))

        # Basis 4
        self.basis4_x = QW.QLineEdit()
        self.basis4_x.editingFinished.connect(
            lambda: self.update_basis_val(3, 0, self.basis4_x.text()))

        self.basis4_y = QW.QLineEdit()
        self.basis4_y.editingFinished.connect(
            lambda: self.update_basis_val(3, 1, self.basis4_y.text()))

        self.basis4_z = QW.QLineEdit()
        self.basis4_z.editingFinished.connect(
            lambda: self.update_basis_val(3, 2, self.basis4_z.text()))

        self.basis4_enable = QW.QCheckBox()
        self.basis4_enable.setChecked(False)
        self.basis4_enable.stateChanged.connect(
            lambda: self.hide_basis_widgets(3))

        # Basis 5
        self.basis5_x = QW.QLineEdit()
        self.basis5_x.editingFinished.connect(
            lambda: self.update_basis_val(4, 0, self.basis5_x.text()))

        self.basis5_y = QW.QLineEdit()
        self.basis5_y.editingFinished.connect(
            lambda: self.update_basis_val(4, 1, self.basis5_y.text()))

        self.basis5_z = QW.QLineEdit()
        self.basis5_z.editingFinished.connect(
            lambda: self.update_basis_val(4, 2, self.basis5_z.text()))

        self.basis5_enable = QW.QCheckBox()
        self.basis5_enable.setChecked(False)
        self.basis5_enable.stateChanged.connect(
            lambda: self.hide_basis_widgets(4))

        self.basis_widgets = [[
            self.basis1_x,
            self.basis1_y,
            self.basis1_z,
            self.basis1_enable],
            [
            self.basis2_x,
            self.basis2_y,
            self.basis2_z,
            self.basis2_enable],
            [
            self.basis3_x,
            self.basis3_y,
            self.basis3_z,
            self.basis3_enable],
            [
            self.basis4_x,
            self.basis4_y,
            self.basis4_z,
            self.basis4_enable],
            [
            self.basis5_x,
            self.basis5_y,
            self.basis5_z,
            self.basis5_enable]]
        self.layout_basis = QW.QGridLayout()

        for i in range(5):
            for j in range(4):
                if j != 3:
                    self.basis_widgets[i][j].setText('0')
                    self.basis_widgets[i][j].setEnabled(False)
                    self.basis_widgets[i][j].setValidator(
                        QG.QDoubleValidator(decimals=2))
                self.layout_basis.addWidget(self.basis_widgets[i][j], i, j)

    def update_lattice(self):
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
        self.plot_lattice(a1=a1, a2=a2, a3=a3, basis=basis)

    def update_lattice_name(self, text):
        self.lattice_config['lattice'] = text
        for le in self.params:
            le.setEnabled(False)
        for le in self.needed_params[text]:
            le.setEnabled(True)
        self.update_lattice()

    def update_config_parameter(self, param, text):
        # This function updates the relevant parameter in the lattice_config
        # dict. But only if the text is a float!
        try:
            self.lattice_config[param] = float(text)
        except ValueError:
            pass
        self.update_lattice()

    def plot_lattice(self, a1, a2, a3, basis):
        # This function takes the values from lattice_config and uses them to
        # update the plot
        self.static_ax.clear()
        self.static_fig, self.static_ax = Lattice(a1=a1, a2=a2, a3=a3,
                                                  basis=basis,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False)
        self.static_canvas.activateWindow()

    def update_basis_val(self, basis_no, coord_no, val):
        self.basis[basis_no, coord_no] = float(val)
        self.update_basis()

    def hide_basis_widgets(self, basis_no):
        checkbox = self.basis_widgets[basis_no][3]
        for le in self.basis_widgets[basis_no][:3]:
            le.setEnabled(checkbox.isChecked())
        self.update_basis()

    def update_basis(self):
        enabled_basis_atoms = []
        for i in self.basis_widgets:
            enabled_basis_atoms.append(i[3].isChecked())
        new_basis = np.vstack(([0, 0, 0], self.basis[enabled_basis_atoms]))
        self.lattice_config['basis'] = new_basis
        self.update_lattice()

    def quit(self):
        sys.exit()


def main():
    qapp = QW.QApplication(sys.argv)
    app = full_window()
    app.show()
    sys.exit(qapp.exec_())


def not_main():
    qapp = QW.QApplication(sys.argv)
    app = full_window()
    app.show()
    return qapp, app


if __name__ == "__main__":
    main()
