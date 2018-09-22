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
        self.setWindowTitle(self.title)

        # A list of names for available lattice presets
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

        # A dictionary of the default values for lattice plotting
        self.default_config = {
            'a1': d[0],
            'a2': d[1],
            'a3': d[2],
            'preset_basis': d[3],
            'colors': d[4],
            'sizes': d[5],
            'enabled_user_basis': np.empty((1, 3)),
            'user_basis': np.zeros((5, 3)),
            'lim_type': d[6],
            'grid_type': None,
            'max_': d[8],
            'a': 1,
            'b': 1.2,
            'c': 1.5,
            'theta': 80,
            'beta': 90,
            'gamma': 90,
            'lattice': 'simple cubic'
        }
        # Needed parameters for each lattice (a, b, c, theta, beta, gamma)
        self.needed_params = {
            'simple cubic': [0],
            'bcc': [0],
            'fcc': [0],
            'base centred cubic': [0],
            'tetragonal': [0, 1],
            'tetragonal body centred': [0, 1],
            'tetragonal face centred': [0, 1],
            'orthorhombic': [0, 1, 2],
            'orthorhombic body centred': [0, 1, 2],
            'orthorhombic face centred': [0, 1, 2],
            'orthorhombic base centred': [0, 1, 2],
            'simple monoclinic': [0, 1, 2, 3],
            'base centred monoclinic': [0, 1, 2, 3],
            'hexagonal': [0],
            'triclinic': [0, 1, 2, 3, 4, 5],
            'rhombohedral': [0],
            'diamond': [0],
            'wurtzite': [0, 1],
            'zincblende': [0]
        }
        # Copy of the default config. This is what the user'll actually change
        self.lattice_config = self.default_config.copy()

        # We create the options and add it to our main layout (it also creates
        # the basis fiels)
        self.create_options()
        main_layout.addLayout(self.layout_options)

        # Enable only the needed parameter fields.
        for n in self.needed_params[self.lattice_config['lattice']]:
            self.param_fields[n].setEnabled(True)

        # Create the default plot and return the figure and axis objects for
        # it. Then create the FigureCanvas, add them all to the layout and add
        # a toolbar. Lastly enable mouse support for Axes3D
        self.static_fig, self.static_ax = Lattice(returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        main_layout.addWidget(self.static_canvas)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def create_options(self):
        self.parameter_names = ['a', 'b', 'c', 'theta', 'beta', 'gamma']
        self.param_labels = []
        self.param_fields = []
        self.layout_options = QW.QVBoxLayout()
        # Create the "show plot" button
        self.button_show = QW.QPushButton("Update plot", self)
        self.button_show.clicked.connect(self.update_lattice)

        # Create the lattice chooser dropdown
        self.lattice_chooser = QW.QComboBox(self)
        self.lattice_chooser.addItems(self.lattices)
        self.lattice_chooser.activated[str].connect(self.update_lattice_name)

        # Create the parameter layout
        self.layout_parameters = QW.QFormLayout()
        self.layout_parameters.addRow(self.button_show, self.lattice_chooser)

        for name in self.parameter_names:
            # Create all the parameter labels and fields.

            self.param_labels.append(QW.QLabel(name, self))
            field = QW.QLineEdit()

            # Only allow floats and 2 decimals to input
            field.setValidator(QG.QDoubleValidator(decimals=2))

            # Populate with default values
            field.setText(str(self.lattice_config[name]))
            field.setEnabled(False)

            # Pass both parameter name and value to update_config_parameter
            field.returnPressed.connect(
                lambda name=name, el=field:
                    self.update_config_parameter(name, el.text()))

            # Add the parameter field to the list
            self.param_fields.append(field)

        # When everything has been created we add the parameters to the layout,
        # along with the labels
        for n in range(len(self.param_labels)):
            self.layout_parameters.addRow(self.param_labels[n],
                                          self.param_fields[n])

        self.layout_options.addLayout(self.layout_parameters)
        self.create_basis()
        self.layout_options.addLayout(self.layout_basis)

    def create_basis(self):
        # Basis-stuff
        self.layout_basis = QW.QGridLayout()
        no_basis = 5
        no_coords = 3
        self.basis_coord_widgets = np.empty((no_basis, no_coords),
                                            dtype=object)
        self.basis_check_widgets = []
        # Create all the basis coordinate widgets
        for i in range(no_basis):
            for j in range(no_coords):
                # A QLineEdit for each of the basis atoms coordinates
                el = QW.QLineEdit()
                el.setText(str(0))
                el.setEnabled(False)
                el.setValidator(
                    QG.QDoubleValidator(decimals=2))
                # Pass basis and coordinate number, along with value, to
                # update_basis_val
                el.editingFinished.connect(
                    lambda i=i, j=j, el=el:
                        self.update_basis_val(i, j, el.text()))
                # Add the QLineEdit to the array of basis coordinate widgets
                # and to the layout
                self.basis_coord_widgets[i, j] = el
                self.layout_basis.addWidget(el, i, j)

            # Add a checkbox for each basis atom
            check = QW.QCheckBox()
            check.setChecked(False)

            # For some reason this doesn't work when we do it in a loop...
            # check.stateChanged.connect(
            #     lambda i=i: self.hide_basis_widgets(i))

            # Add the checkbox to the list of widgets, and the layout.
            self.basis_check_widgets.append(check)
            self.layout_basis.addWidget(check, i, no_coords + 1)

        # It's ugly but it works. We make the checkbox to stuff
        self.basis_check_widgets[0].stateChanged.connect(
            lambda: self.hide_basis_widgets(0))
        self.basis_check_widgets[1].stateChanged.connect(
            lambda: self.hide_basis_widgets(1))
        self.basis_check_widgets[2].stateChanged.connect(
            lambda: self.hide_basis_widgets(2))
        self.basis_check_widgets[3].stateChanged.connect(
            lambda: self.hide_basis_widgets(3))
        self.basis_check_widgets[4].stateChanged.connect(
            lambda: self.hide_basis_widgets(4))

    def update_lattice(self):
        # Grab a new lattice based on the parameters in lattice_config
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

        # Update primitive lattice vectors and (preset) basis.
        self.lattice_config.update(dict(zip(('a1', 'a2', 'a3', 'preset_basis'),
                                            (a1, a2, a3, basis))))
        self.plot_lattice()

    def update_lattice_name(self, text):
        # Set the lattice name
        self.lattice_config['lattice'] = text

        # disable all fields and enable only those needed. It's the easiest
        # way, if maybe a bit redundant
        for le in self.param_fields:
            le.setEnabled(False)
        for n in self.needed_params[text]:
            self.param_fields[n].setEnabled(True)

        self.update_lattice()

    def update_config_parameter(self, param, text):
        # This function updates the relevant parameter in the lattice_config
        # dict, but only if the text is a float!
        try:
            self.lattice_config[param] = float(text)
        except ValueError:
            pass
        self.update_lattice()

    def plot_lattice(self):
        # This function takes the values from lattice_config and uses them to
        # update the plot.

        # Clear the axes
        self.static_ax.clear()

        # Grab lattice vectors and basis(es) from lattice_config
        a1 = self.lattice_config['a1']
        a2 = self.lattice_config['a2']
        a3 = self.lattice_config['a3']
        preset_basis = self.lattice_config['preset_basis']
        user_basis = self.lattice_config['enabled_user_basis']
        basis = np.vstack((preset_basis, user_basis))

        # Plot the new lattice
        self.static_fig, self.static_ax = Lattice(a1=a1, a2=a2, a3=a3,
                                                  basis=basis,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False)

        # Remember to have the canvas draw it!
        self.static_canvas.draw()

    def update_basis_val(self, basis_no, coord_no, val):
        self.lattice_config['user_basis'][basis_no, coord_no] = float(val)
        self.update_basis()

    def hide_basis_widgets(self, basis_no):
        # enable or disable basis coord widgets and update the basis
        checkbox = self.basis_check_widgets[basis_no]
        for le in self.basis_coord_widgets[basis_no]:
            le.setEnabled(checkbox.isChecked())
        self.update_basis()

    def update_basis(self):
        # We get a list of basis atoms that are enabled
        enabled_basis_atoms = []
        for i in self.basis_check_widgets:
            enabled_basis_atoms.append(i.isChecked())

        # update the enabled_user_basis config and plot the lattice with the
        # new basis
        new_basis = self.lattice_config['user_basis'][enabled_basis_atoms]
        self.lattice_config['enabled_user_basis'] = new_basis
        self.plot_lattice()

    def tester(self, i):
        print(i)

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
