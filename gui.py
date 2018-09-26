import sys
from cmp import *
from itertools import compress

from PyQt5 import QtWidgets as QW, QtGui as QG, QtCore as QC

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
            'user_colors': ['e'] * 5,
            'enabled_user_colors': [],
            'preset_colors': ['e'] * 4,
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
            'beta': 70,
            'gamma': 60,
            'lattice': 'simple cubic',
            'max_preset_basis': 4,
            'max_basis': 5
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

        # A list of allowed colors
        self.allowed_colors = ['b', 'r', 'g', 'y', 'k', 'w', 'm', 'c', 'e']

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
            field.setValidator(QG.QDoubleValidator(decimals=3))

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
        # self.create_preset_basis(
        #     n_basis=self.lattice_config['max_preset_basis'])
        # self.layout_options.addWidget(QW.QLabel('Preset specified basis'))
        # self.layout_options.addLayout(self.layout_preset_basis)
        self.create_user_basis()
        # self.layout_options.addWidget(QW.QLabel('Basis'))
        self.layout_options.addLayout(self.layout_basis)

    # def create_preset_basis(self, n_basis):
    #     # So far the largest number of atoms in a preset basis is 4.
    #     self.layout_preset_basis = QW.QGridLayout()
    #     names = ['x', 'y', 'z', 'color']
    #     for n, name in enumerate(names):
    #         label = QW.QLabel(name)
    #         label.setAlignment(QC.Qt.AlignCenter)
    #         self.layout_preset_basis.addWidget(label, 0, n)
    #     n_coords = 3
    #     self.preset_basis_coord_widgets = np.empty((n_basis, n_coords),
    #                                                dtype=object)
    #     self.preset_basis_color_widgets = np.empty(n_basis, dtype=object)
    #     for i in range(n_basis):
    #         for j in range(n_coords):
    #             el = QW.QLineEdit()
    #             el.setEnabled(False)
    #             if i == 0:
    #                 el.setText('0')
    #             self.preset_basis_coord_widgets[i, j] = el
    #             self.layout_preset_basis.addWidget(el, i + 1, j)
    #         el = QW.QLineEdit()
    #         el.setEnabled(True)
    #         el.setMaxLength(1)
    #         el.setText('e')
    #         el.editingFinished.connect(
    #             lambda i=i, el=el: self.update_basis_color(
    #                 'preset', i, el.text()))
    #         self.preset_basis_color_widgets[i] = el
    #         self.layout_preset_basis.addWidget(el, i + 1, n_coords)

    def create_user_basis(self):
        # Basis-stuff
        self.layout_basis = QW.QGridLayout()
        n_basis = self.lattice_config['max_basis']
        n_coords = 3
        names = ['x', 'y', 'z', 'color']
        for n, name in enumerate(names):
            label = QW.QLabel(name)
            label.setAlignment(QC.Qt.AlignCenter)
            self.layout_basis.addWidget(label, 0, n)
        self.basis_coord_widgets = np.empty((n_basis, n_coords),
                                            dtype=object)
        self.basis_color_widgets = []
        self.basis_check_widgets = []
        # Create all the basis coordinate widgets
        for i in range(n_basis):
            for j in range(n_coords):
                # A QLineEdit for each of the basis atoms coordinates
                el = QW.QLineEdit()
                el.setText(str(0))
                el.setEnabled(False)
                el.setValidator(
                    QG.QDoubleValidator(decimals=3))
                # Pass basis and coordinate number, along with value, to
                # update_basis_val
                el.editingFinished.connect(
                    lambda i=i, j=j, el=el:
                        self.update_basis_val(i, j, el.text()))
                # Add the QLineEdit to the array of basis coordinate widgets
                # and to the layout
                self.basis_coord_widgets[i, j] = el
                self.layout_basis.addWidget(el, i + 1, j)

            # Add a color lineedit for each basis atom
            el = QW.QLineEdit()
            el.setText(str('e'))
            el.setMaxLength(1)
            el.editingFinished.connect(
                lambda i=i, el=el:
                    self.update_basis_color('user', i, el.text()))
            self.basis_color_widgets.append(el)
            self.layout_basis.addWidget(el, i + 1, n_coords)

            # Add a checkbox for each basis atom
            check = QW.QCheckBox()
            check.setChecked(False)

            # For some reason this doesn't work when we do it in a loop...
            # check.stateChanged.connect(
            #     lambda i=i: self.hide_basis_widgets(i))

            # Add the checkbox to the list of widgets, and the layout.
            self.basis_check_widgets.append(check)
            self.layout_basis.addWidget(check, i + 1, n_coords + 2)

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

    def update_lattice_name(self, text):
        # Set the lattice name
        self.lattice_config['lattice'] = text

        # disable all fields and enable only those needed. It's the easiest
        # way, if maybe a bit redundant
        for le in self.param_fields:
            le.setEnabled(False)
        for n in self.needed_params[text]:
            self.param_fields[n].setEnabled(True)

        # We should also load the default values for the parameters
        for n, param in enumerate(self.parameter_names):
            self.lattice_config[param] = self.default_config[param]
            self.param_fields[n].setText(str(self.lattice_config[param]))

        # And then we update the lattice
        self.update_lattice()

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
        print(basis)
        n_basis = np.atleast_2d(basis).shape[0]
        # Update primitive lattice vectors and (preset) basis.
        self.lattice_config.update(dict(zip(('a1', 'a2', 'a3', 'user_basis'),
                                            (a1, a2, a3, basis))))
        if n_basis > 1:
            # basis_padding = np.zeros((
            #     self.lattice_config['max_basis'] - n_basis))
            # basis = np.hstack((basis, basis_padding))
            self.start_lattice_with_basis()
        else:
            # If the user chose a preset without a basis, we just give it an
            # empty basis, and disable those that aren't needed
            basis = np.zeros((self.lattice_config['max_basis'], 3))
            self.lattice_config['user_basis'] = basis
            self.start_lattice_without_basis()
        # self.update_preset_basis_widgets()
        self.plot_lattice()

    def update_preset_basis_widgets(self):
        basis = np.atleast_2d(self.lattice_config['user_basis'])
        n_basis = basis.shape[0]
        for el in self.basis_coord_widgets[:n_basis].flatten():
            el.setText("{0:.3f}".format(coord))

        for el in self.basis_coord_widgets[n_basis:].flatten():
            el.setText('')

    def update_config_parameter(self, param, text):
        # This function updates the relevant parameter in the lattice_config
        # dict, but only if the text is a float!
        try:
            self.lattice_config[param] = float(text)
        except ValueError:
            pass
        self.update_lattice()

    def update_basis_color(self, type_, num, text):
        colors = self.lattice_config['{}_colors'.format(type_)]
        text = text.lower()
        if text not in self.allowed_colors:
            print('you did wrong!')
        else:
            colors[num] = text
            self.plot_lattice()

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
        new_colors = self.lattice_config['user_colors']
        new_colors = list(compress(new_colors, enabled_basis_atoms))
        self.lattice_config['enabled_user_colors'] = new_colors
        self.plot_lattice()

    def start_lattice_with_basis(self):
        # This function runs when the user chooses a preset which includes a
        # basis (like wurtzite). It hides the unneeded basis widgets (like all
        # the checkboxes), and disables the ones that are needed, so the user
        # can't change them.
        # n_basis = np.atleast_2d(self.lattice_config['user_basis']).shape[0]
        # max_basis = self.lattice_config['max_basis']
        # unneeded_coords = self.basis_coord_widgets[n_basis:max_basis]
        # for el in unneeded_coords.flatten():
        #     el.setHidden(True)
        for i in self.basis_check_widgets:
            i.setHidden(True)
        # for i in self.basis_color_widgets[n_basis:max_basis]:
        #     i.setHidden(True)

    def start_lattice_without_basis(self):
        # This function runs when the user chooses a preset which does not
        # include a basis (ie a pure Bravais lattice). It shows all basis
        # widgets again and loads an empty basis
        n_basis = self.lattice_config['max_basis']
        basis = np.zeros((n_basis, 3))
        self.lattice_config['user_basis'] = basis
        # for el in self.basis_coord_widgets.flatten():
        #     el.setHidden(False)
        for el in self.basis_check_widgets:
            el.setHidden(False)
            el.setEnabled(False)
        # for el in self.basis_color_widgets:
        #     el.setHidden(False)

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

        user_colors = self.lattice_config['enabled_user_colors']
        preset_colors = self.lattice_config['preset_colors']
        n_preset_basis = np.atleast_2d(preset_basis).shape[0]
        preset_colors = preset_colors[:n_preset_basis]
        colors = preset_colors + user_colors
        # List comprehension to change 'e' to 'xkcd:cement'
        colors = ['xkcd:cement' if i == 'e' else i for i in colors]

        # Plot the new lattice
        self.static_fig, self.static_ax = Lattice(a1=a1, a2=a2, a3=a3,
                                                  basis=basis,
                                                  colors=colors,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False)

        # Remember to have the canvas draw it!
        self.static_canvas.draw()

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
