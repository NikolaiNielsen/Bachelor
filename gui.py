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
        self.layout_main = QW.QHBoxLayout(self._main)
        self.lattice_window = lattice_window()
        self.layout_main.addWidget(self.lattice_window)

        # A shortcut to close the app.
        self.closer = QW.QShortcut(QG.QKeySequence('Ctrl+Q'), self, self.quit)
        self.title = "Crystal Structure"
        self.setWindowTitle(self.title)
        bar = self.menuBar()
        programs = bar.addMenu('Programs')
        lattices = programs.addAction('Crystal Structure')
        lattices.triggered.connect(self.create_lattice)
        scattering = programs.addAction('Scattering')
        scattering.triggered.connect(self.create_scattering)
        about = bar.addAction('About')
        about.triggered.connect(self.about_section)

    def create_lattice(self):
        # First we delete the current main layout
        self.delete_layout(self.layout_main)

        # Next we create the new one
        self.lattice_window = lattice_window()
        self.layout_main.addWidget(self.lattice_window)

    def create_scattering(self):
        self.delete_layout(self.layout_main)
        self.scattering_window = scattering_window()
        self.layout_main.addWidget(self.scattering_window)

    def about_section(self):
        author = 'Created by Nikolai Plambech Nielsen, lpk331@alumni.ku.dk.\n'
        author2 = 'The whole code can be found on '
        author3 = 'https://github.com/NikolaiNielsen/bachelor'
        author = author + author2 + author3
        msg = QW.QMessageBox()
        msg.setText('Visualizations of concepts in Condensed Matter Physics')
        msg.setWindowTitle('About')
        msg.setInformativeText(author)
        msg.exec_()

    def delete_layout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                else:
                    self.delete_layout(item.layout())

    def tester(self, arg1, arg2):
        print(arg1, arg2)

    def quit(self):
        sys.exit()


class lattice_window(QW.QMainWindow):
    def __init__(self):
        super().__init__()
        self._main = QW.QWidget()
        self.setCentralWidget(self._main)
        self.layout_main = QW.QHBoxLayout(self._main)
        self.create_variables()
        # We create the options and add it to our main layout (it also creates
        # the basis fiels)
        self.create_options()
        self.layout_main.addLayout(self.layout_options)
        self.create_plot_window()
        self.layout_main.addWidget(self.static_canvas)

    def create_plot_window(self):
        # Create the default plot and return the figure and axis objects for
        # it. Then create the FigureCanvas, add them all to the layout and add
        # a toolbar. Lastly enable mouse support for Axes3D
        self.static_fig, self.static_ax = Lattice(returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def create_variables(self):
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
                         'zincblende',
                         'conventional fcc',
                         'conventional bcc']
        self.colors = [
            'xkcd:cement', 'red', 'blue', 'green', 'cyan',
            'magenta', 'black', 'yellow']
        # A dictionary of the default values for lattice plotting
        self.default_config = {
            'a1': d[0],
            'a2': d[1],
            'a3': d[2],
            'preset_basis': d[3],
            'user_colors': ['xkcd:cement'] * 5,
            'enabled_user_colors': ['xkcd:cement'],
            'preset_colors': ['xkcd:cement'] * 4,
            'sizes': d[5],
            'enabled_user_basis': np.zeros((1, 3)),
            'user_basis': np.zeros((5, 3)),
            'lim_type': d[6],
            'grid_type': None,
            'max_': d[8],
            'a': 1,
            'b': 1.2,
            'c': 1.5,
            'alpha': 80,
            'beta': 70,
            'gamma': 60,
            'lattice': 'simple cubic',
            'max_preset_basis': 4
        }
        # Needed parameters for each lattice (a, b, c, alpha, beta, gamma)
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
            'zincblende': [0],
            'conventional fcc': [0],
            'conventional bcc': [0]
        }
        self.presets_with_basis = {
            'wurtzite': 4,
            'diamond': 2,
            'zincblende': 2,
            'conventional fcc': 4,
            'conventional bcc': 2
        }
        # Copy of the default config. This is what the user'll actually change
        self.lattice_config = self.default_config.copy()

    def create_options(self):
        self.parameter_names = ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        self.parameter_text = ['a', 'b', 'c', 'alpha (degrees)',
                               'beta (degrees)', 'gamma (degrees)']
        self.parameter_tooltips = ['', '', '', 'Angle between side 1 and 3',
                                   'Angle between side 1 and 2',
                                   'Angle between side 2 and 3']
        self.param_labels = []
        self.param_fields = []
        self.layout_options = QW.QVBoxLayout()
        # Create the "show plot" button
        self.button_show = QW.QPushButton("Update plot", self)
        self.button_show.clicked.connect(self.update_lattice)

        # Create the lattice chooser dropdown
        self.lattice_chooser = QW.QComboBox(self)
        self.lattice_chooser.addItems(self.lattices)
        sep_index = len(self.lattices) - len(self.presets_with_basis)
        self.lattice_chooser.insertSeparator(sep_index)
        self.lattice_chooser.activated[str].connect(self.update_lattice_name)

        # Create the parameter layout
        self.layout_parameters = QW.QFormLayout()
        self.layout_parameters.addRow(self.button_show, self.lattice_chooser)

        for n, name in enumerate(self.parameter_names):
            # Create all the parameter labels and fields.
            label = QW.QLabel(self.parameter_text[n], self)
            label.setToolTip(self.parameter_tooltips[n])
            self.param_labels.append(label)
            field = QW.QLineEdit()
            field.setToolTip(self.parameter_tooltips[n])

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
        # Enable only the needed parameter fields.
        for n in self.needed_params[self.lattice_config['lattice']]:
            self.param_fields[n].setEnabled(True)
        self.create_user_basis()

    def create_preset_basis(self, n_basis):
        # So far the largest number of atoms in a preset basis is 4.
        self.layout_preset_basis = QW.QVBoxLayout()
        self.basis_title = QW.QLabel('Basis coordinates')
        self.basis_title.setAlignment(QC.Qt.AlignCenter)

        font = QG.QFont()
        font.setBold(True)
        self.basis_title.setFont(font)
        self.layout_preset_basis.addWidget(self.basis_title)

        self.layout_preset_basis_grid = QW.QGridLayout()
        names = ['x', 'y', 'z', 'color']
        for n, name in enumerate(names):
            label = QW.QLabel(name)
            label.setAlignment(QC.Qt.AlignCenter)
            self.layout_preset_basis_grid.addWidget(label, 0, n)
        n_coords = 3
        self.preset_basis_coord_widgets = np.empty((n_basis, n_coords),
                                                   dtype=object)
        self.preset_basis_color_widgets = []
        for i in range(n_basis):
            for j in range(n_coords):
                el = QW.QLineEdit()
                el.setEnabled(False)
                if i == 0:
                    el.setText('0')
                self.preset_basis_coord_widgets[i, j] = el
                self.layout_preset_basis_grid.addWidget(el, i + 1, j)
            el = QW.QComboBox()
            el.addItems(self.colors)
            el.activated[str].connect(
                lambda i=i, el=el: self.update_basis_color(
                    'preset', self.preset_basis_color_widgets.index(el), i))
            self.preset_basis_color_widgets.append(el)
            self.layout_preset_basis_grid.addWidget(el, i + 1, n_coords)
        self.layout_preset_basis.addLayout(self.layout_preset_basis_grid)
        self.current_basis_layout = self.layout_preset_basis
        self.layout_options.addLayout(self.layout_preset_basis)

    def create_user_basis(self):
        # Basis-stuff
        font = QG.QFont()
        font.setBold(True)
        self.layout_basis = QW.QVBoxLayout()
        self.basis_title = QW.QLabel('Basis coordinates')
        self.basis_title.setAlignment(QC.Qt.AlignCenter)
        self.basis_title.setFont(font)
        self.layout_basis.addWidget(self.basis_title)
        self.layout_basis_grid = QW.QGridLayout()
        n_basis = 5
        n_coords = 3
        names = ['x', 'y', 'z', 'color']
        for n, name in enumerate(names):
            label = QW.QLabel(name)
            label.setAlignment(QC.Qt.AlignCenter)
            self.layout_basis_grid.addWidget(label, 0, n)
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
                # We want the first to be enabled
                el.setEnabled(i == 0)
                el.setValidator(QG.QDoubleValidator(decimals=2))
                # Pass basis and coordinate number, along with value, to
                # update_basis_val
                el.editingFinished.connect(
                    lambda i=i, j=j, el=el:
                        self.update_basis_val(i, j, el.text()))
                # Add the QLineEdit to the array of basis coordinate widgets
                # and to the layout
                self.basis_coord_widgets[i, j] = el
                self.layout_basis_grid.addWidget(el, i + 1, j)

            # Add a color lineedit for each basis atom
            el = QW.QComboBox()
            el.addItems(self.colors)
            # Okay, for some reason i is the combo-box text. I don't know why.
            # So we're gonna do a slight hack to find the proper "i". We're
            # gonna index the list of color widgets
            el.activated[str].connect(
                lambda i=i, el=el: self.update_basis_color(
                    'user', self.basis_color_widgets.index(el), i))

            self.basis_color_widgets.append(el)
            self.layout_basis_grid.addWidget(el, i + 1, n_coords)

            # Add a checkbox for each basis atom
            check = QW.QCheckBox()
            check.setChecked(i == 0)

            # For some reason this doesn't work when we do it in a loop...
            # check.stateChanged.connect(
            #     lambda i=i: self.hide_basis_widgets(i))

            # Add the checkbox to the list of widgets, and the layout.
            self.basis_check_widgets.append(check)
            self.layout_basis_grid.addWidget(check, i + 1, n_coords + 2)

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

        # We also reset the basis and colors:
        self.lattice_config = self.default_config.copy()

        self.current_basis_layout = self.layout_basis
        self.layout_basis.addLayout(self.layout_basis_grid)
        self.layout_options.addLayout(self.layout_basis)

    def update_lattice(self):
        # Grab a new lattice based on the parameters in lattice_config
        a = self.lattice_config['a']
        b = self.lattice_config['b']
        c = self.lattice_config['c']
        alpha = self.lattice_config['alpha'] * np.pi / 180
        beta = self.lattice_config['beta'] * np.pi / 180
        gamma = self.lattice_config['gamma'] * np.pi / 180
        name = self.lattice_config['lattice']
        (a1, a2, a3), basis, _ = lattices.chooser(lattice_name=name,
                                                  a=a, b=b, c=c,
                                                  alpha=alpha,
                                                  beta=beta,
                                                  gamma=gamma)

        # Update primitive lattice vectors and (preset) basis.
        self.lattice_config.update(dict(zip(('a1', 'a2', 'a3', 'preset_basis'),
                                            (a1, a2, a3, basis))))
        if name in self.presets_with_basis:
            self.update_preset_basis_widgets()
        self.plot_lattice()

    def update_lattice_name(self, text):
        # Delete current basis layout.
        self.delete_layout(self.current_basis_layout)
        if text in self.presets_with_basis:
            # We have a preset with a basis, so we delete the user basis and
            # load the preset basis
            self.create_preset_basis(self.presets_with_basis[text])
            self.current_basis_layout = self.layout_preset_basis
        else:
            self.create_user_basis()
            self.current_basis_layout = self.layout_basis

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

    def update_preset_basis_widgets(self):
        basis = np.atleast_2d(self.lattice_config['preset_basis'])
        for n_atom, atom in enumerate(basis):
            for n_coord, coord in enumerate(atom):
                el = self.preset_basis_coord_widgets[n_atom, n_coord]
                el.setText("{0:.3f}".format(coord))

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
        colors[num] = text
        self.update_basis()

    def plot_lattice(self):
        # This function takes the values from lattice_config and uses them to
        # update the plot.

        # Clear the axes
        self.static_ax.clear()

        # Grab lattice vectors and basis(es) from lattice_config
        a1 = self.lattice_config['a1']
        a2 = self.lattice_config['a2']
        a3 = self.lattice_config['a3']

        # Grab the basis and colors
        if self.lattice_config['lattice'] in self.presets_with_basis:
            # We are dealing with a preset with basis
            basis = self.lattice_config['preset_basis']
            n_basis = np.atleast_2d(basis).shape[0]
            colors = self.lattice_config['preset_colors']
            colors = colors[:n_basis]
        else:
            colors = self.lattice_config['enabled_user_colors']
            basis = self.lattice_config['enabled_user_basis']

        # Plot the new lattice
        self.static_fig, self.static_ax = Lattice(a1=a1, a2=a2, a3=a3,
                                                  basis=basis,
                                                  colors=colors,
                                                  fig=self.static_fig,
                                                  ax=self.static_ax,
                                                  returns=True,
                                                  plots=False,
                                                  checks=False)

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
        new_colors = self.lattice_config['user_colors']
        new_colors = list(compress(new_colors, enabled_basis_atoms))
        self.lattice_config['enabled_user_colors'] = new_colors
        self.plot_lattice()

    def delete_layout(self, layout):
        if layout is not None:
            while layout.count():
                item = layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                else:
                    self.delete_layout(item.layout())


class scattering_window(lattice_window):
    def __init__(self):
        super().__init__()

    def create_variables(self):
        self.lattices = ['simple cubic', 'conventional fcc',
                         'conventional bcc']
        self.colors = [
            'xkcd:cement', 'red', 'blue', 'green', 'cyan',
            'magenta', 'black', 'yellow']
        # A dictionary of the default values for lattice plotting
        self.default_config = {
            'a1': d[0],
            'a2': d[1],
            'a3': d[2],
            'preset_basis': d[3],
            'user_colors': ['xkcd:cement'] * 5,
            'form_factors': [1] * 5,
            'enabled_form_factors': [1],
            'enabled_user_colors': ['xkcd:cement'],
            'preset_colors': ['xkcd:cement'] * 4,
            'sizes': d[5],
            'enabled_user_basis': np.zeros((1, 3)),
            'user_basis': np.zeros((5, 3)),
            'lattice': 'simple cubic',
            'max_preset_basis': 4
        }
        self.presets_with_basis = {
            'wurtzite': 4,
            'diamond': 2,
            'zincblende': 2,
            'conventional fcc': 4,
            'conventional bcc': 2
        }
        self.lattice_config = self.default_config.copy()
        self.current = 'user'

    def create_plot_window(self):
        # Create the default plot and return the figure and axis objects for
        # it. Then create the FigureCanvas, add them all to the layout and add
        # a toolbar. Lastly enable mouse support for Axes3D
        self.static_fig, self.static_ax, self.static_ax2 = Scattering(
            returns=True, plots=False)
        self.static_canvas = FigureCanvas(self.static_fig)
        self.addToolBar(NavigationToolbar(self.static_canvas, self))
        self.static_ax.mouse_init()

    def create_options(self):
        self.layout_options = QW.QVBoxLayout()
        # Create the lattice chooser dropdown
        self.lattice_chooser = QW.QComboBox(self)
        self.lattice_chooser.addItems(self.lattices)
        sep_index = len(self.lattices) - len(self.presets_with_basis)
        self.lattice_chooser.insertSeparator(sep_index)
        self.lattice_chooser.activated[str].connect(self.update_lattice_name)
        self.layout_options.addWidget(self.lattice_chooser)
        self.create_user_basis()
        self.add_form_factors()

    def update_lattice_name(self, text):
        # Delete current basis layout.
        self.delete_layout(self.current_basis_layout)
        if text in self.presets_with_basis:
            # We have a preset with a basis, so we delete the user basis and
            # load the preset basis
            self.create_preset_basis(self.presets_with_basis[text])
        else:
            self.create_user_basis()
        self.lattice_config['lattice'] = text
        self.add_form_factors()

        # And then we update the lattice
        self.update_lattice()

    def update_lattice(self):
        # Grab a new lattice based on the parameters in lattice_config
        a = 1
        name = self.lattice_config['lattice']
        _, basis, _ = lattices.chooser(lattice_name=name, a=a)

        # Update primitive lattice vectors and (preset) basis.
        self.lattice_config['preset_basis'] = basis
        if name in self.presets_with_basis:
            self.update_preset_basis_widgets()
        self.plot_lattice()

    def add_form_factors(self):
        # This method runs whenever a basis has been created, to add the form
        # factors. I do this because then I can reuse as much code as possible
        # First we find out whether we're using a preset or user basis
        place = 4
        if self.lattice_config['lattice'] in self.presets_with_basis:
            lattice_name = self.lattice_config['lattice']
            n_basis = self.presets_with_basis[lattice_name]
            self.current_basis_layout = self.layout_preset_basis_grid
            move_checkboxes = False
        else:
            n_basis = 5
            self.current_basis_layout = self.layout_basis_grid
            move_checkboxes = True

        self.lattice_config['form_factors'] = [1] * n_basis

        self.form_factor_fields = []
        label = QW.QLabel('Form Factors')
        label.setAlignment(QC.Qt.AlignCenter)
        self.current_basis_layout.addWidget(label, 0, place)
        for i in range(n_basis):
            el = QW.QLineEdit()
            el.setText('1')
            if i and move_checkboxes:
                el.setEnabled(False)
            el.setValidator(QG.QDoubleValidator(decimals=2))
            el.editingFinished.connect(
                lambda i=i, el=el: self.update_form_factor(i, el.text()))
            self.form_factor_fields.append(el)
            self.current_basis_layout.addWidget(el, i + 1, place)
            if move_checkboxes:
                el = self.basis_check_widgets[i]
                self.current_basis_layout.addWidget(el, i + 1, place + 1)

    def update_form_factor(self, i, text):
        self.lattice_config['form_factors'][i] = float(text)
        self.update_basis()

    def hide_basis_widgets(self, basis_no):
        # enable or disable basis coord widgets and update the basis
        checkbox = self.basis_check_widgets[basis_no]
        for le in self.basis_coord_widgets[basis_no]:
            le.setEnabled(checkbox.isChecked())
        self.form_factor_fields[basis_no].setEnabled(checkbox.isChecked())
        self.update_basis()

    def update_basis(self):
        # We get a list of basis atoms that are enabled
        enabled_basis_atoms = []
        for i in self.basis_check_widgets:
            enabled_basis_atoms.append(i.isChecked())
        # update the enabled_user_basis config and plot the lattice with the
        # new basis
        new_basis = self.lattice_config['user_basis'][enabled_basis_atoms]
        new_colors = self.lattice_config['user_colors']
        new_colors = list(compress(new_colors, enabled_basis_atoms))
        form_factors = self.lattice_config['form_factors']
        form_factors = list(compress(form_factors, enabled_basis_atoms))
        self.lattice_config['enabled_user_basis'] = new_basis
        self.lattice_config['enabled_user_colors'] = new_colors
        self.lattice_config['enabled_form_factors'] = form_factors
        self.plot_lattice()

    def plot_lattice(self):
        # This function takes the values from lattice_config and uses them to
        # update the plot.

        # Clear the axes
        self.static_ax.clear()
        self.static_ax2.clear()

        # Grab the basis and colors
        if self.lattice_config['lattice'] in self.presets_with_basis:
            # We are dealing with a preset with basis
            basis = self.lattice_config['preset_basis']
            n_basis = np.atleast_2d(basis).shape[0]
            colors = self.lattice_config['preset_colors']
            colors = colors[:n_basis]
            form_factors = self.lattice_config['form_factors']
        else:
            colors = self.lattice_config['enabled_user_colors']
            basis = self.lattice_config['enabled_user_basis']
            form_factors = self.lattice_config['enabled_form_factors']

        # Plot the new lattice
        self.static_fig, self.static_ax, self.static_ax2 = Scattering(
            basis=basis,
            colors=colors,
            form_factor=form_factors,
            fig=self.static_fig,
            axes=(self.static_ax, self.static_ax2),
            returns=True,
            plots=False)

        # Remember to have the canvas draw it!
        self.static_canvas.draw()


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
