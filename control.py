import numpy as np
import lattices


def lattice_input_sanitization(a1, a2, a3, basis, colors, sizes, grid_type,
                               unit_type, lattice_name,
                               max_, min_, lim_type, verbose):
    # Input sanitization:

    if lattice_name is not None:
        lattice, basis, lattice_type = lattices.chooser(lattice_name,
                                                        verbose=verbose)
        a1, a2, a3 = lattice
        # Classify the lattice
    else:
        a1, a2, a3 = np.array([a1, a2, a3])
        basis = np.array(basis)
        lattice_type = lattices.classifier(a1, a2, a3, basis)

        # Rotate the lattice
        a1, a2, a3, basis = lattices.rotator(a1, a2, a3, basis,
                                             lattice_type, verbose=verbose)
    lattice = np.array([a1, a2, a3])

    # We need the number of basis-vectors.
    # If there is only 1 basis vector, then len(np.shape(basis)) == 1
    # otherwise the length is 2, and the first element is number of basis
    # vectors
    length_basis = np.shape(basis)
    if len(length_basis) == 1:
        n_basis = 1
    elif len(length_basis) > 1:
        n_basis = length_basis[0]

    # Make a list, n_basis long, for the colors and sizes,
    # if they're not specified.
    c_name = colors.__class__.__name__
    if c_name == "str":
        c = colors
        colors = []
        for i in range(n_basis):
            colors.append(c)
    elif c_name == "list" and len(colors) < n_basis:
        c = colors[0]
        colors = []
        for i in range(n_basis):
            colors.append(c)

    s_name = sizes.__class__.__name__
    if s_name == "int" or s_name == "float":
        s = sizes
        sizes = []
        for i in range(n_basis):
            sizes.append(s)
    elif s_name == "list" and len(sizes) < n_basis:
        s = sizes[0]
        sizes = []
        for i in range(n_basis):
            sizes.append(s)

    # Choosing gridline and unit cell type. First the default settings.
    latticelines = lattices.latticelines
    unitcells = lattices.unitcells

    if grid_type is None:
        grid_type = latticelines[lattice_type]

    if unit_type is None:
        unit_type = unitcells[lattice_type]
    else:
        try:
            unit_type = unit_type.lower()
        except AttributeError:
            print('Please input a string for unit_type. Giving primitive')
            unit_type = "primitive"

        if unit_type not in ["primitive", "conventional"]:
            print(("Input either 'primitive' or 'conventional' for type."
                   " Giving 'primitive'"))
            unit_type = "primitive"

    # Ugly hack for fixing bcc involving juggling of limits, so we plot 2 unit
    # cells (conventional) in each direction
    if (lattice_name == "bcc" and max_ == [2, 2, 2] and
            lim_type == "proper" and unit_type == "conventional"):
        max_ = [0, 0, 4]

    return (a1, a2, a3, lattice, basis, colors, sizes, grid_type, unit_type,
            lattice_type, lattice_name, max_, min_, lim_type)


def lattice_plotting(fig, ax, ):
    pass
