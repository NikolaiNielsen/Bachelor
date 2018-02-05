# Let's get the importing out of the way.
# Numpy for calculations, matplotlib for plotting.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# Make the matplotlib plots interactive
get_ipython().magic('matplotlib notebook')


def LatticeGenerator(a1, a2, a3, basis, colors='xkcd:cement', sizes=1):
    size_default = 36

    # We need the number of basis-vectors.
    # If there is only 1 basis vector, then len(np.shape(basis)) == 1
    # otherwise the length is 2, and the first element is number of basis
    # vectors
    length_basis = np.shape(basis)
    if len(length_basis) == 1:
        N_basis = 1
    elif len(length_basis) > 1:
        N_basis = length_basis[0]

    # Make a list, N_basis long, for the colors and sizes,
    # if they're not specified.
    c_name = colors.__class__.__name__
    if c_name == "str":
        c = colors
        colors = []
        for i in range(N_basis):
            colors.append(c)
    elif c_name == "list" and len(colors) != N_basis:
        c = colors[0]
        colors = []
        for i in range(N_basis):
            colors.append(c)

    s_name = sizes.__class__.__name__
    if s_name == "int" or s_name == "float":
        s = sizes
        sizes = []
        for i in range(N_basis):
            sizes.append(s)
    elif s_name == "list" and len(sizes) != N_basis:
        s = sizes[0]
        sizes = []
        for i in range(N_basis):
            sizes.append(s)

    # set the range of lattice vectors to be calculated
    nx_max, ny_max, nz_max = 4, 4, 4
    nx_min, ny_min, nz_min = -2, -2, -2
    # Calculate the amount of atomic positions to be calculated
    numAtoms = ((nx_max + 1 - nx_min) * (ny_max + 1 - ny_min) *
                (nz_max + 1 - nz_min) * N_basis)
    # Make a zero array for all of the atomic positions
    AtomicPositions = np.zeros((numAtoms, 3))
    AtomicColors = []
    AtomicSizes = []

    # Loop over all chosen linear combinations of basis vectors and plot each
    counter = 0
    for nx in range(nx_min, nx_max + 1):
        for ny in range(ny_min, ny_max + 1):
            for nz in range(nz_min, nz_max + 1):
                lattice_position = nx * a1 + ny * a2 + nz * a3

                for n_atom in range(N_basis):
                    AtomicPositions[counter, ] = (lattice_position +
                                                  basis[n_atom, ])
                    AtomicColors.append(colors[n_atom])
                    AtomicSizes.append(size_default * sizes[n_atom])
                    counter += 1

    # We pass splices with all elements for Colors and Sizes, because they
    # change size with the function call
    plotter(a1, a2, a3, AtomicPositions, AtomicColors[:], AtomicSizes[:])
    plotter(a1, a2, a3, AtomicPositions, AtomicColors[:], AtomicSizes[:],
            x_min=0, x_max=1,
            y_min=0, y_max=1,
            z_min=0, z_max=1)


def limiter(AtomicPositions, AtomicColors, AtomicSizes,
            x_min=0, x_max=2,
            y_min=0, y_max=2,
            z_min=0, z_max=2):
    # We need a way to hide points outside the limits
    # The limits
    atoms, dims = np.shape(AtomicPositions)
    rows = []

    # loop over all row ID's
    for rowID in range(atoms):
        # if the atom is outside the limits, we append the row ID to a list
        row = AtomicPositions[rowID, ]
        # (it's actually easier and prettier to check if they're inside)
        inside_x = x_min <= row[0] <= x_max
        inside_y = y_min <= row[1] <= y_max
        inside_z = z_min <= row[2] <= z_max
        inside = inside_x and inside_y and inside_z
        if not inside:
            rows.append(rowID)

    # delete all rows (axis 0 of the array) that are outside the limits
    AtomicPositions = np.delete(AtomicPositions, rows, 0)
    for ID in sorted(rows, reverse=True):
        del AtomicColors[ID]
        del AtomicSizes[ID]

    return AtomicPositions, AtomicColors, AtomicSizes


def plotter(a1, a2, a3, AtomicPositions, AtomicColors, AtomicSizes,
            x_min=0, x_max=2,
            y_min=0, y_max=2,
            z_min=0, z_max=2):
    # Create the figure
    A1 = limiter(AtomicPositions, AtomicColors, AtomicSizes,
                 x_min, x_max, y_min, y_max, z_min, z_max)

    AtomicPositions, AtomicColors, AtomicSizes = A1

    fig = plt.figure()
    ax = fig.gca(projection="3d")

    # Plot atoms. For now a single size and color
    ax.scatter(AtomicPositions[:, 0], AtomicPositions[:, 1],
               AtomicPositions[:, 2], c=AtomicColors,
               s=AtomicSizes)

    # Create grid lines
    g_col = 'k'
    g_w = 0.5
    for nx in range(int(np.ceil(x_min)), int(np.floor(x_max)) + 1):
        for ny in range(int(np.ceil(y_min)), int(np.floor(y_max)) + 1):
            ax.plot(np.array([nx, nx]),
                    np.array([ny, ny]),
                    np.array([np.ceil(z_min), np.floor(z_max)]),
                    c=g_col, linewidth=g_w)

        for nz in range(int(np.ceil(z_min)), int(np.floor(z_max)) + 1):
            ax.plot(np.array([nx, nx]),
                    np.array([np.ceil(z_min), np.floor(z_max)]),
                    np.array([nz, nz]),
                    c=g_col, linewidth=g_w)

    for ny in range(int(np.ceil(y_min)), int(np.floor(y_max)) + 1):
        for nz in range(int(np.ceil(z_min)), int(np.floor(z_max)) + 1):
            ax.plot(np.array([np.ceil(x_min), np.floor(x_max)]),
                    np.array([ny, ny]),
                    np.array([nz, nz]),
                    c=g_col, linewidth=g_w)

    # plot lattice vectors
    ax.quiver(0, 0, 0, a1[0], a1[1], a1[2])
    ax.quiver(0, 0, 0, a2[0], a2[1], a2[2])
    ax.quiver(0, 0, 0, a3[0], a3[1], a3[2])
    ax.text(a1[0] / 2, a1[1] / 2, a1[2] / 2, '$a_1$')
    ax.text(a2[0] / 2, a2[1] / 2, a2[2] / 2, '$a_2$')
    ax.text(a3[0] / 2, a3[1] / 2, a3[2] / 2, '$a_3$')

    # Set limits, orthographic projection (so we get the beautiful hexagons),
    # no automatic gridlines, ticks and labels
    ax.set_xlim([x_min, x_max])
    ax.set_ylim([y_min, y_max])
    ax.set_zlim([z_min, z_max])
    ax.set_proj_type('ortho')
    ax.grid(False)
    ax.set_xticks(list(range(x_min, x_max + 1)))
    ax.set_yticks(list(range(y_min, y_max + 1)))
    ax.set_zticks(list(range(z_min, z_max + 1)))
    ax.set_xlabel('$x$')
    ax.set_ylabel('$y$')
    ax.set_zlabel('$z$')
