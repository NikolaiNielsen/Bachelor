"""
Next:
Reciprocal lattice and lattice planes
Scattering
- Parameters:
  - Incident angle
  - Incident wavelength of light
  - Start with simple cubic, and 90 degree incident angle (specific case)
Distortion of Fermi Surface:
- Transform Hamiltonian to k-space. Transforms from a matrix to just a number
"""

import itertools

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lattices
import scattering

eq = np.isclose

d = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
     np.array([0, 0, 0]), "xkcd:cement", 2, "proper", "latticevectors",
     [0, 0, 0], [2, 2, 2])


def Lattice(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, min_=d[8], max_=d[9], lattice_name=None,
        indices=None, verbose=False):
    """
    Creates, limits and plots the lattice
    """
    if lattice_name is not None:
        lattice, basis, lattice_type = lattices.chooser(lattice_name,
                                                        verbose=verbose)
        a1, a2, a3 = lattice
        # Classify the lattice
    else:
        lattice_type = lattices.classifier(a1, a2, a3, basis)

        # Rotate the lattice
        a1, a2, a3, basis = lattices.rotator(a1, a2, a3, basis,
                                             lattice_type, verbose=verbose)

    # Input sanitization:
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

    # Choosing gridline type. First the default settings.
    latticelines = lattices.latticelines

    # This isn't really pretty, but it works. If the Gridtype is already set,
    # we use that value. If not, and the latticetype is undefined we choose the
    # default type (latticevectors). Otherwise (when the latticetype is
    # defined) we choose the gridtype appropriate for the lattice
    if grid_type is not None:
        pass
    else:
        grid_type = latticelines[lattice_type]
    # set the range of lattice vectors to be calculated
    r_min, r_max, n_min, n_max = lattices.find_limits(lim_type, a1, a2, a3,
                                                      min_, max_)

    (atomic_positions, lattice_coefficients, atomic_colors, atomic_sizes,
     lattice_position) = lattices.generator(a1, a2, a3, basis, colors, sizes,
                                            lim_type, n_min, n_max, r_min,
                                            r_max)
    # Objects to limit to the plot-box
    objects = [atomic_positions, lattice_coefficients, atomic_colors,
               atomic_sizes, lattice_position]
    objects = lattices.limiter(atomic_positions, objects, r_min, r_max)
    (atomic_positions, lattice_coefficients, atomic_colors, atomic_sizes,
     lattice_position) = objects

    if indices is not None:
        if len(indices) != 3:
            print("We need 3 indices! We'll give you (1,1,1)")
            indices = (1, 1, 1)
        d, planes = lattices.reciprocal(a1, a2, a3, indices, r_min, r_max)
        planes = lattices.plane_limiter(planes, r_min, r_max)

    if verbose:
        print("Lattice: {}".format(lattice_type))

    # Create the figure
    fig = plt.figure()
    ax = fig.gca(projection="3d")

    # Plot atoms
    ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
               atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    # Get the relevant gridlines:
    g_col = 'k'
    g_w = 0.5
    pruned_lines = lattices.grid_lines(a1, a2, a3, atomic_positions,
                                       lattice_position, grid_type,
                                       verbose=verbose)
    for line in pruned_lines:
        ax.plot(line[0], line[1], line[2], color=g_col, linewidth=g_w)

    if indices is not None:
        # If we plot the family of lattice planes, we plot the displacement
        # vector and the planes
        ax.quiver(0, 0, 0, d[0], d[1], d[2])
        ax.text(d[0] / 2, d[1] / 2, d[2] / 2, '$d$')
        for p in planes:
            ax.plot_surface(p[0], p[1], p[2], color='xkcd:cement', shade=False,
                            alpha=0.4)
    else:
        # otherwise we plot the lattice vectors
        ax.quiver(0, 0, 0, a1[0], a1[1], a1[2])
        ax.quiver(0, 0, 0, a2[0], a2[1], a2[2])
        ax.quiver(0, 0, 0, a3[0], a3[1], a3[2])
        ax.text(a1[0] / 2, a1[1] / 2, a1[2] / 2, '$a_1$')
        ax.text(a2[0] / 2, a2[1] / 2, a2[2] / 2, '$a_2$')
        ax.text(a3[0] / 2, a3[1] / 2, a3[2] / 2, '$a_3$')

    # Set limits, orthographic projection (so we get the beautiful hexagons),
    # no automatic gridlines, and no axes
    ax.set_aspect('equal')
    ax.set_proj_type('ortho')
    ax.set_xlim([r_min[0], r_max[0]])
    ax.set_ylim([r_min[1], r_max[1]])
    ax.set_zlim([r_min[2], r_max[2]])
    ax.grid(False)
    plt.axis('equal')
    plt.axis('off')

    # make the panes transparent (the plot box)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    plt.show()


def Reciprocal(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, min_=d[8], max_=d[9], lattice_name=None,
        indices=(1, 1, 1), verbose=False):
    Lattice(a1, a2, a3, basis, colors, sizes, lim_type, grid_type, min_, max_,
            lattice_name, indices, verbose)
