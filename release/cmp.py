import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lattices
import gui

d = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
     np.array([0, 0, 0]), "xkcd:cement", 2, "proper", "latticevectors",
     [2, 2, 2])


eq = np.isclose


def Lattice(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, max_=d[8], lattice_name=None,
        unit_type=None, indices=None, arrows=True, grid=True,
        verbose=False, returns=False, fig=None, ax=None, plots=True,
        rounder=True, checks=True):
    """
    Creates, limits and plots the lattice
    """

    # Minimum coefficients for lattice vectors
    min_ = [0, 0, 0]
    # Settings for grid lines
    g_col = 'k'
    g_w = 0.5
    # Number of points per plane
    num_plane_points = 20

    # Input sanitization
    if lattice_name is not None:
        lattice, basis, lattice_type = lattices.chooser(lattice_name,
                                                        verbose=verbose)
        a1, a2, a3 = lattice
        # Classify the lattice
    else:
        a1, a2, a3 = np.array([a1, a2, a3])
        basis = np.array(basis)
        lattice_type = lattices.classifier(a1, a2, a3, basis)

        if checks:
            # Rotate the lattice, but only if we actually need to check it (for
            # example when we're not dealing with the gui)
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
    if (lattice_type in ['bcc', 'tetragonal body centred',
                         'orthorhombic body centred'] and
        max_ == [2, 2, 2] and lim_type == "proper" and (
            unit_type == "conventional")):
        max_ = [0, 0, 4]
    elif ('base centred' in lattice_type and max_ == [2, 2, 2] and
          lim_type == "proper" and unit_type == "conventional"):
        max_ = [0, 4, 2]

    # set the range of lattice vectors to be calculated
    r_min, r_max, n_min, n_max = lattices.find_limits(lim_type, a1, a2, a3,
                                                      min_, max_,
                                                      unit_type=unit_type)
    if rounder:
        r_min, r_max = np.around([r_min, r_max], decimals=5)
    if verbose:
        print("Limits found as: (type: {})".format(lim_type))
        print("r_min, r_max, n_min, n_max:")
        print(r_min, r_max, n_min, n_max)

    # if we plot the conventional cell we want to give r_min and r_max to
    # limiter. If not we want to give n_min and n_max
    if unit_type == "conventional":
        lim_min, lim_max = r_min, r_max
    else:
        lim_min, lim_max = n_min, n_max

    objects = lattices.generator(a1, a2, a3, basis, colors, sizes,
                                 n_min, n_max)
    if rounder:
        atomic_positions = np.around(objects[0], decimals=5)
        objects = [atomic_positions] + [i for i in objects[1:]]
    if verbose:
        print("Number of atoms and lattice points before limiting:")
        print(objects[0].size / 3, np.sum(objects[-1]))
    # Objects to limit to the plot-box
    (atomic_positions, lattice_coefficients, atomic_colors, atomic_sizes,
     lattice_position) = lattices.limiter(points=objects[0],
                                          objects=objects,
                                          min_=lim_min,
                                          max_=lim_max,
                                          unit_type=unit_type,
                                          lattice=lattice,
                                          verbose=verbose)
    if verbose:
        print("Number of atoms and lattice points AFTER limiting:")
        print(objects[0].size / 3, np.sum(objects[-1]))

    if indices is not None:
        if len(indices) != 3:
            print("We need 3 indices! We'll give you (1,1,1)")
            indices = (1, 1, 1)
        d, planes = lattices.reciprocal(a1, a2, a3, indices, r_min, r_max,
                                        points=num_plane_points)
        planes = lattices.plane_limiter(planes, r_min, r_max)

    if verbose:
        print("Lattice: {}".format(lattice_type))

    # Create the figure
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.gca(projection="3d")

    # Plot atoms
    ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
               atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    # Get the relevant gridlines:
    pruned_lines = lattices.grid_lines(a1, a2, a3, atomic_positions,
                                       lattice_position, grid_type,
                                       verbose=verbose)
    if grid:
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
    elif arrows:
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

    plot_max = np.amax(r_max)
    plot_min = np.amin(r_min)
    ax.set_xlim([plot_min, plot_max])
    ax.set_ylim([plot_min, plot_max])
    ax.set_zlim([plot_min, plot_max])
    ax.grid(False)
    if not verbose:
        ax.axis('off')

    # make the panes transparent (the plot box)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    if returns:
        return fig, ax

    if plots:
        plt.show()


if __name__ == "__main__":
    gui.main()
