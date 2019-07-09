import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from functools import partial
from scipy.spatial import Delaunay

import lattices
import scattering
import band_structure
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
        rounder=True, checks=True, limit=False):
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
        d, planes = lattices.reciprocal(a1, a2, a3, indices, r_min, r_max)
        # planes = lattices.plane_limiter(planes, r_min, r_max)

    if verbose:
        print("Lattice: {}".format(lattice_type))

    # Create the figure
    if fig is None:
        fig = plt.figure()
    if ax is None:
        ax = fig.gca(projection="3d")

    # Get the relevant gridlines:
    pruned_lines = lattices.grid_lines(a1, a2, a3, atomic_positions,
                                       lattice_position, grid_type,
                                       verbose=verbose)
    if grid:
        for line in pruned_lines:
            ax.plot(line[0], line[1], line[2], color=g_col, linewidth=g_w)

    # Plot atoms
    ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
               atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    if indices is not None:
        # If we plot the family of lattice planes, we plot the displacement
        # vector and the planes
        start = a1 + a2 + a3
        ax.quiver(*start, *d)
        ax.text(*(start+d/2), '$d$')

        # If the planes are vertical, then points in xy are colinear, and
        # triangulation cannot happen with these. As such, for plot_trisurf to
        # work we need to triangulate with either x,z or y,z. Here we find out
        # which we want to triangulate with, by accounting for the direction of
        # the displacement vector.
        dhat = d/lattices.mag(d)
        xhat, yhat, zhat = np.eye(3)
        in_xy_plane = eq(dhat.dot(zhat), 0)
        if in_xy_plane:
            # If d is along x, then we triangulate with y, otherwise we
            # triangulate with x
            along_x = eq(dhat.dot(xhat), 1)
            coords = [1, 2] if along_x else [0, 2]

        for points in planes:
            x, y, z = points.T
            if in_xy_plane:
                simp = Delaunay(points[:, coords]).simplices
                ax.plot_trisurf(x, y, simp, z, color='xkcd:cement',
                                shade=False, alpha=0.4)
            else:
                ax.plot_trisurf(x, y, z, color='xkcd:cement', shade=False,
                                alpha=0.4)

    elif arrows:
        # otherwise we plot the lattice vectors
        ax.quiver(0, 0, 0, *a1)
        ax.quiver(0, 0, 0, *a2)
        ax.quiver(0, 0, 0, *a3)
        ax.text(*(a1/2), '$a_1$')
        ax.text(*(a2/2), '$a_2$')
        ax.text(*(a3/2), '$a_3$')

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

    if plots:
        plt.show()

    if returns:
        return fig, ax


Reciprocal = partial(Lattice, indices=(1, 1, 1))


def plot_reciprocal(a1, a2, a3, fig=None, ax=None, indices=(1, 1, 1),
                    grid=False, verbose=False, returns=False, limtype=0):

    if indices is None:
        n_min = np.array([-1, -1, -1])
        n_max = np.array([1, 1, 1])
    else:
        h, k, ell = indices

        # we plot between -h and h with 1 lattice point as padding.
        if limtype == 0:
            n_min = np.array((min(h, -h), min(k, -k), min(ell, -ell))) - 1
            n_max = np.array((max(h, -h), max(k, -k), max(ell, -ell))) + 1
        elif limtype == 1:
            n_min = np.array((min(h, 0), min(k, 0), min(ell, 0))) - 1
            n_max = np.array((max(h, 0), max(k, 0), max(ell, 0))) + 1
        elif limtype == 2:
            n_min = np.array((min(h, -h), min(k, -k), min(ell, -ell)))
            n_max = np.array((max(h, -h), max(k, -k), max(ell, -ell)))
        elif limtype == 3:
            n_min = np.array((min(h, 0), min(k, 0), min(ell, 0)))
            n_max = np.array((max(h, 0), max(k, 0), max(ell, 0)))
    # we want black colors, and small dots
    colors = ['k']
    sizes = [0.5]

    # we use primitive unit cell, and lattice lines along lattice vectors
    unit_type = 'primitive'
    grid_type = 'latticevectors'

    # First the scaling factor for the reciprocal lattice
    scale = a1.dot(np.cross(a2, a3))
    # Then the reciprocal lattice
    b1 = 2 * np.pi * np.cross(a2, a3) / scale
    b2 = 2 * np.pi * np.cross(a3, a1) / scale
    b3 = 2 * np.pi * np.cross(a1, a2) / scale

    # Create the array of lattice vectors and basis
    lattice = np.array([b1, b2, b3])
    basis = np.array([0, 0, 0])

    objects = lattices.generator(b1, b2, b3, basis, colors, sizes,
                                 n_min, n_max)

    atomic_positions = np.around(objects[0], decimals=5)
    objects = [atomic_positions] + [i for i in objects[1:]]

    (atomic_positions, lattice_coefficients, atomic_colors, atomic_sizes,
     lattice_position) = lattices.limiter(points=objects[0],
                                          objects=objects,
                                          min_=n_min,
                                          max_=n_max,
                                          unit_type=unit_type,
                                          lattice=lattice,
                                          verbose=verbose)

    if fig is None:
        recip_fig = plt.figure(figsize=(4, 4))
    else:
        recip_fig = fig

    if ax is None:
        recip_ax = recip_fig.gca(projection="3d")
    else:
        recip_ax = ax

    pruned_lines = lattices.grid_lines(b1, b2, b3, atomic_positions,
                                       lattice_position, grid_type,
                                       verbose=verbose)
    g_col = 'k'
    g_w = 0.5
    if grid:
        for line in pruned_lines:
            recip_ax.plot(line[0], line[1], line[2], color=g_col,
                          linewidth=g_w)

    recip_ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
                     atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    if indices is not None:
        # And the normal vector for the (hkl)-family of planes.
        G = h * b1 + k * b2 + ell * b3
        # plot some arrows
        recip_ax.quiver(0, 0, 0, *G)
        recip_ax.text(*(G/2), f'$G, ({h},{k},{ell})$')

    # making the axes prettier
    recip_ax.set_aspect('equal')
    recip_ax.set_proj_type('ortho')

    # Get the max limit, so we can make plot box cubic
    xlim = recip_ax.get_xlim()
    ylim = recip_ax.get_ylim()
    zlim = recip_ax.get_zlim()
    limits = np.array(list(zip(xlim, ylim, zlim)))
    plot_max = np.max(limits[1])
    plot_min = np.min(limits[0])
    recip_ax.set_xlim([plot_min, plot_max])
    recip_ax.set_ylim([plot_min, plot_max])
    recip_ax.set_zlim([plot_min, plot_max])
    recip_ax.grid(False)
    if not verbose:
        recip_ax.axis('off')

    # make the panes transparent (the plot box)
    recip_ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    recip_ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    recip_ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    if returns:
        return recip_fig, recip_ax
    else:
        plt.show()


def rotatefig(event, fig1, ax1, canvas2, ax2):
    if event.button == 1:
        elev, azim = ax1.elev, ax1.azim
        ax2.view_init(elev, azim)
        canvas2.draw()


def Scattering(lattice_name='simple cubic',
               basis=None,
               k_in=np.array([0, 0, -1.5]),
               form_factor=None,
               highlight=None,
               show_all=True,
               normalize=True,
               verbose=False,
               returns=False,
               return_indices=False,
               colors=None,
               laue_scale=1,
               figs=None,
               axes=None,
               plots=True):

    min_, max_ = (-2, -2, -1), (2, 2, 1)
    g_col = 'k'
    g_w = 0.5
    g_a = 0.6
    size_default = 36
    point_sizes = 2
    point_sizes *= size_default
    plane_z = 25
    beam_end_z = max_[2]
    unit_cell_type = "conventional"
    lim_type = "proper"
    outgoing_length = 10

    # input sanitization for the lattice/basis
    lattice_name = lattice_name.lower()
    if basis is not None:
        a1, a2, a3 = np.eye(3, dtype=int)
        basis = np.array(basis)
    else:
        lattice_name = lattice_name.lower()
        if lattice_name == "bcc":
            lattice_name = "conventional bcc"
        elif lattice_name == "fcc":
            lattice_name = "conventional fcc"
        elif lattice_name == "simple cubic":
            pass
        else:
            print("Allowed inputs: 'simple cubic', 'bcc', 'fcc'.")
            return
        lattice, basis, _ = lattices.chooser(lattice_name, verbose=verbose)
        a1, a2, a3 = lattice

    grid_type = lattices.latticelines[lattice_name]

    # Getting the number of atoms in the basis
    length_basis = np.shape(basis)
    if len(length_basis) == 1:
        n_basis = 1
    elif len(length_basis) > 1:
        n_basis = length_basis[0]

    if form_factor is None:
        form_factor = [1] * n_basis

    if colors is None:
        form_fact_array = np.array(form_factor)
        if (form_fact_array == form_factor[0]).all():
            atom_colors = ["xkcd:cement"] * n_basis
        else:
            atom_colors = ["xkcd:cement"] + (["xkcd:cornflower blue"] *
                                             (n_basis - 1))
    else:
        atom_colors = colors

    atom_sizes = [1] * n_basis

    # Normalizing wave vector (multiplying by k0 = 2Pi/a)
    k_in = np.array(k_in)

    if k_in[2] > 0:
        k_in[2] = -k_in[2]
        print(("the z-coordinate of k_in should be negative. "
               "Flipping it: k_in = {}".format(k_in)))

    k_title = np.copy(k_in)
    if normalize:
        k_in = k_in * 2 * np.pi

    # Calculating stuff for plotting the crystal
    r_min, r_max, n_min, n_max = lattices.find_limits(lim_type, a1, a2, a3,
                                                      min_, max_,
                                                      unit_cell_type)
    objects = lattices.generator(a1, a2, a3, basis, atom_colors,
                                 atom_sizes, n_min, n_max)
    objects = lattices.limiter(objects[0], objects, r_min, r_max,
                               unit_cell_type)
    (atomic_positions, _, atomic_colors, atomic_sizes,
     lattice_position) = objects

    pruned_lines = lattices.grid_lines(a1, a2, a3, atomic_positions,
                                       lattice_position, grid_type,
                                       verbose=verbose)

    # Create the neutron beam display vector
    k_disp = k_in / lattices.mag(k_in)
    lambda_ = 2 * np.pi / lattices.mag(k_in)

    # Scattering stuff
    intensities, k_out, indices = scattering.calc_scattering(a1, a2, a3, basis,
                                                             form_factor,
                                                             k_in)
    points = scattering.projection(k_out, p0=np.array([0, 0, plane_z]))

    # Plotting the basic
    if figs is None:
        macro_fig = plt.figure(figsize=(8, 8))
        micro_fig = plt.figure(figsize=(5, 5))
    else:
        macro_fig, micro_fig = figs
    if axes is None:
        macro_ax = macro_fig.add_axes([0.05, 0.05, 0.9, 0.9], projection='3d')
        micro_ax = micro_fig.add_axes([0.05, 0.05, 0.9, 0.9], projection='3d')
    else:
        macro_ax, micro_ax = axes

    # Plot atoms
    macro_ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
                     atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    micro_ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
                     atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    for line in pruned_lines:
        macro_ax.plot(line[0], line[1], line[2],
                      color=g_col, linewidth=g_w, alpha=g_a)
        micro_ax.plot(line[0], line[1], line[2],
                      color=g_col, linewidth=g_w, alpha=g_a)

    # Plotting the beam: First we create the beam display vector
    macro_ax.quiver(0, 0, beam_end_z, k_disp[0], k_disp[1], k_disp[2],
                    color='b', lw=2, pivot='tip', length=lambda_ * laue_scale)
    micro_ax.quiver(0, 0, beam_end_z, k_disp[0], k_disp[1], k_disp[2],
                    color='b', lw=2, pivot='tip', length=lambda_ * laue_scale)

    if intensities.size == 0:
        print("There is no scattering for this choice of k_in")

    else:
        # I assume the points are unique, now that I have deleted the ones
        # pointing into the crystal

        # Normalize intensities
        intensities /= np.amax(intensities)
        # Create the color array
        colors = np.zeros((intensities.size, 4))
        colors[:, 3] = intensities

        if highlight is not None:
            # Checking for proper highlighting
            hi_index = np.array(highlight)
            num_ints = hi_index.shape
            extra = 0
            if num_ints != (3,):
                print("We need 3 and only 3 indices! Highlighting nothing")
            else:
                indices_index = np.where((indices == hi_index).all(axis=1))[0]
                if indices_index.shape != (1,):
                    print("There is no scattering along {}".format(highlight))
                else:
                    # We have highlighting!
                    d, planes = lattices.reciprocal(a1, a2, a3, hi_index,
                                                    r_min - extra,
                                                    r_max + extra,
                                                    points=20)
                    planes = lattices.plane_limiter(planes, r_min - extra,
                                                    r_max + extra)
                    # We change the color of highlighted point and plot the
                    # family of planes
                    high_intensity = intensities[indices_index]
                    colors[indices_index] = [1, 0, 0, high_intensity]
                    for p in planes:
                        micro_ax.plot_surface(p[0], p[1], p[2], color="r",
                                              shade=False, alpha=0.2)
                    # We also plot the outgoing line corresponding to this
                    # scattering. First we get the point (and squeeze it, to
                    # make it 1D again)
                    p = np.squeeze(points[indices_index, :])
                    start = np.array([0, 0, beam_end_z])
                    ray = p - start
                    line = np.array([start, start + outgoing_length * ray])
                    macro_ax.plot(line[:, 0], line[:, 1], line[:, 2],
                                  color='r',
                                  alpha=0.3, ls='--',
                                  lw=g_w * 2)
                    # Plotting outgoing vector and Laue condition
                    k_out_high = np.squeeze(k_out[indices_index])
                    G_high = k_in - k_out_high
                    vecs = np.array([k_out_high, k_out_high, G_high])
                    vecs_disp = vecs / lattices.mag(k_in)
                    starts = np.array([start,
                                       start - k_disp * lambda_ * laue_scale,
                                       (start - vecs_disp[2] * lambda_ *
                                        laue_scale)])
                    micro_ax.quiver(starts[:, 0],
                                    starts[:, 1],
                                    starts[:, 2],
                                    vecs_disp[:, 0],
                                    vecs_disp[:, 1],
                                    vecs_disp[:, 2],
                                    color=['r', 'r', 'g'],
                                    alpha=0.5,
                                    lw=1,
                                    length=lambda_ * laue_scale)

        ranges = (np.amax(points, axis=0) - np.amin(points, axis=0))[:-1]

        # Plotting detection plane
        abs_extra = 0.0
        rel_extra = 0.1
        def_ = 2
        x_min = np.amin(points[:, 0]) * (1 + rel_extra) - abs_extra
        x_max = np.amax(points[:, 0]) * (1 + rel_extra) + abs_extra
        y_min = np.amin(points[:, 1]) * (1 + rel_extra) - abs_extra
        y_max = np.amax(points[:, 1]) * (1 + rel_extra) + abs_extra

        # We want the detection plane to be square:
        max_ = max(x_max, y_max)
        min_ = min(x_min, y_min)

        x_range = np.array([min(min_, -def_), max(max_, def_)])
        y_range = np.array([min(min_, -def_), max(max_, def_)])
        x, y = np.meshgrid(x_range, y_range)
        z = plane_z * np.ones(x.shape)
        macro_ax.plot_surface(x, y, z, color='k', alpha=0.2)

        # plotting intersections
        macro_ax.scatter(points[:, 0], points[:, 1], plane_z, color=colors)
        for i in range(len(indices)):
            x, y = points[i, 0:2]
            s = indices[i]
            c = colors[i, :-1]
            macro_ax.text(x, y, plane_z, s, color=c, va='bottom', ha='center')

        # Setting limits for the second figure
        det_max_x = np.amax(x_range)
        det_min_x = np.amin(x_range)
        det_max_y = np.amax(y_range)
        det_min_y = np.amin(y_range)
        det_max = max(det_max_x, det_max_y)
        det_min = min(det_min_x, det_min_y)
        # ax2.set_xlim(det_min, det_max)
        # ax2.set_ylim(det_min, det_max)

        if show_all:
            # Plotting outgoing vectors
            n = k_out.shape[0]
            k_plot = k_out / lattices.mag(k_in)
            start_point = np.array((0, 0, beam_end_z))
            start_points = np.repeat(np.atleast_2d(start_point), n, axis=0)
            macro_ax.quiver(start_points[:, 0],
                            start_points[:, 1],
                            start_points[:, 2],
                            k_plot[:, 0],
                            k_plot[:, 1],
                            k_plot[:, 2],
                            color='g',
                            alpha=0.5,
                            lw=g_w,
                            length=lambda_)

            # plotting outgoing lines
            for p in points:
                ray = p - start_point
                line = np.array([start_point,
                                 start_point + outgoing_length * ray])
                macro_ax.plot(line[:, 0], line[:, 1], line[:, 2],
                              color='k', alpha=0.3, ls='--',
                              lw=g_w)

    macro_ax.set_aspect('equal')
    macro_ax.set_proj_type('ortho')
    macro_ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    macro_ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    macro_ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    micro_ax.set_aspect('equal')
    micro_ax.set_proj_type('ortho')
    micro_ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    micro_ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    micro_ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Some limit trickery. We make the plot box cubic:
    plot_max = np.amax(r_max)
    plot_min = np.amin(r_min)
    plane_max = max_
    plane_min = min_
    macro_max = max([plot_max, plane_max, plane_z])
    macro_min = min([plot_min, plane_min, plane_z])
    macro_ax.set_xlim(macro_min, macro_max)
    macro_ax.set_ylim(macro_min, macro_max)
    macro_ax.set_zlim(macro_min, macro_max)

    macro_ax.grid(False)
    macro_ax.axis('off')

    micro_ax.set_xlim(plot_min, plot_max)
    micro_ax.set_ylim(plot_min, plot_max)
    micro_ax.set_zlim(plot_min, plot_max)

    micro_ax.grid(False)
    micro_ax.axis('off')

    tit = (r'Scattering on a cubic lattice. $k_{in} = (2\pi/a)\cdot$' +
           '{}'.format(k_title) + f'\nForm factors: {form_factor}')
    tit2 = (r'Scattering on a cubic lattice. $k_{in} = $' +
            '{}'.format(k_title))
    if normalize:
        macro_ax.set_title(tit)
    else:
        macro_ax.set_title(tit2)
    if plots:
        plt.show()
    return_list = []
    if returns:
        return_list += [macro_fig, micro_fig, macro_ax, micro_ax]
    if return_indices:
        return_list.append(indices)
    if returns or return_indices:
        return return_list


def Band_structure(V0=0, n_k=51, G_range=list(range(-3, 4)),
                   potential="harmonic", edges=False,
                   E_F=None,
                   returns=False,
                   plots=True):

    # First some input sanitization
    potentials = {"harmonic": band_structure.VG_cos,
                  "dirac": band_structure.VG_dirac}
    class_name = potential.__class__.__name__
    if class_name not in ["function", "str"]:
        print(("Please input either 'dirac', 'harmonic' or the name of your "
               "own potential function. Giving harmonic."))
        potential = band_structure.VG_cos
    elif class_name == "str":
        if potential.lower() not in ["dirac", "harmonic"]:
            print(("Please input either 'dirac' or 'harmonic' as the "
                   "potential. Giving harmonic."))
            potential = band_structure.VG_cos
        potential = potentials[potential.lower()]

    # We calculate the band structure
    o = band_structure.calc_band_structure(V0=V0, n_k=n_k, G_range=G_range,
                                           potential=potential)
    kxs, kys, band, max_E = o
    if E_F is None:
        E_F = max_E
    E_F_mat = max_E * np.ones((n_k, n_k))
    max_k = np.amax(kxs)
    min_k = np.amin(kxs)

    # Create the figure
    fig = plt.figure(figsize=(10, 4))
    ax = fig.gca(projection="3d")
    ax.set_position([0.05, 0, 0.5, 1])

    # Optional plotting of Fermi surface in main axes.
    ax.contour(kxs, kys, band, E_F, colors='r', linewidths=3)

    # Plotting of the main event: the band structure
    ax.plot_surface(kxs, kys, band, alpha=0.9)
    ax.plot_surface(kxs, kys, E_F_mat, alpha=0.2)
    ax.set_xlim([min_k, max_k])
    ax.set_ylim([min_k, max_k])
    ax.set_xlabel(r'$k_x/k_0$')
    ax.set_ylabel(r'$k_y/k_0$')
    ax.set_zlabel(r'$E/E_0$')
    ax.set_title(('Band structure of square lattice. ' +
                  '$V_0/E_0 = {}$. $E_F = {}$'.format(V0, np.round(E_F, 3))))

    # optional plotting of the edges
    if edges:
        # First we get the edges and wave vectors
        edge1 = band[0, :]
        edge2 = band[-1, :]
        edge3 = band[:, 0]
        edge4 = band[:, -1]
        k = kxs[0, :]
        start = -1 / 2 * np.ones(k.shape)
        end = 1 / 2 * np.ones(k.shape)

        # And plot them.
        ax.plot(k, start, edge1, c='k')
        ax.plot(k, end, edge2, c='k')
        ax.plot(start, k, edge3, c='k')
        ax.plot(end, k, edge4, c='k')

    # Plotting of the second set of axes
    ax2 = plt.axes([0.7, 0.2, 0.25, 0.6])
    ax2.contour(kxs, kys, band, E_F)
    ax2.set_xlabel(r'$k_x/k_0$')
    ax2.set_ylabel(r'$k_y/k_0$')
    ax2.set_title('Fermi surface')

    if plots:
        plt.show()

    if returns:
        return fig, ax, ax2


if __name__ == "__main__":
    gui.main()
