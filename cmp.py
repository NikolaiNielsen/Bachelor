import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lattices
import scattering
import band_structure
import control

eq = np.isclose

d = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]), 
     np.array([0, 0, 0]), "xkcd:cement", 2, "proper", "latticevectors",
     [2, 2, 2])


def Lattice(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, max_=d[8], lattice_name=None,
        unit_type=None, indices=None, arrows=True, grid=True,
        verbose=False, returns=False, fig=None, ax=None, plots=True):
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
    (a1, a2, a3, lattice, basis, colors, sizes, grid_type, unit_type,
     lattice_type, lattice_name, max_, min_,
     lim_type) = control.lattice_input_sanitization(a1, a2, a3,
                                                    basis, colors,
                                                    sizes,
                                                    grid_type,
                                                    unit_type,
                                                    lattice_name,
                                                    max_, min_,
                                                    lim_type,
                                                    verbose)
    # set the range of lattice vectors to be calculated
    r_min, r_max, n_min, n_max = lattices.find_limits(lim_type, a1, a2, a3,
                                                      min_, max_,
                                                      unit_type=unit_type)
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


def Reciprocal(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, max_=d[8], lattice_name=None,
        unit_type=None, indices=(1, 1, 1), arrows=True, grid=True,
        verbose=False, returns=False):
    if returns:
        fig, ax = Lattice(a1=a1, a2=a2, a3=a3,
                          basis=basis,
                          colors=colors,
                          sizes=sizes,
                          lim_type=lim_type,
                          grid_type=grid_type,
                          max_=max_,
                          unit_type=unit_type,
                          lattice_name=lattice_name,
                          indices=indices,
                          arrows=arrows,
                          grid=grid,
                          verbose=verbose,
                          returns=True)
        return fig, ax
    else:
        Lattice(a1=a1, a2=a2, a3=a3,
                basis=basis,
                colors=colors,
                sizes=sizes,
                lim_type=lim_type,
                grid_type=grid_type,
                max_=max_,
                unit_type=unit_type,
                lattice_name=lattice_name,
                indices=indices,
                arrows=arrows,
                grid=grid,
                verbose=verbose)


def Scattering(lattice_name='simple cubic',
               basis=None,
               k_in=np.array([0, 0, -1.5]),
               form_factor=None,
               highlight=None,
               show_all=False,
               normalize=True,
               verbose=False,
               returns=False,
               colors=None,
               laue_scale=1):

    min_, max_ = (-2, -2, -1), (2, 2, 1)
    g_col = 'k'
    g_w = 0.5
    g_a = 0.6
    size_default = 36
    point_sizes = 2
    point_sizes *= size_default
    plane_z = 3.5
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
    (atomic_positions, lattice_coefficients, atomic_colors, atomic_sizes,
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

    # Plotting the basics
    detector_screen_position = [0.7, 0.2, 0.25, 0.625]
    fig = plt.figure(figsize=(10, 4))
    ax = fig.gca(projection="3d")
    ax.set_position([0, 0, 0.7, 1])

    # Create second set of axes for detection screen
    ax2 = plt.axes(detector_screen_position)
    ax2.tick_params(axis="both", labelbottom=False, labelleft=False)

    # Plot atoms
    ax.scatter(atomic_positions[:, 0], atomic_positions[:, 1],
               atomic_positions[:, 2], c=atomic_colors, s=atomic_sizes)

    for line in pruned_lines:
        ax.plot(line[0], line[1], line[2],
                color=g_col, linewidth=g_w, alpha=g_a)

    # Plotting the beam: First we create the beam display vector
    ax.quiver(0, 0, beam_end_z, k_disp[0], k_disp[1], k_disp[2],
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
                        ax.plot_surface(p[0], p[1], p[2], color="r",
                                        shade=False, alpha=0.2)
                    # We also plot the outgoing line corresponding to this
                    # scattering. First we get the point (and squeeze it, to
                    # make it 1D again)
                    p = np.squeeze(points[indices_index, :])
                    start = np.array([0, 0, beam_end_z])
                    ray = p - start
                    line = np.array([start, start + outgoing_length * ray])
                    ax.plot(line[:, 0], line[:, 1], line[:, 2], color='r',
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
                    ax.quiver(starts[:, 0],
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
        ax2.scatter(points[:, 0], points[:, 1], c=colors)
        for i in range(len(indices)):
            x, y = points[i, 0:2] - 0.05 * ranges
            s = indices[i]
            c = colors[i, :-1]
            ax2.text(x, y, s, color=c, va='top', ha='right')

        # Plotting detection plane
        abs_ex = 0.0
        rel_ex = 0.1
        def_y = 2
        def_x = 2
        x_min = np.amin(points[:, 0]) * (1 + rel_ex) - abs_ex
        x_max = np.amax(points[:, 0]) * (1 + rel_ex) + abs_ex
        y_min = np.amin(points[:, 1]) * (1 + rel_ex) - abs_ex
        y_max = np.amax(points[:, 1]) * (1 + rel_ex) + abs_ex
        x_range = np.array([min(x_min, -def_x), max(x_max, def_x)])
        y_range = np.array([min(y_min, -def_y), max(y_max, def_y)])
        x, y = np.meshgrid(x_range, y_range)
        z = plane_z * np.ones(x.shape)
        ax.plot_surface(x, y, z, color='k', alpha=0.2)

        # plotting intersections
        ax.scatter(points[:, 0], points[:, 1], plane_z, color=colors)

        # Setting limits for the second figure
        det_max_x = np.amax(x_range)
        det_min_x = np.amin(x_range)
        det_max_y = np.amax(y_range)
        det_min_y = np.amin(y_range)
        det_max = max(det_max_x, det_max_y)
        det_min = min(det_min_x, det_min_y)
        ax2.set_xlim(det_min, det_max)
        ax2.set_ylim(det_min, det_max)

        if show_all:
            # Plotting outgoing vectors
            n = k_out.shape[0]
            k_plot = k_out / lattices.mag(k_in)
            start_point = np.array((0, 0, beam_end_z))
            start_points = np.repeat(np.atleast_2d(start_point), n, axis=0)
            ax.quiver(start_points[:, 0],
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
                ax.plot(line[:, 0], line[:, 1], line[:, 2],
                        color='k', alpha=0.3, ls='--',
                        lw=g_w)

    ax.set_aspect('equal')
    ax.set_proj_type('ortho')
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))

    # Some limit trickery. We make the plot box cubic:
    plot_max = np.amax(r_max)
    plot_min = np.amin(r_min)
    ax.set_xlim(plot_min, plot_max)
    ax.set_ylim(plot_min, plot_max)
    ax.set_zlim(plot_min, plot_max)

    ax.grid(False)
    ax.axis('off')

    tit = (r'Scattering on a cubic lattice. $k_{in} = (2\pi/a)\cdot$' +
           '{}'.format(k_title))
    tit2 = (r'Scattering on a cubic lattice. $k_{in} = $' +
            '{}'.format(k_title))
    if normalize:
        ax.set_title(tit)
    else:
        ax.set_title(tit2)
    ax2.set_title('Detection screen.\nForm factors: {}'.format(form_factor))
    if returns:
        return fig, ax, ax2
    plt.show()


def Band_structure(V0=0, n_k=51, G_range=list(range(-3, 4)),
                   potential="harmonic", edges=False,
                   E_F=None,
                   returns=False):

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

    if returns:
        return fig, ax, ax2
    plt.show()
