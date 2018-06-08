import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lattices
import scattering
import band_structure

eq = np.isclose

d = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
     np.array([0, 0, 0]), "xkcd:cement", 2, "proper", "latticevectors",
     [2, 2, 2])


def Lattice(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, max_=d[8], lattice_name=None,
        type_=None, indices=None, verbose=False, returns=False):
    """
    Creates, limits and plots the lattice
    """

    # Minimum coefficients for lattice vectors
    min_ = [0, 0, 0]

    num_plane_points = 20
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

    # Choosing gridline and unit cell type. First the default settings.
    latticelines = lattices.latticelines
    unitcells = lattices.unitcells

    if grid_type is None:
        grid_type = latticelines[lattice_type]

    if type_ is None:
        type_ = unitcells[lattice_type]
    else:
        try:
            type_ = type_.lower()
        except AttributeError:
            print('Please input a string for type_. Giving primitive')
            type_ = "primitive"

        if type_ not in ["primitive", "conventional"]:
            print(("Input either 'primitive' or 'conventional' for type."
                   " Giving 'primitive'"))
            type_ = "primitive"

    # set the range of lattice vectors to be calculated
    r_min, r_max, n_min, n_max = lattices.find_limits(lim_type, a1, a2, a3,
                                                      min_, max_, type_=type_)
    
    # if we plot the conventional cell we want to give r_min and r_max to
    # limiter. If not we want to give n_min and n_max
    if type_ == "conventional":
        lim_min, lim_max = r_min, r_max
    else:
        lim_min, lim_max = n_min, n_max

    objects = lattices.generator(a1, a2, a3, basis, colors, sizes,
                                 n_min, n_max)
    # Objects to limit to the plot-box
    objects = lattices.limiter(objects[0], objects, lim_min, lim_max)
    (atomic_positions, lattice_coefficients, atomic_colors, atomic_sizes,
     lattice_position) = objects

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
    ax.axis('off')

    # make the panes transparent (the plot box)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    
    if returns:
        return fig, ax
    
    plt.show()


def Reciprocal(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        lim_type=d[6], grid_type=None, max_=d[8], lattice_name=None, type_=None,
        indices=(1, 1, 1), verbose=False, returns=False):
    if returns:
        fig, ax = Lattice(a1=a1, a2=a2, a3=a3,
                          basis=basis,
                          colors=colors,
                          sizes=sizes,
                          lim_type=lim_type,
                          grid_type=grid_type,
                          max_=max_,
                          type_=type_,
                          lattice_name=lattice_name,
                          indices=indices,
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
                type_=type_,
                lattice_name=lattice_name,
                indices=indices,
                verbose=verbose)


def Scattering(lattice_name='simple cubic',
               basis=None,
               k_in=np.array([0, 0, -1.5]),
               scattering_length=np.array([1, 1, 1, 1]),
               highlight=None,
               show_all=False,
               normalize=True,
               verbose=False,
               returns=False):

    k_in = np.array(k_in)
    point_sizes = 2
    lattice_name = lattice_name.lower()
    min_, max_ = (-2, -2, -1), (2, 2, 1)
    grid_type = lattices.latticelines[lattice_name]
    lim_type = "proper"
    atom_colors = ["xkcd:cement",
                   "xkcd:cornflower blue",
                   "xkcd:cornflower blue",
                   "xkcd:cornflower blue"]
    atom_sizes = [1, 1, 1, 1]
    g_col = 'k'
    g_w = 0.5
    g_a = 0.6
    size_default = 36
    point_sizes *= size_default
    plane_z = 2
    beam_end_z = max_[2]
    unit_cell_type = "conventional"
    # input sanitization
    if basis is not None:
        a1, a2, a3 = np.eye(3, dtype=int)
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
    
    # Normalizing wave vector (multiplying by k0 = 2Pi/a)
    k_title = np.copy(k_in)
    if normalize:
        k_in *= 2 * np.pi
    
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
                                                             scattering_length,
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
              color='b', lw=2, pivot='tip', length=lambda_)

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
                    d, planes = lattices.reciprocal(a1, a2, a3, hi_index,
                                                    r_min - extra,
                                                    r_max + extra,
                                                    points=20)
                    planes = lattices.plane_limiter(planes, r_min - extra,
                                                    r_max + extra)
                    high_intensity = intensities[indices_index]
                    colors[indices_index] = [1, 0, 0, high_intensity]
                    for p in planes:
                        ax.plot_surface(p[0], p[1], p[2], color="r",
                                        shade=False, alpha=0.2)

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
        ax.scatter(points[:, 0], points[:, 1], plane_z)

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
                line_x = [start_point[0], p[0]]
                line_y = [start_point[1], p[1]]
                line_z = [start_point[2], p[2]]
                ax.plot(line_x, line_y, line_z, color='k', alpha=0.3, ls='--',
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
    ax.set_title(r'Scattering on a cubic lattice. $k_{in} = (2\pi/a)\cdot$' + '{}'.format(k_title))
    ax2.set_title('Detection screen.\nScattering lengths: {}'.format(scattering_length))
    if returns:
        return fig, ax, ax2
    plt.show()


def Band_structure(V0=0, n_k=101, G_range=[-1, 0, 1], potential="harmonic"):

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

    o = band_structure.calc_band_structure(V0=V0, n_k=n_k, G_range=G_range,
                                           potential=potential)
    kxs, kys, band, max_E = o
    max_E_mat = max_E * np.ones((n_k, n_k))
    max_k = np.amax(kxs)
    min_k = np.amin(kxs)

    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.plot_surface(kxs, kys, band)
    ax.plot_surface(kxs, kys, max_E_mat, alpha=0.5)
    ax.set_xlim([min_k, max_k])
    ax.set_ylim([min_k, max_k])
    plt.show()
