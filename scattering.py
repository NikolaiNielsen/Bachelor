# Calculations for scattering
import itertools
import numpy as np
import lattices

eq = np.isclose


def calc_scattering(a1, a2, a3, basis, scattering_length, k_in):
    """
    Calculates the scattering off a given lattice, for a given set of
    scattering lengths and incident wavevector
    """

    basis = np.atleast_2d(basis)
    n_basis = np.shape(basis)[0]

    # First the scaling factor for the reciprocal lattice
    scale = a1.dot(np.cross(a2, a3))
    # Then the reciprocal lattice
    b1 = 2 * np.pi * np.cross(a2, a3) / scale
    b2 = 2 * np.pi * np.cross(a3, a1) / scale
    b3 = 2 * np.pi * np.cross(a1, a2) / scale
    lattice = np.array([b1, b2, b3])

    # Then we need to create all the reciprocal lattice vectors
    h_min, j_min, l_min = [-5] * 3
    h_max, j_max, l_max = [5] * 3

    # We subtract 1, because we exclude the (0,0,0) index
    num_G = ((h_max + 1 - h_min) * (j_max + 1 - j_min) *
             (l_max + 1 - l_min)) - 1
    G_array = np.zeros((num_G, 3))

    # Create the range of miller indices
    h = range(h_min, h_max + 1)
    j = range(j_min, j_max + 1)
    ell = range(l_min, l_max + 1)

    # Create the array of indices and exclude the (0,0,0)
    indices = np.asarray(list(itertools.product(h, j, ell)))
    not_zeros = np.sum(indices == 0, 1) != 3
    indices = indices[not_zeros]

    # Create the array of reciprocal lattice vectors
    G_array = indices @ lattice

    # Next we calculate the structure factors (we assume neutron scattering)
    structure_factors = np.zeros(num_G, dtype=complex)
    counter = 0
    for G in G_array:
        for n in range(n_basis):
            atom = basis[n, ]
            structure_factors[counter] += (scattering_length[n] *
                                           np.exp(1j * G.dot(atom)))
        counter += 1
    intensity = structure_factors * np.conj(structure_factors)
    intensity = intensity.real.astype(float)
    # Then we get the scattered wavevectors:
    k_out = G_array + k_in

    # And select only those that actually satisfy |k_in| = |k_out|
    mag_equal = eq(lattices.mag(k_out), lattices.mag(k_in))

    # delete all indices, intensities and k_out, that don't satisfy the above
    indices = indices[mag_equal]
    intensity = intensity[mag_equal]
    k_out = k_out[mag_equal]
    return intensity, k_out, indices


def smart_cubic_k(indices, k0=None, k1=None, k2=None, theta=None, phi=None):
    """
    Creates a "smart" wave vector, which is guaranteed to scatter. It needs a
    desired family of lattice planes and at least two components (in cartesian
    coordinates) or two angles (in spherical polar).
    """

    indices = np.array(indices)
    scaler = np.pi * indices.dot(indices)

    # Some input sanitization
    if indices.shape != (3,):
        print("We need 3 indices and only 3 indices!")
        return

    # Next we check which of the components are specifies. We prefer angles
    angles = np.array([theta is not None, phi is not None])
    ks = np.array([k0 is None, k1 is None, k2 is None])
    if angles.all():
        # We check if the values are acceptable
        if 0 <= theta <= np.pi / 2 and 0 <= phi < 2 * np.pi:
            r = -scaler / (np.sin(theta) * (indices[0] * np.cos(phi) +
                                            indices[1] * np.sin(phi)) +
                           indices[2] * np.cos(theta))
            if r < 0:
                r = -r
            # Also we flip the coordinates, as it otherwise points upwards
            k0 = -r * np.sin(theta) * np.cos(phi)
            k1 = -r * np.sin(theta) * np.sin(phi)
            k2 = -r * np.cos(theta)

        else:
            print("You need 0 <= theta < pi/2 and 0 <= phi <= 2*pi")
    elif np.sum(ks) == 1:
        # Only one unspecified component, so only one if-statement will execute
        if ks[0]:
            k0 = - (scaler + indices[1] * k1 + indices[2] * k2) / indices[0]
        if ks[1]:
            k1 = - (scaler + indices[0] * k0 + indices[2] * k2) / indices[1]
        if ks[2]:
            k2 = - (scaler + indices[0] * k0 + indices[1] * k1) / indices[2]

        if k2 > 0:
            # we flip the sign on k2 to make sure it's negative.
            k2 = -k2
    else:
        print("You need to specify only 2 of the k-coordinates, or 2 angles")
        k0, k1, k2 = 0, 0, np.pi

    return np.array([k0, k1, k2])


def projection(k_array, r=np.array([0, 0, 1]), n=np.array([0, 0, 1]), p0=np.
               array([0, 0, 5])):
    """
    Calculates the projections of the scattered vector onto the detector plate.
    """
    # We assume k_array is a 2D array, with shape (N,3) where N is the number
    # of scattered wave vectors
    num = n.dot(p0 - r)
    den = k_array @ n
    d = num / den
    # Next we make sure that there are enough d's: make it 2d:
    D = np.vstack((d, d, d)).T
    p = r + k_array * D
    return p


def prune_scattered_points(points, intensities, indices):
    """
    We prune the array of points, so we're left with unique points. The
    intensities are adjusted accordingly, and the indices are made into a list
    of strings, with one element per unique point
    """
    # First we get the unique points, the index of the new points from the old
    # array, and a list of the inverse ids (list of n elements, where n is the
    # number of points. Each element is the id of the corresponding unique
    # point in the new array)
    unique_points, ids, inverse = np.unique(points, axis=0, return_index=True,
                                            return_inverse=True)

    # next we create an array of zeros for the new intensities
    new_intensities = np.zeros(ids.shape)

    # And a empty list with a guarantied number of elements
    new_indices = [None] * new_intensities.size

    for i in range(intensities.size):
        # The id of the new, unique point
        new_id = inverse[i]
        # A tuple of the old index row
        index = tuple(indices[i])
        # We add the old intensity from a given point, to the new unique points
        # intensity
        new_intensities[new_id] += intensities[i]

        # We populate the list of new indices, taking care of whether or not an
        # element of the list is empty
        if new_indices[new_id] is None:
            new_indices[new_id] = "{}".format(index)
        else:
            new_indices[new_id] = "{}\n{}".format(new_indices[new_id], index)

    return unique_points, new_intensities, new_indices


def projection_sphere(k_array, r=5, o=np.array([0, 0, 1])):
    """
    calculates the (positive) intersection between the line defined by o and
    k_array, and the sphere with radius r and origin o
    """
    d = r / lattices.mag(k_array)
    D = np.vstack((d, d, d)).T
    p = o + k_array * D
    return p
