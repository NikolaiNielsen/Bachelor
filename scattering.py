# Lattice plotting
import itertools
import numpy as np
import lattices
eq = np.isclose


def scattering(a1, a2, a3, basis, scattering_length, k_in):
    basis = np.atleast_2d(basis)
    length_basis = np.shape(basis)
    if len(length_basis) == 1:
        n_basis = 1
    elif len(length_basis) > 1:
        n_basis = length_basis[0]

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
    l = range(l_min, l_max + 1)

    # Create the array of indices and exclude the (0,0,0)
    indices = np.asarray(list(itertools.product(h, j, l)))
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


def smart_cubic_k(indices, k0=None, k1=None, k2=None):
    indices = np.array(indices)
    scaler = np.pi * indices.dot(indices)

    # Some input sanitization
    if indices.shape != (3,):
        print("We need 3 indices and only 3 indices!")
        return

    # More sanitization
    ks = np.array([k0 is None, k1 is None, k2 is None])
    if np.sum(ks) != 1:
        print("You need to specify 2 and only 2 of the k-coordinates!")
        return

    # Now that we've made sure only 1 component is unspecified, only one of
    # these if statements will execute
    if ks[0]:
        k0 = - (scaler + indices[1] * k1 + indices[2] * k2) / indices[0]
    if ks[1]:
        k1 = - (scaler + indices[0] * k0 + indices[2] * k2) / indices[1]
    if ks[2]:
        k2 = - (scaler + indices[0] * k0 + indices[1] * k1) / indices[2]

    return np.array([k0, k1, k2])


def projection(k_array, r=np.array([0, 0, 0]), n=np.array([0, 0, 1]), p0=np.
               array([0, 0, 5])):

    # We assume k_array is a 2D array, with shape (N,3) where N is the number
    # of scattered wave vectors
    num = n.dot(p0 - r)
    den = k_array @ n
    d = num / den
    # Next we make sure that there are enough d's: make it 2d:
    D = np.vstack((d, d, d)).T
    p = r + k_array * D
    return p
