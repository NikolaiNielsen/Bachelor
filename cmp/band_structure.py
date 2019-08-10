import itertools
import numpy as np
import lattices


def VG_dirac(coeff, V0=1):
    return V0


def VG_cos(coeff, V0=1):
    allowed = np.array([[0, 1], [0, -1], [-1, 0], [1, 0]])
    if (allowed == coeff).all(axis=1).any():
        VG = V0 / 2
    else:
        VG = 0
    return VG


def potential_matrix(range_=[-1, 0, 1], potential=VG_dirac, V0=1,
                     coeff_matrix=False):
    sets = list(itertools.product(range_, range_,))
    coefficient_array = np.array(sets)

    num_el = len(range_)**2
    rows = np.arange(num_el)
    cols = np.arange(num_el)

    V_mat = np.zeros((num_el, num_el))

    for row in rows:
        psi = coefficient_array[row]
        for col in cols:
            V = coefficient_array[col]
            test = psi - V
            V_mat[row, col] = potential(test, V0)

    if coeff_matrix:
        V_mat_char = np.chararray((num_el, num_el), unicode=True, itemsize=7)
        for row in rows:
            psi = coefficient_array[row]
            for col in cols:
                V = coefficient_array[col]
                test = psi - V
                V_mat_char[row, col] = "{}".format(test)
        return V_mat, V_mat_char

    return V_mat


def calc_band_structure(V0=100, n_k=101, G_range=[-1, 0, 1],
                        potential=VG_dirac):
    kx = np.linspace(-1 / 2, 1 / 2, n_k)
    ky = np.linspace(-1 / 2, 1 / 2, n_k)
    num_Gs = (len(G_range))**2
    # First we create the relevant matrix for some k:
    b1 = np.array([1, 0])
    b2 = np.array([0, 1])

    ms = np.array(list(itertools.product(G_range, G_range)))
    recip = np.array([b1, b2])
    Gs = ms @ recip
    E = np.zeros((num_Gs, n_k, n_k))
    VG_mat = potential_matrix(range_=G_range, potential=potential, V0=V0)
    kxs, kys = np.meshgrid(kx, ky)
    for i in range(n_k):
        for j in range(n_k):
            k = np.array([kx[i], ky[j]])
            Diag = np.diag(lattices.mag(k - Gs)**2) / 2
            Full = Diag + VG_mat
            Eigs = np.linalg.eigvalsh(Full)
            E[:, i, j] = Eigs

    band_to_return = E[0]

    flattened = band_to_return.flatten()
    sorted_ = np.sort(flattened)
    index = np.arange(sorted_.size)
    highest = np.ceil((n_k**2) / 2).astype('int')
    allowed = sorted_[index < highest]
    max_E = np.amax(allowed)

    return kxs, kys, band_to_return, max_E
