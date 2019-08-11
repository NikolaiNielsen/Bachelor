import numpy as np


def calc_positions(a1, a2, a3, n1, n2, n3, basis):
    """
    Finish writing this function.
    The function must calculate the lattice point (R) and the absolute position
    of each atom in the basis, for this particular lattice point.

    Inputs:
    - a1, a2, a3: 3-component vector. The primitive lattice vectors
    - n1, n2, n3: integers. The lattice point coefficients.
    - basis: (n, 3) array. The positions of the atoms in the basis. Examples:
             "basis[0, 2]" refers to the first atoms z-component,
             "basis[2, 1]" refers to the third atoms y-component,
             "basis[1, 0]" refers to the second atoms x-component.

    Required outputs:
    - lattice_point: 3-component vector. The lattice vector for this specific
                     lattice point
    - atomic_positions: (n, 3) array. The absolute positions of the atoms in
                        the basis, for this specific lattice point.

    Sample inputs:
    a1 = [1, 0, 0],
    a2 = [0, 1, 0],
    a3 = [0, 0, 1],
    n1 = 2,
    n2 = 1,
    n3 = 3,
    basis = [[  0,   0,   0],
             [0.5, 0.5, 0.5]]

    Sample outputs:
    lattice_point = [2, 1, 3]
    atomic_positions = [[  2,   1,   3],
                        [2.5, 1.5, 3.5]]
    """

    # Mock up:
    # Calculate the current lattice point
    lattice_point = np.array((0, 0, 0))

    # Get the number of atoms in the basis
    n = basis.shape[0]

    # Calculate the atomic positions
    atomic_positions = np.zeros((n, 3))

    # Sample solution

    lattice_point = n1*a1 + n2*a2 + n3*a3
    # Reshape the lattice_point, so we can easily add it to the basis
    lattice_add = lattice_point.reshape(1, 3)
    atomic_positions = basis + lattice_add

    return lattice_point, atomic_positions


def calc_reciprocal_lattice(a1, a2, a3, m1, m2, m3):
    """
    Finish this function.
    The function must calculate the primitive lattice vectors for the
    reciprocal lattice, along with the specific reciprocal lattice point given
    by the coefficients.

    Inputs:
    - a1, a2, a3: 3-component vectors. Primitive lattice vectors of the direct
                  lattice.
    - m1, m2, m3: integers. Coefficients for the reciprocal lattice point.

    Outputs:
    - b1, b2, b3: 3-component vectors. Primitive lattice vectors of the
                  reciprocal lattice.
    - G: 3-component vector. The reciprocal lattice point.

    Sample inputs:
    a1 = [1, 0, 0]
    a2 = [0, 1, 0]
    a3 = [0, 0, 1]
    m1 = 3
    m2 = 1
    m3 = 0

    Sample outputs:
    b1 = [2*pi, 0, 0]
    b2 = [0, 2*pi, 0]
    b3 = [0, 0, 2*pi]
    G = [6*pi, 2*pi, 0]
    """

    # Mock up:
    b1 = np.array((0, 0, 0))
    b2 = np.array((0, 0, 0))
    b3 = np.array((0, 0, 0))
    G = np.array((0, 0, 0))

    # Sample Solution

    # Calculate the scale factor
    scale = (2 * np.pi) / np.dot(a1, np.cross(a2, a3))
    # Calculate the primitive lattice vectors
    b1 = scale * np.cross(a2, a3)
    b2 = scale * np.cross(a3, a1)
    b3 = scale * np.cross(a1, a2)
    # Calculate the lattice point
    G = m1*b1 + m2*b2 + m3*b3

    return b1, b2, b3, G
