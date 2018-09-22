# In this file you should create a function called "lattice_position", which
# takes in the following arguments:
# - The primitive lattice vectors: a1, a2, a3
# - The coefficients for the current lattice point: nx, ny, nz
# - The array containing the positions of the basis atoms: basis
#
# The function must return a numpy array, containing the position of the atoms
# for this lattice point. The size of the numpy array should be the same size
# as the "basis"-array

import numpy as np


def lattice_position(a1, a2, a3, nx, ny, nz, basis):
    position = nx * a1 + ny * a2 + nz * a3
    atomic_positions = []
    for atom in basis:
        atomic_positions.append(position + atom)
    atomic_positions = np.array(atomic_positions)
    return position, atomic_positions
