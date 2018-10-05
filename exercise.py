# In this file you should create a function called "lattice_position", which
# takes in the following arguments:
# - The primitive lattice vectors: a1, a2, a3
# - The coefficients for the current lattice point: nx, ny, nz
# - The array containing the positions of the basis atoms: basis
#
# The function must return a numpy array, containing the position of the atoms
# for this lattice point. The size of the numpy array should be the same size
# as the "basis"-array
#
# Sample input:
# a1 = [1, 0, 0], a2 = [0, 1, 0], a3 = [0, 0, 1]
# nx = 2, ny = 1, nz = 3
# basis = [[  0,   0,   0],
#          [0.5, 0.5, 0.5]]
#
# Sample output:
# atomic_positions = [[  2,   1,   3],
#                     [2.5, 1.5, 3.5]]

import numpy as np


def lattice_position(a1, a2, a3, nx, ny, nz, basis):
    position = nx * a1 + ny * a2 + nz * a3
    atomic_positions = []
    for atom in basis:
        atomic_positions.append(position + atom)
    atomic_positions = np.array(atomic_positions)
    return position, atomic_positions
