"""
Next:
Reciprocal lattice and lattice planes
Scattering
- Parameters:
  - Incident angle
  - Incident wavelength of light
  - Start with simple cubic, and 90 degree incident angle (specific case)
Distortion of Fermi Surface:
- Transform Hamiltonian to k-space. Transforms from a matrix to just a number
"""

import itertools

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import lattice

eq = np.isclose

d = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
     np.array([0, 0, 0]), "xkcd:cement", 2, "proper", "latticevectors",
     [0, 0, 0], [2, 2, 2])


def Lattice(
        a1=d[0], a2=d[1], a3=d[2], basis=d[3], colors=d[4], sizes=d[5],
        LimType=d[6], GridType=None, Mins=d[8], Maxs=d[9], Lattice=None,
        verbose=False):

    lattice.creator(a1, a2, a3, basis, colors, sizes, LimType,
                    GridType, Mins, Maxs, Lattice, verbose)
