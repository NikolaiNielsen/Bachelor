# This is the top-level file for the project. This file contains the main
# functions for each of the subjects dealt with in the thesis.
#
# For now there are three functions:
# - Lattice: This function takes in a given set of lattice vectors and a basis
#   and constructs the corresponding crystal structure.
# - Reciprocal: This function does the same as Lattice, except it will also
#   plot a given family of lattice planes along with the lattice
# - Scattering: This function will plot a given crystal structure (simple
#   cubic, FCC or BCC are allowed), and given an incoming wave vector will
#   calculate the corresponding scattering with said crystal. These scattered
#   wave vectors will then be "detected" on a screen, simulating a neutron-
#   scattering experiment.
#
# Usage of the functions are as follows.
# --------------------
# Lattice
# --------------------
# The desired crystal structure can be specified in two ways:
# - The name of the Bravais Lattice Type can be specified as the lattice_name
#   argument
# - The lattice vectors (a1, a2, a3) can be specified as 3 ndarrays, containing
#   3 elements each. The basis can be specified as an ndarray, containing n
#   rows and 3 columns
#
# There is also a list of optional arguments that can be passed to the Lattice
# function:
# - colors: a list/tuple with n strings specifying the colors of each of the
#   atoms in the basis. A single string can also be given, if only one color is
#   desired
# - sizes: a list/tuple of n numbers specifying the relative size of the atoms.
#   2 is the default value.
# - grid_type: either "latticevectors" or "axes". "latticevectors" plots
#   gridlines along the lattice vectors. "axes" tries to plot gridlines
#   parallel to the x, y and z-axes (requires a crystal structure that can be
#   specified as a orthorhombic lattice with basis)
# - min_, max_: 2 list/tuples specifying the minimum/maximum amount of lattice
#   points to plot in each direction of the lattice vectors.
# - lim_type: provides 3 different ways of limiting the plotted crystal
#   structure: "individual", "sum" and "proper". The default is "proper", which
#   calculates the 8 lattice points given by linear combinations of the
#   lattice-vectors and the coefficients specified in min_ and max_. This forms
#   a parallelipiped, where a orthorhombic plot-box is calculated to fit these
#   8 vertices inside.
# - verbose: True/False boolean value. Makes the functions print a lot of
#   information helpful for debugging. Only partially implemented
# - indices: a list/tuple of 3 elements, corresponding to the Miller indices h,
#   j and l. If specified this plots the corresponding (family of) lattice
#   planes. For the Lattice function, this argument defaults to None, meaning
#   no planes are plotted
#
# --------------------
# Reciprocal
# --------------------
# This function is actually just a copy of the Lattice function, with a default
# argument for indices of (1, 1, 1). Otherwise everything else is identical
#
# --------------------
# Scattering
# --------------------
# 
#
#
#
