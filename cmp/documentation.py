# cmp.py is the main file for the project. cmp.py contains the main functions
# for each of the subjects dealt with in the thesis. This file
# (documentation.py) contains documentation on how to use these functions and
# provides annotated examples. This file can be run as is, and the examples
# given will be executed
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
#   specified as a orthorhombic lattice with basis). Defaults to the preffered
#   type of gridlines for a given Bravais lattice type
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
# This function simulates neutron scattering on a simple cubic lattice with a
# basis. It plots two figures in one window. One holds the crystal, a
# representation of the incoming neutron beam, the detection screen and the
# points as detected by this screen. The second figure shows just the detection
# plane and the points detected by it, along with the Miller indices giving
# rise to this particular scattering event. The alpha value of the points on
# the screen corresponds to the normalized intensity of the scattering events.

# The function needs 3 pieces of information to work properly. The basis,
# incoming wave vector and scattering lengths for the atoms in the basis.

# The basis can be specified in the following ways:
# - lattice_name: a string, specifying the Bravais Lattice Type (currently it
#   only accepts "simple cubic", "bcc" and "fcc")
# - basis: an ndarray, with n rows and 3 columns, specifying the absolute
#   position of the atoms within the unit cell (assuming lattice spacing of 1)

# The two other requires arguments are
# - k_in: an ndarray with 3 numbers. The script simulate a neutron beam
#   incident on the top of the material. As such, the z-component should be
#   negative. If positive, the script will flip the sign.
# - scattering_length: an ndarray with n numbers, where n is the number of
#   atoms in the basis. Due to normalization of the intensity, only the
#   relative values of elements in scattering_length matter. As such it is
#   recommended to input values near unity to minimize any rounding error. By
#   default all scattering lengths are equal

# Other optional arguments that can be passed are as follows:
# - highlight: a list/tuple of 3 elements. Corresponds to the Miller indices h,
#   j and l. This will highlight the given scattering event associated with
#   these Miller indices (if present)
# - show_all: True/False boolean. If set to true this will plot all outgoing
#   wave vectors (with lengths equal to their wavelength), and lines going from
#   the scattering point out to the points on the detection plane.
# - verbose: True/False boolean. Prints information useful in debugging. Also
#   only partially implemented
#
# --------------------
# SETUP
# --------------------
# Before using these scripts, the numpy and matplotlib packages should be
# installed in your python environment. I suggest installing the "anaconda"
# python stack (https://anaconda.com), which installs python and pretty much
# any package you need for any type of scientific computing (including
# matplotlib and numpy, of course).
#
# When using these scripts, make sure that the files are in the python path. I
# recommend doing this by using this folder as your working directory (ie,
# launch the python shell from this folder, or create a python script in this
# folder which includes the commands you want to execute)
#
# If you are using jupyter notebook, create one in this folder, and include a
# line with "%matplotlib notebook" (without quotes). This will make the plot
# windows appear in the notebook, interactively, and not as separate windows
#
# --------------------
# EXAMPLE USE
# --------------------
# For all of these examples to work, the following line needs to be present at
# the top of the document:
from cmp import *
# The following line plots an FCC lattice (primitive unit cell) using the
# lattice_name argument:
Lattice(lattice_name="fcc")
# An FCC lattice, specifying the lattice vectors and basis manually
# (conventinal unit cell):
Lattice(a1=np.array([1, 0, 0]),
        a2=np.array([0, 1, 0]),
        a3=np.array([0, 0, 1]),
        basis=np.array([[0, 0, 0],
                        [0.5, 0.5, 0],
                        [0.5, 0, 0.5],
                        [0, 0.5, 0.5]]))
# A simple cubic lattice with the (1, 1, 1) planes plotted
Reciprocal(lattice_name="simple cubic",
           indices=(1, 1, 1))
# Scattering on a BCC lattice, with k_in at normal incidence, (magnitude of
# 2*pi), highlighting the (0, 0, 2) planes and plotting all information. All
# scattering lengths equal (which is the default)
Scattering(lattice_name="bcc",
           k_in=np.array([0, 0, -2 * np.pi]),
           highlight=(0, 0, 2),
           show_all=True)

# --------------------
# GLOSSARY
# --------------------
# "np" is the numpy package. The standard package for doing numerical work.
# A "list" is a comma-separated list of elements, enclosed by square brackets:
a = [1, 2, 3]
# A "tuple" is like a list, but enclosed by parentheses:
b = (1, 2, 3)
# An "ndarray" is the multidimensional array type specified by numpy. This
# corresponds to the default MATLAB arrays. It is invoked by passing a
# list/tuple of list/tuples to the np.array() function. If I want to create a
# unit vector in the x-direction I type
x = np.array([1, 0, 0])
# If I want to specify the basis for a BCC lattice (conventional unit cell) I
# need to pass a list of lists. Each sub-list is a row in the resulting array.
# Note the double square brackets: ([[], []])
basis = np.array([[0, 0, 0],
                  [0.5, 0.5, 0.5]])
