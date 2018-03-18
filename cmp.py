# Let's get the importing out of the way.
# Numpy for calculations, matplotlib for plotting.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import itertools

# alias for isclose
eq = np.isclose
# Lattice plotting
# Defaults for LatticeCreator
d = (np.array([1, 0, 0]), np.array([0, 1, 0]), np.array([0, 0, 1]),
     np.array([0, 0, 0]), "xkcd:cement", 2, "wdynamic", "latticevectors",
     [0, 0, 0], [2, 2, 2])


def LatticeCreator(a1=d[0], a2=d[1], a3=d[2],
                   basis=d[3], colors=d[4], sizes=d[5],
                   LimType=d[6], GridType=d[7], Mins=d[8], Maxs=d[9],
                   Lattice=None):
    """
    Creates and limits the lattice
    """
    if Lattice is not None:
        lattice, basis = LatticeChooser(Lattice)
        a1, a2, a3 = lattice

    size_default = 36
    # Input sanitization:
    # We need the number of basis-vectors.
    # If there is only 1 basis vector, then len(np.shape(basis)) == 1
    # otherwise the length is 2, and the first element is number of basis
    # vectors
    length_basis = np.shape(basis)
    if len(length_basis) == 1:
        N_basis = 1
    elif len(length_basis) > 1:
        N_basis = length_basis[0]

    # Make a list, N_basis long, for the colors and sizes,
    # if they're not specified.
    c_name = colors.__class__.__name__
    if c_name == "str":
        c = colors
        colors = []
        for i in range(N_basis):
            colors.append(c)
    elif c_name == "list" and len(colors) < N_basis:
        c = colors[0]
        colors = []
        for i in range(N_basis):
            colors.append(c)

    s_name = sizes.__class__.__name__
    if s_name == "int" or s_name == "float":
        s = sizes
        sizes = []
        for i in range(N_basis):
            sizes.append(s)
    elif s_name == "list" and len(sizes) < N_basis:
        s = sizes[0]
        sizes = []
        for i in range(N_basis):
            sizes.append(s)
    # set the range of lattice vectors to be calculated
    r_min, r_max, n_min, n_max = FindLimits(LimType, a1, a2, a3, Mins, Maxs)
    # Calculate the amount of atomic positions to be calculated
    numAtoms = ((n_max[0] + 1 - n_min[0]) * (n_max[1] + 1 - n_min[1]) *
                (n_max[2] + 1 - n_min[2]) * N_basis)

    # Make a zero array for all of the atomic positions. numAtoms in one
    # direction and 3 in the other (coordinates)
    AtomicPositions = np.zeros((numAtoms, 3))
    # Empty lists for colors, sizes and whether or not they're lattice points
    AtomicColors = []
    AtomicSizes = []
    LatticePosition = []

    # Loop over all chosen linear combinations of basis vectors and plot each
    counter = 0
    for nx in range(n_min[0], n_max[0] + 1):
        for ny in range(n_min[1], n_max[1] + 1):
            for nz in range(n_min[2], n_max[2] + 1):
                lattice_position = nx * a1 + ny * a2 + nz * a3
                for n_atom in range(N_basis):
                    AtomicPositions[counter, ] = (lattice_position +
                                                  basis[n_atom, ])
                    AtomicColors.append(colors[n_atom])
                    AtomicSizes.append(size_default * sizes[n_atom])
                    if (AtomicPositions[counter, ] == lattice_position).all():
                        LatticePosition.append(True)
                    else:
                        LatticePosition.append(False)
                    counter += 1

    # Another way to do this is to use itertools.product to create all
    # permutations of -2, ..., 4 with repeat of 3, and then use np.asarray() to
    # convert this into a numpy array. The "problem" is that this doesn't allow
    # one to have nx_max = / = ny_max, etc. All ranges must be equal.
    # I should check to see which is fastest.
    # Strike that above problem. Just pass it a list for each coordinate with
    # the range and use no repeat.
    # AtomicCoefficients = np.asarray(list(itertools.product(x, y, z)))
    # Where x, y, z is list of integers from nx_min to nx_max etc.
    # This would yield list of coefficients (nx, ny, nz), then we just multiply
    # the first dimension by a1, the second by a2 and so on. But not now

    atoms, dims = np.shape(AtomicPositions)
    # Get the rows with the function above
    rows = Limiter(AtomicPositions, r_min, r_max)
    # delete all rows (axis 0 of the array) that are outside the limits
    AtomicPositions = np.delete(AtomicPositions, rows, 0)
    # Go through the list of rows to delete in reverse order, and delete what's
    # needed from colors and sizes
    for ID in sorted(rows, reverse=True):
        del AtomicColors[ID]
        del AtomicSizes[ID]
        del LatticePosition[ID]
    LatticePlotter(a1, a2, a3, AtomicPositions, AtomicColors, AtomicSizes,
                   LatticePosition, GridType, r_min, r_max)


# Let's first define some things
def mag(a):
    """
    Returns magnitude of vector or each row of an array
    """
    # Return magnitude of vector
    if len(a.shape) == 1:
        return np.linalg.norm(a)
    # Return magnitude of each row of an array.
    else:
        return np.linalg.norm(a, axis=1)


def LatticeClassifier(a1, a2, a3, basis):
    """
    test all bravais lattice types (primitive unit cells for all, conventional
    for fcc and bcc). It works by first checking how many of the lattice
    vectors have an equal magnitude, and then checking the angles between the
    lattice vectors. The angles are checked below with the extensive amount of
    boolean values.
    """

    # Create a lattice array and get the magnitude of the lattice vectors
    lattice = np.array([a1, a2, a3])
    mag_lattice = mag(lattice)
    mag_a1, mag_a2, mag_a3 = mag_lattice

    # Angles between lattice vectors
    cos12 = a1.dot(a2) / (mag_a1 * mag_a2)
    cos31 = a1.dot(a3) / (mag_a1 * mag_a3)
    cos23 = a2.dot(a3) / (mag_a2 * mag_a3)

    # bool of equality between lattice vector magnitudes
    mag_12Eq = eq(mag_a1, mag_a2)
    mag_31Eq = eq(mag_a1, mag_a3)
    mag_23Eq = eq(mag_a2, mag_a3)

    # If all magnitudes are equal, then there's equality between all magnitudes
    mag_AllEq = mag_12Eq and mag_31Eq and mag_23Eq

    # We check each of the permutations where only 1 pair is equal to each
    # other in magnitude
    mag_Only12Eq = (mag_12Eq and (not mag_31Eq) and (not mag_23Eq))
    mag_Only31Eq = (mag_31Eq and (not mag_12Eq) and (not mag_23Eq))
    mag_Only23Eq = (mag_23Eq and (not mag_31Eq) and (not mag_12Eq))
    # XOR the above permutations together to make sure only one is true
    mag_Only2Eq = mag_Only12Eq ^ mag_Only31Eq ^ mag_Only23Eq

    # Check for orthogonality
    ortho12 = eq(0, np.dot(a1, a2))
    ortho31 = eq(0, np.dot(a1, a3))
    ortho23 = eq(0, np.dot(a2, a3))

    # all are orthogonal
    ortho = ortho12 and ortho31 and ortho23

    # only two pairs are orthogonal
    ortho2 = ((ortho12 and ortho31) ^ (ortho12 and ortho23) ^
              (ortho31 and ortho23))

    # The three possible permutations for hexagonal lattice. Hexagonal can both
    # have 3 equal sides or only two
    hexa3 = eq(0.5, cos12) and eq(0, cos31) and eq(0, cos23)
    hexa2 = eq(0.5, cos31) and eq(0, cos23) and eq(0, cos12)
    hexa1 = eq(0.5, cos23) and eq(0, cos12) and eq(0, cos31)
    hexa = hexa1 ^ hexa2 ^ hexa3

    fcc = eq(0.5, cos12) and eq(0.5, cos31) and eq(0.5, cos23)

    # rhombohedral have all angles equal to each other, but is not fcc
    rhombo = (eq(cos12, cos31) and eq(cos12, cos23) and
              eq(cos31, cos23) and not fcc)

    # the three bcc permutations available
    bcc3 = (eq(0, cos12) and eq(np.sqrt(3) / 3, cos31) and
            eq(np.sqrt(3) / 3, cos23))
    bcc2 = (eq(0, cos31) and eq(np.sqrt(3) / 3, cos23) and
            eq(np.sqrt(3) / 3, cos12))
    bcc1 = (eq(0, cos23) and eq(np.sqrt(3) / 3, cos12) and
            eq(np.sqrt(3) / 3, cos31))
    bcc = bcc1 ^ bcc2 ^ bcc3

    # The three tetragonal body centred permutations
    tbc3 = eq(0, cos12) and eq(cos23, cos31) and eq(cos23, mag_a2 /
                                                    (2 * mag_a3))
    tbc2 = eq(0, cos31) and eq(cos12, cos23) and eq(cos12, mag_a1 /
                                                    (2 * mag_a2))
    tbc1 = eq(0, cos23) and eq(cos31, cos12) and eq(cos31, mag_a3 /
                                                    (2 * mag_a1))
    tbc = tbc1 ^ tbc2 ^ tbc3

    # The three tetragonal face centred permutations
    tfc1 = (eq(cos12, cos31) and eq(cos12, mag_a1 / (2 * mag_a2)) and
            eq(cos23, (2 * mag_a2**2 - mag_a1**2) / (2 * mag_a3**2)))
    tfc2 = (eq(cos31, cos23) and eq(cos31, mag_a3 / (2 * mag_a1)) and
            eq(cos12, (2 * mag_a1**2 - mag_a3**2) / (2 * mag_a2**2)))
    tfc3 = (eq(cos23, cos12) and eq(cos23, mag_a2 / (2 * mag_a3)) and
            eq(cos31, (2 * mag_a3**2 - mag_a2**2) / (2 * mag_a1**2)))
    tfc = tfc1 ^ tfc2 ^ tfc3

    # Tetragonal base centred
    tBase3 = eq(cos12, np.sqrt(2) / 2) and eq(cos31, cos23) and eq(0, cos23)
    tBase2 = eq(cos31, np.sqrt(2) / 2) and eq(cos23, cos12) and eq(0, cos12)
    tBase1 = eq(cos23, np.sqrt(2) / 2) and eq(cos12, cos31) and eq(0, cos31)
    tBase = tBase1 ^ tBase2 ^ tBase3

    # Base centred monoclinic has 6 different permutations, and can have either
    # no sides equal, two sides equal or all sides equal. With two or three
    # sides equal it has a 2D triangular lattice, where each 2D lattice is
    # displaced with respect to the other.
    BaseMono3 = (eq(cos12, mag_a1 / (2 * mag_a2)) and
                 eq(cos23, a1.dot(a3) / (2 * mag_a2 * mag_a3)))
    BaseMono2 = (eq(cos31, mag_a3 / (2 * mag_a1)) and
                 eq(cos12, a3.dot(a2) / (2 * mag_a1 * mag_a2)))
    BaseMono1 = (eq(cos23, mag_a2 / (2 * mag_a3)) and
                 eq(cos31, a2.dot(a1) / (2 * mag_a3 * mag_a1)))
    BaseMono4 = (eq(cos31, mag_a1 / (2 * mag_a3)) and
                 eq(cos23, a1.dot(a2) / (2 * mag_a3 * mag_a2)))
    BaseMono5 = (eq(cos23, mag_a3 / (2 * mag_a2)) and
                 eq(cos12, a3.dot(a1) / (2 * mag_a2 * mag_a1)))
    BaseMono6 = (eq(cos12, mag_a2 / (2 * mag_a1)) and
                 eq(cos31, a2.dot(a3) / (2 * mag_a1 * mag_a3)))
    BaseMono = (BaseMono1 ^ BaseMono2 ^ BaseMono3 ^
                BaseMono4 ^ BaseMono5 ^ BaseMono6)

    # Orthorhombic body centred
    obc1 = (eq(cos12, 0) and eq(cos31, mag_a1 / (2 * mag_a3)) and
            eq(cos23, mag_a2 / (2 * mag_a3)))
    obc2 = (eq(cos31, 0) and eq(cos23, mag_a3 / (2 * mag_a2)) and
            eq(cos12, mag_a1 / (2 * mag_a2)))
    obc3 = (eq(cos23, 0) and eq(cos12, mag_a2 / (2 * mag_a1)) and
            eq(cos31, mag_a3 / (2 * mag_a1)))
    obc = obc1 ^ obc2 ^ obc3

    # Just need one here. Permutations lead to the exact same expression
    ofc = (eq(cos12, (mag_a1**2 + mag_a2**2 - mag_a3**2) /
                     (2 * mag_a1 * mag_a2)) and
           eq(cos23, (-mag_a1**2 + mag_a2**2 + mag_a3**2) /
                     (2 * mag_a2 * mag_a3)) and
           eq(cos31, (mag_a1**2 - mag_a2**2 + mag_a3**2) /
                     (2 * mag_a3 * mag_a1)))

    # oBase False positives BaseMono, since the dot product in the
    # angle formulas all give 0
    # Again we have 6 possible permutations
    oBase1 = eq(cos12, mag_a1 / (2 * mag_a2)) and ortho23 and ortho31
    oBase2 = eq(cos31, mag_a3 / (2 * mag_a1)) and ortho12 and ortho23
    oBase3 = eq(cos23, mag_a2 / (2 * mag_a3)) and ortho31 and ortho12
    oBase4 = eq(cos23, mag_a3 / (2 * mag_a2)) and ortho12 and ortho31
    oBase5 = eq(cos31, mag_a1 / (2 * mag_a3)) and ortho23 and ortho12
    oBase6 = eq(cos12, mag_a2 / (2 * mag_a1)) and ortho31 and ortho23
    oBase = oBase1 ^ oBase2 ^ oBase3 ^ oBase4 ^ oBase5 ^ oBase6

    # Triclinic has no angles equal
    tri = ((not eq(cos12, cos23)) and (not eq(cos23, cos31)) and
           (not eq(cos31, cos12)))

    # We start with an undetermined lattice type
    LatticeType = "undetermined"
    if mag_AllEq:
        # Side lengths are equal. Lattice types:
        # Cubic,
        # rhombohedral,
        # face centered cubic,
        # hexagonal with a=|a_1|,
        # Tetragonal Body centered with b=+-sqrt(2)a
        # base centred monoclinic, b = +-sqrt(3)a, c = a
        # conventional fcc
        # conventional bcc
        if ortho:
            # Let's detect the conventional unit cells of fcc and bcc. Requires
            # more than one basis vector
            if len(basis.shape) == 1:
                LatticeType = "simple cubic"
            else:
                # we exclude the first basis vector ([0,0,0])
                reduced_basis = basis[1:]
                mag_basis = mag(reduced_basis)
                # bcc has one basis vector in the reduced basis
                if reduced_basis.shape[0] == 1:
                    # Make sure the reduced basis is a vector
                    reduced_basis = reduced_basis.flatten()
                    # check if the vector has the right length
                    lengthEq = eq(np.sqrt(3) / 2, mag_basis / mag_a1)
                    # calculate angles
                    angles = lattice.dot(reduced_basis) / (mag_basis * mag_a1)
                    # make sure the angles are all the right magnitude
                    anglesEq = np.all(eq(angles, np.sqrt(3) / 3))
                    # if lengths and angles are correct, it's bcc
                    if lengthEq and anglesEq:
                        LatticeType = "conventional bcc"

                # fcc has 3 basis vectors in the reduced basis
                elif reduced_basis.shape[0] == 3:
                    # check if all length ratios are sqrt(2)/2
                    lengthEq = eq(np.sqrt(2) / 2, mag_basis / mag_a1).all()

                    # Calculate angles between lattice vectors and basis
                    # vectors
                    normalizer = np.outer(mag_lattice, mag_basis)
                    angles = lattice.dot(reduced_basis) / normalizer
                    # Calculate rank of matrix
                    rank = np.linalg.matrix_rank(angles)

                    # get angles that are sqrt(2)/2 and 0 respectively
                    num_sqrt2 = eq(angles, np.sqrt(2) / 2)
                    num_0 = eq(angles, 0)

                    # Check whether there are 2 or 1 of these per row
                    # respectively
                    num_sqrt2_true = (num_sqrt2.sum(axis=1) ==
                                      np.array([2, 2, 2])).all()
                    num_0_true = (num_0.sum(axis=1) ==
                                  np.array([1, 1, 1])).all()
                    anglesEq = rank == 3 and num_0_true and num_sqrt2_true
                    if lengthEq and anglesEq:
                        LatticeType = "conventional fcc"

        elif hexa:
            LatticeType = "hexagonal 1"
        elif fcc:
            LatticeType = "fcc"
        elif tbc:
            LatticeType = "tetragonal body centred"
        elif BaseMono:
            LatticeType = "base centred monoclinic 1"
        elif rhombo:
            LatticeType = "rhombohedral"
        else:
            pass

    elif mag_Only2Eq:
        # Only two lengths are equal. Possible lattices
        # BCC
        # Tetragonal bc
        # tetragonal fc
        # tetragonal cubic
        # Cubic Base Centered
        # Base centered monoclinic (b=+-sqrt(3) * a)
        # Hexagonal
        # simple monoclinic, a=b or a=c or b=c
        if bcc:
            LatticeType = "bcc"
        # tbc actually gives a false positive for regular bcc.
        elif tbc:
            LatticeType = "tetragonal body centred"
        elif tfc:
            LatticeType = "tetragonal face centred"
        elif tBase:
            LatticeType = "base centred cubic"
        elif ortho:
            LatticeType = "tetragonal"
        elif hexa:
            LatticeType = "hexagonal 2"
        elif BaseMono:
            LatticeType = "base centred monoclinic 2"
        elif ortho2:
            LatticeType = "simple monoclinic"
        else:
            pass
    else:
        # no side lengths are equal. Possible lattices:
        # Orthorhombic
        # OBC
        # OFC
        # Orthorhombic base centered
        # Tetragonal base centered (primitive), b != a and b != sqrt(2)a
        # Simple Monoclinic
        # Base centered Monoclinic
        # Triclinic
        if ortho:
            LatticeType = "orthorhombic"
        elif obc:
            LatticeType = "orthorhombic body centred"
        elif ofc:
            LatticeType = "orthorhombic face centred"
        elif tBase:
            LatticeType = "tetragonal base centred"
        elif oBase:
            LatticeType = "orthorhombic base centred"
        elif BaseMono:
            LatticeType = "base centred monoclinic 3"
        elif ortho2:
            LatticeType = "simple monoclinic"
        elif tri:
            LatticeType = "triclinic"
        else:
            pass
    return LatticeType


def RotMatrix(v=np.array([1, 1, 1]), theta=np.pi / 4):
    """
    Generates the rotation matrix for rotation about a given vector with a
    given angle. See https://en.wikipedia.org/wiki/Rotation_matrix
    """
    # Make sure we have a unit vector
    v = v / mag(v)
    # Create the cross product matrix
    vCross = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    # Tensor product
    vTens = np.tensordot(v, v, 0)
    # Return rotation matrix
    return (np.cos(theta) * np.identity(3) + np.sin(theta) *
            vCross + (1 - np.cos(theta)) * vTens)


def LatticeChooser(lattice_name="simple cubic"):
    """
    Outputs the chosen lattice and basis
    """
    # Let's just sanitize the input
    lattice_name = lattice_name.lower()
    L = {}
    B = {}
    a, b, c, theta = 1, 1.5, 2, 80 * np.pi / 180
    # Create the relevant lattices (transposed - using row vectors)
    # Simple cubic
    lcubic = np.array([[a, 0, 0], [0, a, 0], [0, 0, a]])
    L["simple cubic"] = lcubic
    # BCC
    lbcc = np.array([[a, 0, 0], [0, a, 0], [a / 2, a / 2, a / 2]])
    L["bcc"] = lbcc
    # FCC
    lfcc = np.array([[a / 2, a / 2, 0], [a / 2, 0, a / 2], [0, a / 2, a / 2]])
    L["fcc"] = lfcc
    # Base Centered Cubic
    lcbase = np.array([[a, 0, 0], [a / 2, a / 2, 0], [0, 0, a]])
    L["base centred cubic"] = lcbase
    # Tetragonal
    ltetra = np.array([[a, 0, 0], [0, a, 0], [0, 0, b]])
    L["tetragonal"] = ltetra
    # Tetragonal Body Centred
    ltbc = np.array([[a, 0, 0], [0, a, 0], [a / 2, a / 2, b / 2]])
    L["tetragonal body centred"] = ltbc
    # Tetragonal Face Centred
    ltfc = np.array([[a / 2, a / 2, 0], [a / 2, 0, b / 2], [0, a / 2, b / 2]])
    L["tetragonal face centred"] = ltfc
    # tetragonal base centred
    ltbase = np.array([[a / 2, a / 2, 0], [0, a, 0], [0, 0, b]])
    L["tetragonal base centred"] = ltbase
    # Orthorhombic
    lortho = np.array([[a, 0, 0], [0, b, 0], [0, 0, c]])
    L["orthorhombic"] = lortho
    # Orthorhombic body centred
    lobc = np.array([[a, 0, 0], [0, b, 0], [a / 2, b / 2, c / 2]])
    L["orthorhombic body centred"] = lobc
    # Orthorhombic Face centred
    lofc = np.array([[a / 2, b / 2, 0], [a / 2, 0, c / 2], [0, b / 2, c / 2]])
    L["orthorhombic face centred"] = lofc
    # Orthorhombic base centred
    lobase = np.array([[a, 0, 0], [a / 2, b / 2, 0], [0, 0, c]])
    L["orthorhombic base centred"] = lobase
    # simple monoclic
    lsmono = np.array([[a, 0, 0], [0, b, 0],
                       [c * np.cos(theta), 0, c * np.sin(theta)]])
    L["simple monoclinic"] = lsmono
    # base centred monoclinic
    lbcmono = np.array([[a, 0, 0], [a / 2, b / 2, 0],
                        [c * np.cos(theta), 0, c * np.sin(theta)]])
    L["base centred monoclinic 3"] = lbcmono
    # Base centred monoclinic (2)
    lbcmono2 = np.array([[a, 0, 0], [a / 2, np.sqrt(3) * a / 2, 0],
                         [c * np.cos(theta), 0, c * np.sin(theta)]])
    L["base centred monoclinic 2"] = lbcmono2
    # Base centred monoclinic (3)
    lbcmono3 = np.array([[a, 0, 0], [a / 2, np.sqrt(3) * a / 2, 0],
                         [a * np.cos(theta), 0, a * np.sin(theta)]])
    L["base centred monoclinic 1"] = lbcmono3
    # Hexagonal 1
    lhexa1 = np.array([[a, 0, 0], [a / 2, np.sqrt(3) * a / 2, 0], [0, 0, a]])
    L["hexagonal 1"] = lhexa1
    # Hexagonal 2
    lhexa2 = np.array([[a, 0, 0], [a / 2, np.sqrt(3) * a / 2, 0], [0, 0, b]])
    L["hexagonal 2"] = lhexa2
    # Triclinc stuff
    gamma = 70 * np.pi / 180
    beta = 60 * np.pi / 180
    cx = c * np.cos(beta)
    cy = c * (np.cos(theta) - np.cos(beta) * np.cos(gamma)) / np.sin(gamma)
    cz = np.sqrt(c**2 - cx**2 - cy**2)
    ltri = np.array([[a, 0, 0], [b * np.cos(gamma), b * np.sin(gamma), 0],
                     [cx, cy, cz]])
    L["triclinic"] = ltri
    # Rhombohedral
    lrhombo = np.array([[a, b, b], [b, a, b], [b, b, a]])
    L["rhombohedral"] = lrhombo

    # conventional fcc
    L["conventional fcc"] = lcubic
    B["conventional fcc"] = a * np.array([[0.5, 0.5, 0],
                                          [0.5, 0, 0.5],
                                          [0, 0.5, 0.5]])

    L["conventional bcc"] = lcubic
    B["conventional bcc"] = a * np.array([0.5, 0.5, 0.5])

    try:
        lattice = L[lattice_name]
    except KeyError:
        print("You did da dumdum, and I now give you simple cubic")
        lattice = L["simple cubic"]

    basisOrigin = np.array([0, 0, 0])

    try:
        basis = B[lattice_name]
    except KeyError:
        basis = np.array([])

    if basis.shape[0] > 0:
        basis = np.vstack((basisOrigin, basis))
    elif basis.shape[0] == 0:
        basis = np.hstack((basisOrigin, basis))
    else:
        print("something went horribly wrong")

    return lattice, basis


def LatticeTester(verbose=False):
    """
    Tests all the lattices for lattice detection, with permutations and
    rotations
    """
    # List of all the lattices
    lattices = ["simple cubic", "fcc", "bcc", "conventional fcc",
                "conventional bcc", "base centred cubic", "tetragonal",
                "tetragonal body centred", "tetragonal face centred",
                "tetragonal base centred", "orthorhombic",
                "orthorhombic base centred", "orthorhombic body centred",
                "orthorhombic face centred", "simple monoclinic",
                "base centred monoclinic 1", "base centred monoclinic 2",
                "base centred monoclinic 3", "hexagonal 1", "hexagonal 2",
                "triclinic", "rhombohedral"]

    # Create the rotation matrix
    R = RotMatrix()
    for name in lattices:
        # Create the lattice
        lattice, basis = LatticeChooser(name)
        # rotate the lattice and basis
        lattice = (R@lattice.T).T
        basis = (R@basis.T).T
        for perm in itertools.permutations([0, 1, 2]):
            # permute the lattice
            a1, a2, a3 = lattice[list(perm)]

            # next we classify it
            LatticeType = LatticeClassifier(a1, a2, a3, basis)

            if verbose:
                print("Lattice: {}. Classification: {}. Permutation {}".format(
                      name,
                      LatticeType,
                      perm))
            else:
                if name != LatticeType:
                    s = "L: {}, C: {}, P: {}".format(name, LatticeType, perm)
                    print(s)
    if verbose:
        print("Test done.")
    else:
        print("Test done. If nothing printed, all were succesfully classified")


def FindLimits(LimType, a1, a2, a3, Min=[0, 0, 0], Max=[2, 2, 2]):
    """
    Calculates the limits on the coordinates (the plot box), and the limits on
    the basis vector ranges.
    """

    # If we have hard limit types we just pass the min and max as coordinate
    # limits and calculate limits of the basis vector ranges based on
    # coordinate limits
    # THIS METHOD IS UGLY BECAUSE WE DIVIDE BY 0 AND GET NAN AND INF.
    # I NEED TO CLEAN IT UP
    if LimType.lower() in "hard":
        r_min, r_max = np.array(Min), np.array(Max)
        lattice = np.array((a1, a2, a3))
        # Ideally we want n1_max * a1 = r_max. But they're not (generally)
        # parallel. So what we do is we dissolve into components, so
        # r_max_x / a1_x = nx_max and similarly for the other components
        # But we want nx_max to be an integer, so we take the ceiling value
        # (round up) to allow for "spillage" and then we take the maximal
        # value of the three coordinate ratios (again, for extra spillage).
        # We do the same procedure (but with floor and minimum) for the n_min
        # Let's get rid of those pesky infinities and nans
        # (from x / 0 and 0 / 0 respectively). Just replace by 0
        max_quot = r_max / lattice
        min_quot = r_min / lattice
        max_quot[np.isnan(max_quot) + np.isinf(max_quot)] = 0
        min_quot[np.isnan(min_quot) + np.isinf(min_quot)] = 0
        n_max = np.amax(np.ceil(max_quot), 0)
        n_min = np.amin(np.floor(min_quot), 0)

    # For dynamic limits we pass Min and Max as limits of basis vector range
    # and calculate coordinate limit based on basis vector range
    elif LimType.lower() in "dynamic":
        n_min, n_max = np.array(Min), np.array(Max)
        lattice = np.array((a1, a2, a3))
        # Take the max value for each of the cardinal directions, for the
        # three scaled lattice vectors (so x_max is max x value of
        # Max[0] * a1, Max[1] * a2 and Max[2] * a3).
        # this can be done by multiplying the transposed lattice matrix by the
        # n_max vector, then taking max value
        max_vects = lattice.T * n_max
        r_max = np.amax(max_vects, 0)
        # Similar for minimums:
        min_vects = lattice.T * n_min
        r_min = np.amin(min_vects, 0)
    # Different type of coordinate limits. Take r_max as sum(lattice * max)
    # Works well for orthogonal or near-orthogonal lattice vectors
    elif LimType.lower() in "wdynamic":
        n_min, n_max = np.array(Min), np.array(Max)
        lattice = np.array([a1, a2, a3])
        r_max = np.sum(lattice.T * n_max, 0)
        r_min = np.sum(lattice.T * n_min, 0)
    else:
        print("You did poop...")
    # And lastly we return the relevant arrays, with n_min / max -+ some value
    # to allow for "spillage". The value is the maximal value of the Max array.
    # Also, let's make sure n_min / max are arrays of integers. Don't worry,
    # they've already been rounded
    return (r_min, r_max, n_min.astype('int') - np.max(Max),
            n_max.astype('int') + np.max(Max))


def Limiter(l, r_min=np.array([0, 0, 0]), r_max=np.array([2, 2, 2])):
    """
    A function to highlight points that are outside the limits of the plot
    """
    rows = []
    num, _ = np.shape(l)
    # loop over all row ID's
    for rowID in range(num):
        # if the atom is outside the limits, we append the row ID to a list
        row = l[rowID, ]
        # (it's actually easier and prettier to check if they're inside)
        inside_x = r_min[0] <= row[0] <= r_max[0]
        inside_y = r_min[1] <= row[1] <= r_max[1]
        inside_z = r_min[2] <= row[2] <= r_max[2]
        inside = inside_x and inside_y and inside_z
        if not inside:
            rows.append(rowID)
    return rows


def rotatefcc(a1, a2, a3, basis):
    """
    Rotate the (primitive) fcc lattice for easier gridline plotting
    """
    a1prop = np.array([1 / 2, 1 / 2, 0])
    a2prop = np.array([1 / 2, 0, 1 / 2])
    a3prop = np.array([0, 1 / 2, 1 / 2])

    # First we orient a1 along (0.5,0.5,0)
    a1cross = np.cross(a1, a1prop)
    theta = np.arcsin(mag(a1cross) / (mag(a1) * mag(a1prop)))
    r1 = RotMatrix(a1cross, theta)
    a1, a2, a3, basis = rotate(a1, a2, a3, basis, r1)

    # Next we try orienting a2 along a2prop. Then we see if a3 is aligned
    # properly. We need to rotate along a1, such that this is kept fixed

    # To find the angle we need to rotate it with, we use the vector rejection
    # of a2 and a2prop, as these are perpendicular to a1
    theta, r2 = RotMatrixAlong(a1prop, a2, a2prop)
    A1, A2, A3, Basis = rotate(a1, a2, a3, basis, r2)
    if eq(A2, a2prop).all():
        # We rotated properly
        a1, a2, a3, basis = A1, A2, A3, Basis
    else:
        # Rotate the other way!
        r2 = RotMatrix(a1prop, -theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r2)

    # Now we check if a3 is right. If not we rotate a2 along a3prop instead
    if eq(a3, a3prop).all():
        pass
    else:
        theta, r2 = RotMatrixAlong(a1prop, a2, a3prop)
        A1, A2, A3, Basis = rotate(a1, a2, a3, basis, r2)
        if eq(A2, a3prop).all():
            # We rotated properly
            a1, a2, a3, basis = A1, A2, A3, Basis
        else:
            # Rotate the other way!
            r2 = RotMatrix(a1prop, -theta)
            a1, a2, a3, basis = rotate(a1, a2, a3, basis, r2)

    return a1, a2, a3, basis


def RotMatrixAlong(a, b, c):
    """
    creates the rotation matrix which rotates b about a, such that it coincides
    with c
    """

    # First we need the relevant vector rejections
    brej = b - b.dot(a) / (mag(a)**2) * a
    crej = c - c.dot(a) / (mag(a)**2) * a

    # Next we get the angle between the rejections
    theta = np.arccos(brej.dot(crej) / (mag(brej) * mag(crej)))

    # and the relevant rotation matrix
    R = RotMatrix(a, theta)
    return theta, R


def rotateHex(a1, a2, a3, basis):
    """
    Rotator for the hexagonal structure
    """
    mag_a1 = mag(a1)
    mag_a2 = mag(a2)
    mag_a3 = mag(a3)
    cos12 = a1.dot(a2) / (mag_a1 * mag_a2)
    cos31 = a1.dot(a3) / (mag_a1 * mag_a3)
    cos23 = a2.dot(a3) / (mag_a2 * mag_a3)
    x = np.array([1, 0, 0])
    z = np.array([0, 0, 1])

    # Rotate the lattice according to which vectors form the triangular lattice
    if eq(0.5, cos12):
        # a1 and a2 form triangular lattice. Align a1 along x
        a1cross = np.cross(a1, x)
        theta = np.arcsin(mag(a1cross) / mag_a1)
        r1 = RotMatrix(a1cross, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r1)

        # now rotate a3 so it's parallel to z, with x as rotation axis
        theta = np.arccos(a3.dot(z) / (mag_a3))
        r3 = RotMatrix(x, theta)
        A1, A2, A3, Basis = rotate(a1, a2, a3, basis, r3)
        # Alright, so I don't know how to get the proper rotation direction, so
        # we just check if it's rotated properly
        if eq(A2[2], 0):
            # We rotated correctly
            a1, a2, a3, basis = A1, A2, A3, Basis
        else:
            # We didn't rotate correctly, so we rotate the other way twice
            r3 = RotMatrix(x, -theta)
            a1, a2, a3, basis = rotate(a1, a2, a3, basis, r3)

    elif eq(0.5, cos23):
        # a2 and a3 form triangular lattice. Align a2 along x
        a2cross = np.cross(a2, x)
        theta = np.arcsin(mag(a2cross) / mag_a2)
        r2 = RotMatrix(a2cross, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r2)

        # now rotate a1 so it's parallel to z, with x as rotation axis
        theta = np.arccos(a1.dot(z) / (mag_a1))
        r1 = RotMatrix(x, theta)
        A1, A2, A3, Basis = rotate(a1, a2, a3, basis, r1)
        # Alright, so I don't know how to get the proper rotation direction, so
        # we just check if it's rotated properly It's rotated properly if the
        # z-coordinate of A3 is 0
        if eq(A3[2], 0):
            # We rotated correctly
            a1, a2, a3, basis = A1, A2, A3, Basis
        else:
            # We didn't rotate correctly, so we rotate the other way twice
            r1 = RotMatrix(x, -theta)
            a1, a2, a3, basis = rotate(a1, a2, a3, basis, r1)

    elif eq(0.5, cos31):
        # a1 and a3 form triangular lattice. Align a1 along x
        a1cross = np.cross(a1, x)
        theta = np.arcsin(mag(a1cross) / mag_a1)
        r1 = RotMatrix(a1cross, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r1)

        # now rotate a2 so it's parallel to z, with x as rotation axis
        theta = np.arccos(a2.dot(z) / (mag_a2))
        r2 = RotMatrix(x, theta)
        A1, A2, A3, Basis = rotate(a1, a2, a3, basis, r2)
        # Alright, so I don't know how to get the proper rotation direction, so
        # we just check if it's rotated properly
        if eq(A3[2], 0):
            # We rotated correctly
            a1, a2, a3, basis = A1, A2, A3, Basis
        else:
            # We didn't rotate correctly, so we rotate the other way twice
            r3 = RotMatrix(x, -theta)
            a1, a2, a3, basis = rotate(a1, a2, a3, basis, r3)
    else:
        print('something went wrong')

    return a1, a2, a3, basis


def rotate(a1, a2, a3, basis, R):
    """
    Rotates the lattice
    """
    return R@a1, R@a2, R@a3, (R@basis.T).T


def rotator(a1, a2, a3, basis, latticetype=None, verb=False):
    """
    Rotates the lattice to make plotting gridlines easier
    """
    # We remember, that |a x b| = |a| |b| sin(theta)
    eq = np.isclose
    x = np.array([1, 0, 0])
    y = np.array([0, 1, 0])
    # Check for orthogonality
    ortho12 = eq(0, np.dot(a1, a2))
    ortho31 = eq(0, np.dot(a1, a3))
    ortho23 = eq(0, np.dot(a2, a3))

    if verb:
        print("Before:")
        print(a1)
        print(a2)
        print(a3)
        print(basis)

    if "hexagonal" in latticetype:
        a1, a2, a3, basis = rotateHex(a1, a2, a3, basis)
    elif latticetype == "fcc":
        a1, a2, a3, basis = rotatefcc(a1, a2, a3, basis)
    elif ortho12:
        # We choose a1 to align along x
        a1cross = np.cross(a1, x)
        theta = np.arcsin(mag(a1cross) / mag(a1))
        r1 = RotMatrix(a1cross, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r1)

        # Now we align a2 along y
        # But we gotta make sure we rotate in the right direction
        a2cross = np.cross(a2, y)
        sign = np.sign(a2cross)[0]
        theta = sign * np.arcsin(mag(a2cross) / mag(a2))
        r2 = RotMatrix(x, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r2)
    elif ortho23:
        # We choose a2 to align along x
        a2cross = np.cross(a2, x)
        theta = np.arcsin(mag(a2cross) / mag(a2))
        r2 = RotMatrix(a2cross, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r2)

        # Now we align a3 along y
        a3cross = np.cross(a3, y)
        sign = np.sign(a3cross)[0]
        theta = sign * np.arcsin(mag(a3cross) / mag(a3))
        r3 = RotMatrix(x, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r3)
    elif ortho31:
        # We choose a1 to align along x
        a1cross = np.cross(a1, x)
        theta = np.arcsin(mag(a1cross) / mag(a1))
        r1 = RotMatrix(a1cross, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r1)

        # Now we align a2 along y
        a3cross = np.cross(a3, y)
        sign = np.sign(a3cross)[0]
        theta = sign * np.arcsin(mag(a3cross) / mag(a3))
        r3 = RotMatrix(x, theta)
        a1, a2, a3, basis = rotate(a1, a2, a3, basis, r3)
    else:
        # Well, it doesn't really matter here, if none of them are orthogonal
        # to each other. We'll just use latticevector gridlines and leave this
        # be.
        pass

    # Let's sanitize the coordinates that are close to 0 (say if a1x =
    # 2*10^-10, then we set it equal 0)
    a1[eq(a1, 0)] = 0
    a2[eq(a2, 0)] = 0
    a3[eq(a3, 0)] = 0
    basis[eq(basis, 0)] = 0
    
    if verb:
        print("after")
        print(a1)
        print(a2)
        print(a3)
        print(basis)

    return a1, a2, a3, basis


def CreateLines(points, v1, v2, v3,
                r_min=np.array([0, 0, 0]),
                r_max=np.array([2, 2, 2])):
    """
    Creates lines along vectors and limits these to the given plot box
    """
    # Let's try and create some grid lines along the lattice vectors some
    # observations regarding grid lines: if we go along a lattice vector to a
    # point, then we need gridlines that are along the direction of the other
    # two lattice vectors A naïve approach would be to just plot 3 lines for
    # each atomic position (after limiting). This would create multiple copies
    # of some gridlines, but it's an easy solution. Let's see how long it takes
    # to compute!

    # Create the plotting parameter for the gridlines.
    # They use the equation r = r0+t * a, where:
    # r0 is a fixed point (atomic position)
    # a is a vector giving the direction of the gridline (lattice vector)
    # t is a scaling parameter, creating the points along the line.
    # We use halfinteger steps for t. That way we know that we'll properly hit
    # other atomic positions. If we used linspace this wouldn't be the case

    t = np.arange(-10, 10, 0.5)

    lines = []

    # Create all gridlines needed and append them to the lines-list
    numPoints = np.shape(points)[0]
    for rowID in range(numPoints):
        CurrentPoint = points[rowID, ]
        line1 = CurrentPoint + np.outer(t, v1)
        line2 = CurrentPoint + np.outer(t, v2)
        line3 = CurrentPoint + np.outer(t, v3)
        lines.append([line1, line2, line3])

    # run through each line and clip points outside limits
    pruned_lines = []
    for point in lines:
        for line in point:
            # Get the points outside the plot and delete them
            rows = Limiter(line, r_min, r_max)
            line = np.delete(line, rows, 0)

            # Because we're working with arrays we're passing copies, we need
            # to append the pruned lines to a new list. And let's only add
            # them if there are actually any points to plot
            line_length, _ = np.shape(line)
            if line_length > 0:
                pruned_lines.append(line)
    return pruned_lines


def LatticePlotter(a1, a2, a3, AtomicPositions, AtomicColors, AtomicSizes,
                   LatticePosition, GridType, r_min, r_max):
    """
    Takes the input lattice, adds gridlines and plots everything
    """
    # Create the figure
    fig = plt.figure()
    ax = fig.gca(projection="3d")

    # Plot atoms. For now a single size and color
    ax.scatter(AtomicPositions[:, 0], AtomicPositions[:, 1],
               AtomicPositions[:, 2], c=AtomicColors, s=AtomicSizes)

    # Create grid lines
    g_col = 'k'
    g_w = 0.5
    rangex = range(int(np.ceil(r_min[0])), int(np.floor(r_max[0])) + 1)
    rangey = range(int(np.ceil(r_min[1])), int(np.floor(r_max[1])) + 1)
    rangez = range(int(np.ceil(r_min[2])), int(np.floor(r_max[2])) + 1)
    lowGrid = GridType.lower()
    if lowGrid in "hard":
        for nx in rangex:
            for ny in rangey:
                ax.plot(np.array([nx, nx]), np.array([ny, ny]),
                        np.array([np.ceil(r_min[2]), np.floor(r_max[2])]),
                        c=g_col, linewidth=g_w)

            for nz in rangez:
                ax.plot(np.array([nx, nx]),
                        np.array([np.ceil(r_min[1]), np.floor(r_max[1])]),
                        np.array([nz, nz]),
                        c=g_col, linewidth=g_w)

        for ny in rangey:
            for nz in rangez:
                ax.plot(np.array([np.ceil(r_min[0]), np.floor(r_max[0])]),
                        np.array([ny, ny]), np.array([nz, nz]),
                        c=g_col, linewidth=g_w)

    elif lowGrid in "latticevectors":
        # gridlines along lattice vectors - really messy for non-orthogonal
        # latticevectors
        pruned_lines = CreateLines(AtomicPositions[LatticePosition],
                                   a1, a2, a3, r_min, r_max)
        for line in pruned_lines:
            ax.plot(line[:, 0], line[:, 1], line[:, 2], c=g_col, linewidth=g_w)

    elif lowGrid in "soft":
        # A Way of finding atoms on cartesian axes
        # bool array of atoms with x = 0 and y = 0
        x0 = AtomicPositions[:, 0] == 0
        y0 = AtomicPositions[:, 1] == 0
        z0 = AtomicPositions[:, 2] == 0

        # Get Lattice spacings
        # z-values of atoms on the z-axis
        z_vals = AtomicPositions[x0 * y0, 2]
        # Keep those with z > 0
        z_vals = z_vals[z_vals > 0]
        # Take the minimum as the lattice spacing
        a_z = np.min(z_vals)

        y_vals = AtomicPositions[x0 * z0, 1]
        y_vals = y_vals[y_vals > 0]
        a_y = np.min(y_vals)

        x_vals = AtomicPositions[y0 * z0, 0]
        x_vals = x_vals[x_vals > 0]
        a_x = np.min(x_vals)

        for nx in np.arange(r_min[0], r_max[0] + 1, a_x):
            for ny in np.arange(r_min[1], r_max[1] + 1, a_y):
                ax.plot(np.array([nx, nx]), np.array([ny, ny]),
                        np.array([r_min[2], r_max[2]]), c=g_col, linewidth=g_w)

            for nz in np.arange(r_min[2], r_max[2] + 1, a_z):
                ax.plot(np.array([nx, nx]), np.array([r_min[1], r_max[1]]),
                        np.array([nz, nz]), c=g_col, linewidth=g_w)

        for ny in np.arange(r_min[1], r_max[1] + 1, a_y):
            for nz in np.arange(r_min[2], r_max[2] + 1, a_z):
                ax.plot(np.array([r_min[0], r_max[0]]), np.array([ny, ny]),
                        np.array([nz, nz]), c=g_col, linewidth=g_w)

    else:
        print("No Gridlines Chosen")

    # plot lattice vectors
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
    ax.set_xlim3d([r_min[0], r_max[0]])
    ax.set_ylim3d([r_min[1], r_max[1]])
    ax.set_zlim3d([r_min[2], r_max[2]])
    ax.grid(False)
    plt.axis('equal')
    plt.axis('off')

    # make the panes transparent (the plot box)
    ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
    plt.show()
