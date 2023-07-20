import re
import typing
from collections import deque

import numpy as np
import rdkit.Chem.AllChem as rdkit
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation

# Holds the elements Van der Waals radii in Angstroms.
atom_vdw_radii = {
    "Al": 2,
    "Sb": 2,
    "Ar": 1.88,
    "As": 1.85,
    "Ba": 2,
    "Be": 2,
    "Bi": 2,
    "B": 2,
    "Br": 1.85,
    "Cd": 1.58,
    "Cs": 2,
    "Ca": 2,
    "C": 1.7,
    "Ce": 2,
    "Cl": 1.75,
    "Cr": 2,
    "Co": 2,
    "Cu": 1.4,
    "Dy": 2,
    "Er": 2,
    "Eu": 2,
    "F": 1.47,
    "Gd": 2,
    "Ga": 1.87,
    "Ge": 2,
    "Au": 1.66,
    "Hf": 2,
    "He": 1.4,
    "Ho": 2,
    "H": 1.09,
    "In": 1.93,
    "I": 1.98,
    "Ir": 2,
    "Fe": 2,
    "Kr": 2.02,
    "La": 2,
    "Pb": 2.02,
    "Li": 1.82,
    "Lu": 2,
    "Mg": 1.73,
    "Mn": 2,
    "Hg": 1.55,
    "Mo": 2,
    "Nd": 2,
    "Ne": 1.54,
    "Ni": 1.63,
    "Nb": 2,
    "N": 1.55,
    "Os": 2,
    "O": 1.52,
    "Pd": 1.63,
    "P": 1.8,
    "Pt": 1.72,
    "K": 2.75,
    "Pr": 2,
    "Pa": 2,
    "Re": 2,
    "Rh": 2,
    "Rb": 2,
    "Ru": 2,
    "Sm": 2,
    "Sc": 2,
    "Se": 1.9,
    "Si": 2.1,
    "Ag": 1.72,
    "Na": 2.27,
    "Sr": 2,
    "S": 1.8,
    "Ta": 2,
    "Te": 2.06,
    "Tb": 2,
    "Tl": 1.96,
    "Th": 2,
    "Tm": 2,
    "Sn": 2.17,
    "Ti": 2,
    "W": 2,
    "U": 1.86,
    "V": 2,
    "Xe": 2.16,
    "Yb": 2,
    "Y": 2,
    "Zn": 1.29,
    "Zr": 2,
    "X": 1.0,
    "D": 1.0,
}

# This dictionary gives easy access to the rdkit bond types.
bond_dict = {
    "1": rdkit.rdchem.BondType.SINGLE,
    "am": rdkit.rdchem.BondType.SINGLE,
    "2": rdkit.rdchem.BondType.DOUBLE,
    "3": rdkit.rdchem.BondType.TRIPLE,
    "ar": rdkit.rdchem.BondType.AROMATIC,
}

# A dictionary which matches atomic number to elemental symbols.
periodic_table = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
    8: "O",
    9: "F",
    10: "Ne",
    11: "Na",
    12: "Mg",
    13: "Al",
    14: "Si",
    15: "P",
    16: "S",
    17: "Cl",
    18: "Ar",
    19: "K",
    20: "Ca",
    21: "Sc",
    22: "Ti",
    23: "V",
    24: "Cr",
    25: "Mn",
    26: "Fe",
    27: "Co",
    28: "Ni",
    29: "Cu",
    30: "Zn",
    31: "Ga",
    32: "Ge",
    33: "As",
    34: "Se",
    35: "Br",
    36: "Kr",
    37: "Rb",
    38: "Sr",
    39: "Y",
    40: "Zr",
    41: "Nb",
    42: "Mo",
    43: "Tc",
    44: "Ru",
    45: "Rh",
    46: "Pd",
    47: "Ag",
    48: "Cd",
    49: "In",
    50: "Sn",
    51: "Sb",
    52: "Te",
    53: "I",
    54: "Xe",
    55: "Cs",
    56: "Ba",
    57: "La",
    58: "Ce",
    59: "Pr",
    60: "Nd",
    61: "Pm",
    62: "Sm",
    63: "Eu",
    64: "Gd",
    65: "Tb",
    66: "Dy",
    67: "Ho",
    68: "Er",
    69: "Tm",
    70: "Yb",
    71: "Lu",
    72: "Hf",
    73: "Ta",
    74: "W",
    75: "Re",
    76: "Os",
    77: "Ir",
    78: "Pt",
    79: "Au",
    80: "Hg",
    81: "Tl",
    82: "Pb",
    83: "Bi",
    84: "Po",
    85: "At",
    86: "Rn",
    87: "Fr",
    88: "Ra",
    89: "Ac",
    90: "Th",
    91: "Pa",
    92: "U",
    93: "Np",
    94: "Pu",
    95: "Am",
    96: "Cm",
    97: "Bk",
    98: "Cf",
    99: "Es",
    100: "Fm",
    101: "Md",
    102: "No",
    103: "Lr",
    104: "Rf",
    105: "Db",
    106: "Sg",
    107: "Bh",
    108: "Hs",
    109: "Mt",
    110: "Ds",
    111: "Rg",
    112: "Cn",
    113: "Uut",
    114: "Fl",
    115: "Uup",
    116: "Lv",
    117: "Uus",
    118: "Uuo",
}


T = typing.TypeVar("T")
OneOrMany = typing.Union[T, typing.Iterable[T]]


class Cell:
    """
    Represents an individual cell in a supercell.

    Parameters
    ----------
    id : :class:`list` of :class:`int`
        A 3 member :class:`list` holding the x, y and z index
        of the cell within the supercell.

    fgs : :class:`dict`
        Maps the fgs in the original unit cell to the
        equivalent fgs in the cell.

    Attributes
    ----------
    id : :class:`numpy.ndarray` of :class:`int`
        A 3 member array holding the x, y and z index
        of the cell within the supercell.

    fgs : :class:`dict`
        Maps the fgs in the original unit cell to the
        equivalent fgs in the cell.

    """

    def __init__(self, id_, fgs):
        self.id = np.array(id_)
        self.fgs = fgs


def dedupe(iterable, key=None, seen=None):
    """
    Yields items from `iterable` barring duplicates.

    If `seen` is provided it contains elements which are not to be
    yielded at all.

    Parameters
    ----------
    iterable : :class:`iterable`
        An iterable of elements which are to be yielded, only once.

    key : :class:`callable`
        A function which gets applied to every member of `iterable`.
        The return of this function is checked for duplication rather
        than the member itself.

    seen : :class:`set`, optional
        Holds items which are not to be yielded.

    Yields
    ------
    :class:`object`
        Element in `iterable` which is not in `seen` and has not been
        yielded before.

    """

    if seen is None:
        seen = set()
    for x in iterable:
        val = key(x) if key is not None else x
        if val not in seen:
            seen.add(val)
            yield x


def flatten(iterable, excluded_types=None):
    """
    Transforms an nested iterable into a flat one.

    Examples
    --------
    *Flattening a Deeply Nested Structure*

    .. testcode:: flattening-a-deeply-nested-structure

        import stk

        nested = [[1, 2, 3], [[4], [5], [[6]], 7]]
        print(list(stk.flatten(nested)))

    gives

    .. testoutput:: flattening-a-deeply-nested-structure

        [1, 2, 3, 4, 5, 6, 7]

    *Avoiding Flattening of Some Types*

    If a type is found in `excluded_types` it will not be yielded from.
    For example, ``str`` is in `excluded_types` by default, so

    .. testcode:: avoiding-flattening-of-some-types

        import stk

        nested = ['abcd', ['efgh']]
        print(list(stk.flatten(nested)))

    gives

    .. testoutput:: avoiding-flattening-of-some-types

        ['abcd', 'efgh']

    instead of

    .. code-block::

        ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h']

    Parameters
    ----------
    iterable : :class:`iterable`
        The iterable which is to be flattened

    excluded_types : :class:`set`, optional
        Holds container types which are not be flattened.

    Yields
    ------
    :class:`object`
        A nested element of `iterable`.

    """

    if excluded_types is None:
        excluded_types = {str}

    for x in iterable:
        if hasattr(x, "__iter__") and type(x) not in excluded_types:
            yield from flatten(x, excluded_types)
        else:
            yield x


def mol_from_mae_file(mae_path):
    """
    Creates a ``rdkit`` molecule from a ``.mae`` file.

    Parameters
    ----------
    mol2_file : :class:`str`
        The full path of the ``.mae`` file from which an rdkit molecule
        should be instantiated.

    Returns
    -------
    :class:`rdkit.Mol`
        An ``rdkit`` instance of the molecule held in `mae_file`.

    """

    mol = rdkit.EditableMol(rdkit.Mol())
    conf = rdkit.Conformer()

    with open(mae_path, "r") as mae:
        content = re.split(r"[{}]", mae.read())

    prev_block = deque([""], maxlen=1)
    for block in content:
        if "m_atom[" in prev_block[0]:
            atom_block = block
        if "m_bond[" in prev_block[0]:
            bond_block = block
        prev_block.append(block)

    labels, data_block, *_ = atom_block.split(":::")
    labels = [
        label
        for label in labels.split("\n")
        if not label.isspace() and label != ""
    ]

    data_block = [
        a.split()
        for a in data_block.split("\n")
        if not a.isspace() and a != ""
    ]

    for line in data_block:
        line = [word for word in line if word != '"']
        if len(labels) != len(line):
            raise RuntimeError(
                (
                    "Number of labels does"
                    " not match number of columns"
                    " in .mae file."
                )
            )

        for label, data in zip(labels, line):
            if "x_coord" in label:
                x = float(data)
            if "y_coord" in label:
                y = float(data)
            if "z_coord" in label:
                z = float(data)
            if "atomic_number" in label:
                atom_num = int(data)

        atom_sym = periodic_table[atom_num]
        atom_coord = Point3D(x, y, z)
        atom_id = mol.AddAtom(rdkit.Atom(atom_sym))
        conf.SetAtomPosition(atom_id, atom_coord)

    labels, data_block, *_ = bond_block.split(":::")
    labels = [
        label
        for label in labels.split("\n")
        if not label.isspace() and label != ""
    ]
    data_block = [
        a.split()
        for a in data_block.split("\n")
        if not a.isspace() and a != ""
    ]

    for line in data_block:
        if len(labels) != len(line):
            raise RuntimeError(
                (
                    "Number of labels does"
                    " not match number of "
                    "columns in .mae file."
                )
            )

        for label, data in zip(labels, line):
            if "from" in label:
                atom1 = int(data) - 1
            if "to" in label:
                atom2 = int(data) - 1
            if "order" in label:
                bond_order = str(int(data))
        mol.AddBond(atom1, atom2, bond_dict[bond_order])

    mol = mol.GetMol()
    mol.AddConformer(conf)
    return mol


def normalize_vector(vector):
    """
    Normalizes the given vector.

    A new vector is returned, the original vector is not modified.

    Parameters
    ----------
    vector : :class:`np.ndarray`
        The vector to be normalized.

    Returns
    -------
    :class:`np.ndarray`
        The normalized vector.

    """

    return np.divide(vector, np.linalg.norm(vector))


def remake(molcule: rdkit.Mol) -> rdkit.Mol:
    emol = rdkit.EditableMol(rdkit.Mol())
    for atom in molcule.GetAtoms():
        new = rdkit.Atom(atom.GetAtomicNum())
        new.SetFormalCharge(atom.GetFormalCharge())
        emol.AddAtom(new)

    for bond in molcule.GetBonds():
        emol.AddBond(
            beginAtomIdx=bond.GetBeginAtomIdx(),
            endAtomIdx=bond.GetEndAtomIdx(),
            order=bond.GetBondType(),
        )

    m = emol.GetMol()
    m.AddConformer(molcule.GetConformer())
    return m


def get_projection(start, target):
    """
    Get the projection of `start` onto `target`.

    """

    return (
        target
        * np.dot(
            start,
            target,
        )
        / np.dot(target, target)
    )


def orthogonal_vector(vector):
    ortho = np.array([0.0, 0.0, 0.0])
    for m, val in enumerate(vector):
        if not np.allclose(val, 0, atol=1e-8):
            n = (m + 1) % 3
            break
    ortho[n] = vector[m]
    ortho[m] = -vector[n]
    return ortho


def rotation_matrix(vector1, vector2):
    """
    Returns a rotation matrix which transforms `vector1` to `vector2`.

    Multiplying `vector1` by the rotation matrix returned by this
    function yields `vector2`.

    Parameters
    ----------
    vector1 : :class:`numpy.ndarray`
        The vector which needs to be transformed to `vector2`.

    vector2 : :class:`numpy.ndarray`
        The vector onto which `vector1` needs to be transformed.

    Returns
    -------
    :class:`numpy.ndarray`
        A rotation matrix which transforms `vector1` to `vector2`.

    References
    ----------
    http://tinyurl.com/kybj9ox
    http://tinyurl.com/gn6e8mz

    """

    # Make sure both inputs are unit vectors.
    vector1 = normalize_vector(vector1)
    vector2 = normalize_vector(vector2)

    # Hande the case where the input and output vectors are equal.
    if np.allclose(vector1, vector2, atol=1e-8):
        return np.identity(3)

    # Handle the case where the rotation is 180 degrees.
    if np.allclose(vector1, np.multiply(vector2, -1), atol=1e-8):
        return rotation_matrix_arbitrary_axis(
            angle=np.pi, axis=orthogonal_vector(vector1)
        )

    v = np.cross(vector1, vector2)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    s = np.linalg.norm(v)
    c = np.dot(vector1, vector2)
    i = np.identity(3)
    mult_factor = (1 - c) / np.square(s)

    # Initialize as a scipy Rotation object, which normalizes the
    # matrix and allows for returns as quaternion or alternative
    # type in the future.
    return Rotation.from_matrix(
        i + vx + np.multiply(np.dot(vx, vx), mult_factor)
    ).as_matrix()


def rotation_matrix_arbitrary_axis(angle, axis):
    """
    Returns a rotation matrix of `angle` radians about `axis`.

    Parameters
    ----------
    angle : :class:`float`
        The size of the rotation in radians.

    axis : :class:`numpy.ndarray`
        A 3 element aray which represents a vector. The vector is the
        axis about which the rotation is carried out. Must be of
        unit magnitude.

    Returns
    -------
    :class:`numpy.ndarray`
        A ``3x3`` array representing a rotation matrix.

    """

    a = np.cos(angle / 2)
    b, c, d = axis * np.sin(angle / 2)

    e11 = np.square(a) + np.square(b) - np.square(c) - np.square(d)
    e12 = 2 * (b * c - a * d)
    e13 = 2 * (b * d + a * c)

    e21 = 2 * (b * c + a * d)
    e22 = np.square(a) + np.square(c) - np.square(b) - np.square(d)
    e23 = 2 * (c * d - a * b)

    e31 = 2 * (b * d - a * c)
    e32 = 2 * (c * d + a * b)
    e33 = np.square(a) + np.square(d) - np.square(b) - np.square(c)

    # Initialize as a scipy Rotation object, which normalizes the
    # matrix and allows for returns as quaternion or alternative
    # type in the future.
    return Rotation.from_matrix(
        np.array([[e11, e12, e13], [e21, e22, e23], [e31, e32, e33]])
    ).as_matrix()


def dice_similarity(mol1, mol2, fp_radius=3):
    """
    Return the chemical similarity between two molecules.

    Parameters
    ----------
    mol1 : :class:`.Molecule`
        The first molecule.

    mol2 : :class:`.Molecule`
        The second molecule.

    fp_radius : :class:`int`, optional
        The radius of the Morgan fingerprint used to calculate
        similarity.

    Returns
    -------
    :class:`float`
        The similarity.

    """

    rdkit_mol1 = mol1.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol1)
    fp1 = rdkit.GetMorganFingerprint(
        mol=rdkit_mol1,
        radius=fp_radius,
    )
    rdkit_mol2 = mol2.to_rdkit_mol()
    rdkit.SanitizeMol(rdkit_mol2)
    fp2 = rdkit.GetMorganFingerprint(
        mol=rdkit_mol2,
        radius=fp_radius,
    )
    return rdkit.DataStructs.DiceSimilarity(fp1, fp2)


def vector_angle(vector1, vector2):
    """
    Returns the angle between two vectors in radians.

    Parameters
    ----------
    vector1 : :class:`numpy.ndarray`
        The first vector.

    vector2 : :class:`numpy.ndarray`
        The second vector.

    Returns
    -------
    :class:`float`
        The angle between `vector1` and `vector2` in radians.

    """

    if np.all(np.equal(vector1, vector2)):
        return 0.0

    numerator = np.dot(vector1, vector2)
    denominator = np.linalg.norm(vector1) * np.linalg.norm(vector2)
    # This if statement prevents returns of NaN due to floating point
    # inaccuracy.
    term = numerator / denominator
    if term >= 1.0:
        return 0.0
    if term <= -1.0:
        return np.pi
    return np.arccos(term)


def get_acute_vector(reference, vector):
    if (
        # vector_angle is NaN if reference is [0, 0, 0].
        not np.allclose(reference, [0, 0, 0], atol=1e-5)
        and vector_angle(vector, reference) > np.pi / 2
    ):
        return vector * -1
    return vector


def get_plane_normal(points):
    centroid = points.sum(axis=0) / len(points)
    return np.around(np.linalg.svd(points - centroid)[-1][2, :], 14)


def cap_absolute_value(value, max_absolute_value=1):
    """
    Returns `value` with absolute value capped at `max_absolute_value`.

    Particularly useful in passing values to trignometric functions
    where numerical errors may result in an argument > 1 being passed
    in.

    This code is modified from the pymatgen source code [1]_.

    Parameters
    ----------
    value : :class:`float`
        Value to cap.

    max_absolute_value : :class:`float`, optional
        Absolute value to cap `value` at.
        Defaults to 1.

    Returns
    -------
    :class:`float`
        `value` capped at `max_absolute_value` with sign preserved.

    References
    ----------
    .. [1] https://pymatgen.org/pymatgen.util.num.html

    """

    return max(min(value, max_absolute_value), -max_absolute_value)
