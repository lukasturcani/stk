"""
This module defines general-purpose objects, functions and classes.

Functions, classes, etc. defined here should not depend on any other
part of ``stk``. They must be completely self-sufficient.

"""

import rdkit.Chem.AllChem as rdkit
from rdkit.Geometry import Point3D
import numpy as np
import time
from contextlib import contextmanager
import os
import subprocess as sp
import gzip
import re
from collections import deque
import tarfile
from glob import iglob


# Holds the elements Van der Waals radii in Angstroms.
atom_vdw_radii = {
    'Al': 2, 'Sb': 2, 'Ar': 1.88, 'As': 1.85, 'Ba': 2,
    'Be': 2, 'Bi': 2, 'B': 2, 'Br': 1.85, 'Cd': 1.58,
    'Cs': 2, 'Ca': 2, 'C': 1.7, 'Ce': 2, 'Cl': 1.75,
    'Cr': 2, 'Co': 2, 'Cu': 1.4, 'Dy': 2, 'Er': 2,
    'Eu': 2, 'F': 1.47, 'Gd': 2, 'Ga': 1.87, 'Ge': 2,
    'Au': 1.66, 'Hf': 2, 'He': 1.4, 'Ho': 2, 'H': 1.09,
    'In': 1.93, 'I': 1.98, 'Ir': 2, 'Fe': 2, 'Kr': 2.02,
    'La': 2, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2, 'Mg': 1.73,
    'Mn': 2, 'Hg': 1.55, 'Mo': 2, 'Nd': 2, 'Ne': 1.54,
    'Ni': 1.63, 'Nb': 2, 'N': 1.55, 'Os': 2, 'O': 1.52,
    'Pd': 1.63, 'P': 1.8, 'Pt': 1.72, 'K': 2.75, 'Pr': 2,
    'Pa': 2, 'Re': 2, 'Rh': 2, 'Rb': 2, 'Ru': 2, 'Sm': 2,
    'Sc': 2, 'Se': 1.9, 'Si': 2.1, 'Ag': 1.72, 'Na': 2.27,
    'Sr': 2, 'S': 1.8, 'Ta': 2, 'Te': 2.06, 'Tb': 2,
    'Tl': 1.96, 'Th': 2, 'Tm': 2, 'Sn': 2.17, 'Ti': 2,
    'W': 2, 'U': 1.86, 'V': 2, 'Xe': 2.16, 'Yb': 2,
    'Y': 2, 'Zn': 1.29, 'Zr': 2, 'X': 1.0, 'D': 1.0
}

# This dictionary gives easy access to the rdkit bond types.
bond_dict = {'1': rdkit.rdchem.BondType.SINGLE,
             'am': rdkit.rdchem.BondType.SINGLE,
             '2': rdkit.rdchem.BondType.DOUBLE,
             '3': rdkit.rdchem.BondType.TRIPLE,
             'ar': rdkit.rdchem.BondType.AROMATIC}

# A dictionary which matches atomic number to elemental symbols.
periodic_table = {
    1: 'H', 2: 'He', 3: 'Li', 4: 'Be', 5: 'B', 6: 'C',
    7: 'N', 8: 'O', 9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
    13: 'Al', 14: 'Si', 15: 'P', 16: 'S', 17: 'Cl',
    18: 'Ar', 19: 'K', 20: 'Ca', 21: 'Sc', 22: 'Ti',
    23: 'V', 24: 'Cr', 25: 'Mn', 26: 'Fe', 27: 'Co',
    28: 'Ni', 29: 'Cu', 30: 'Zn', 31: 'Ga', 32: 'Ge',
    33: 'As', 34: 'Se', 35: 'Br', 36: 'Kr', 37: 'Rb',
    38: 'Sr', 39: 'Y', 40: 'Zr', 41: 'Nb', 42: 'Mo',
    43: 'Tc', 44: 'Ru', 45: 'Rh', 46: 'Pd', 47: 'Ag',
    48: 'Cd', 49: 'In', 50: 'Sn', 51: 'Sb', 52: 'Te',
    53: 'I', 54: 'Xe', 55: 'Cs', 56: 'Ba', 57: 'La',
    58: 'Ce', 59: 'Pr', 60: 'Nd', 61: 'Pm', 62: 'Sm',
    63: 'Eu', 64: 'Gd', 65: 'Tb', 66: 'Dy', 67: 'Ho',
    68: 'Er', 69: 'Tm', 70: 'Yb', 71: 'Lu', 72: 'Hf',
    73: 'Ta', 74: 'W', 75: 'Re', 76: 'Os', 77: 'Ir',
    78: 'Pt', 79: 'Au', 80: 'Hg', 81: 'Tl', 82: 'Pb',
    83: 'Bi', 84: 'Po', 85: 'At', 86: 'Rn', 87: 'Fr',
    88: 'Ra', 89: 'Ac', 90: 'Th', 91: 'Pa', 92: 'U',
    93: 'Np', 94: 'Pu', 95: 'Am', 96: 'Cm', 97: 'Bk',
    98: 'Cf', 99: 'Es', 100: 'Fm', 101: 'Md', 102: 'No',
    103: 'Lr', 104: 'Rf', 105: 'Db', 106: 'Sg', 107: 'Bh',
    108: 'Hs', 109: 'Mt', 110: 'Ds', 111: 'Rg', 112: 'Cn',
    113: 'Uut', 114: 'Fl', 115: 'Uup', 116: 'Lv',
    117: 'Uus', 118: 'Uuo'
}


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


class ChargedMolError(Exception):
    def __init__(self, mol_file, msg):
        self.mol_file = mol_file
        self.msg = msg


class MolFileError(Exception):
    def __init__(self, mol_file, msg):
        self.mol_file = mol_file
        self.msg = msg


class PopulationSizeError(Exception):
    def __init__(self, msg):
        self.msg = msg


class LazyAttr:
    """
    A descriptor for creating lazy attributes.

    """

    def __init__(self, func):
        self.func = func

    def __get__(self, obj, cls):
        if obj is None:
            return self
        val = self.func(obj)
        setattr(obj, self.func.__name__, val)
        return val


class MAEExtractor:
    """
    Extracts the lowest energy conformer from a .maegz file.

    Macromodel conformer searches produce -out.maegz files containing
    all of the conformers found during the search and their energies
    and other data.

    Initializing this class with a :class:`.ConstructedMolecule` finds
    the ``-out.maegz`` file of that :class:`.ConstructedMolecule` and
    converts it to a ``.mae`` file. It then creates and additional
    ``.mae`` file holding only the lowest energy conformer found.

    Attributes
    ----------
    maegz_path : :class:`str`
        The path to the ``-out.maegz`` file generated by the macromodel
        conformer search.

    mae_path : :class:`str`
        The path to the ``.mae`` file holding the conformers generated
        by the macromodel conformer search.

    content : :class:`str`
        The content of the ``.mae`` file hodling all the conformers
        from the macromodel conformer search. This holds other data
        such as their energies too.

    energies : :class:`list`
        The :class:`list` has the form

        .. code-block:: python

            energies = [(0, 231.0), (1, 144.4), ...]

        Each :class:`tuple` holds the id and energy of every conformer
        in the ``.mae`` file, respectively.

    min_energy : :class:`float`
        The minimum energy found in the ``.mae`` file.

    path : :class:`str`
        The full path of the ``.mae`` file holding the extracted lowest
        energy conformer.

    """

    def __init__(self, run_name, n=1):
        self.maegz_path = f'{run_name}-out.maegz'
        self.maegz_to_mae()
        self.extract_conformers(n)

    def extract_conformers(self, n):
        """
        Creates ``.mae`` files holding the lowest energy conformers.

        Parameters
        ----------
        n : :class:`int`
            The number of conformers to extract.

        Returns
        -------
        None : :class:`NoneType`

        """

        for i in range(n):
            # Get the id of the lowest energy conformer.
            num = self.lowest_energy_conformers(n)[i][1]
            # Get the structure block corresponding to the lowest
            # energy conformer.
            content = self.content.split("f_m_ct")
            new_mae = "f_m_ct".join([content[0], content[num]])

            # Write the structure block in its own .mae file, named
            # after conformer extracted.
            if n == 1:
                # Write the structure block in its own .mae file, named
                # after conformer extracted.
                new_name = self.mae_path.replace(
                    '.mae',
                    f'EXTRACTED_{num}.mae'
                )
            else:
                new_name = self.mae_path.replace(
                    '.mae',
                    f'EXTRACTED_{num}_conf_{i}.mae'
                )

            with open(new_name, 'w') as mae_file:
                mae_file.write(new_mae)

            if i == 0:
                # Save the path of the newly created file.
                self.path = new_name

    def extract_energy(self, block):
        """
        Extracts the energy value from a ``.mae`` energy data block.

        Parameters
        ----------
        block : :class:`str`
            An ``.mae`` energy data block.

        Returns
        -------
        :class:`float`
            The energy value extracted from `block` or ``None`` if
            one is not found.

        """

        block = block.split(":::")
        for name, value in zip(block[0].split('\n'),
                               block[1].split('\n')):
            if 'r_mmod_Potential_Energy' in name:
                return float(value)

    def lowest_energy_conformers(self, n):
        """
        Returns the id and energy of the lowest energy conformers.

        Parameters
        ----------
        n : :class:`int`
            The number of lowest energy conformers to return.

        Returns
        -------
        :class:`list`
            A :class:`list` of the form

            .. code-block:: python

                returned = [(23, 123.3), (1, 143.89), (12, 150.6), ...]

            Where each :class:`tuple` holds the id and energy of the
            `n` lowest energy conformers, respectively.

        """

        # Open the .mae file holding all the conformers and load its
        # content.
        with open(self.mae_path, 'r') as mae_file:
            self.content = mae_file.read()
            # Split the content across curly braces. This divides the
            # various sections of the .mae file.
            content_split = re.split(r"[{}]", self.content)

        # Go through all the datablocks in the the .mae file. For each
        # energy block extract the energy and store it in the
        # `energies` list. Store the `index`  (conformer id) along with
        # each extracted energy.
        self.energies = []
        prev_block = deque([""], maxlen=1)
        index = 1
        for block in content_split:
            if ("f_m_ct" in prev_block[0] and
               "r_mmod_Potential_Energy" in block):
                energy = self.extract_energy(block)
                self.energies.append((energy, index))
                index += 1

            prev_block.append(block)

        # Selecting the lowest energy n conformers
        confs = sorted(self.energies)[:n]
        # Define the energy of the lowest energy conformer
        self.min_energy = confs[0][0]
        # Return a list with id and energy of the lowest energy
        # conformers.
        return confs

    def maegz_to_mae(self):
        """
        Converts the .maegz file to a .mae file.

        Returns
        -------
        None : :class:`NoneType`

        """

        self.mae_path = self.maegz_path.replace('.maegz', '.mae')
        with gzip.open(self.maegz_path, 'r') as maegz_file:
            with open(self.mae_path, 'wb') as mae_file:
                mae_file.write(maegz_file.read())


def archive_output():
    """
    Places the ``output`` folder into ``old_output``.

    This function assumes that the ``output`` folder is in the current
    working directory.

    Returns
    -------
    None : :class:`NoneType`

    """

    if 'output' not in os.listdir():
        return

    # Make the ``old_output`` folder if it does not exist already.
    if 'old_output' not in os.listdir():
        os.mkdir('old_output')

    # Find out with what number the ``output`` folder should be
    # labelled within ``old_output``.
    num = len(os.listdir('old_output'))
    new_dir = os.path.join('old_output', str(num))
    os.rename('output', new_dir)


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

    For example

    .. code-block:: python

        [[1,2,3], [[4], [5],[[6]], 7]

    becomes

    .. code-block:: python

        [1,2,3,4,5,6,7]

    If a type is found in `excluded_types` it will not be yielded from.
    For example if ``str`` is in `excluded_types`

    .. code-block:: python

        a = ["abcd", ["efgh"]]

    "abcd" and "efgh" are yielded if `a` is passed to `iterable`. If
    `str` was not in `excluded_types` then "a", "b", "c", "d", "e",
    "f", "g" and "h" would all be yielded individually.

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
        if hasattr(x, '__iter__') and type(x) not in excluded_types:
            yield from flatten(x, excluded_types)
        else:
            yield x


def kabsch(coords1, coords2):
    """
    Return a rotation matrix to minimize dstance between 2 coord sets.

    This is essentially an implementation of the Kabsch algorithm.
    Given two sets of coordinates, `coords1` and `coords2`, this
    function returns a rotation matrix. When the rotation matrix is
    applied to `coords1` the resulting coordinates have their rms
    distance to `coords2` minimized.

    Parameters
    ----------
    coords1 : :class:`numpy.ndarray`
        This array represents a matrix hodling coordinates which need
        to be rotated to minimize their rms distance to coordinates in
        `coords2`. The matrix is n x 3. Each row of the matrix holds
        the x, y and z coordinates of one point, respectively. Here
        ``n`` is the number of points.

    coords2 : :class:`numpy.ndarray`
        This array represents a matrix which holds the coordinates of
        points the distance to which should be minimized. The matrix is
        ``n x 3``. Each row of the matrix holds the x, y and z
        coordinates of one point, respectively. Here ``n`` is the
        number of points.

    Returns
    -------
    :class:`numpy.ndarray`
        A rotation matrix. This will be a ``3 x 3`` matrix.

    References
    ----------
    http://nghiaho.com/?page_id=671
    https://en.wikipedia.org/wiki/Kabsch_algorithm

    """

    h = np.dot(coords1, coords2.T)
    u, s, v = np.linalg.svd(h)

    if int(np.linalg.det(v)) < 0:
        v[:, 2] = -v[:, 2]

    return np.dot(v, u)


def kill_macromodel():
    """
    Kills any applications left open as a result running MacroModel.

    Applications that are typically left open are
    ``jserver-watcher.exe`` and ``jservergo.exe``.

    Returns
    -------
    None : :class:`NoneType`

    """

    if os.name == 'nt':
        # In Windows, use the ``Taskkill`` command to force a close on
        # the applications.
        sp.run(["Taskkill", "/IM", "jserver-watcher.exe", "/F"],
               stdout=sp.PIPE, stderr=sp.PIPE)
        sp.run(["Taskkill", "/IM", "jservergo.exe", "/F"],
               stdout=sp.PIPE, stderr=sp.PIPE)

    if os.name == 'posix':
        sp.run(["pkill", "jservergo"],
               stdout=sp.PIPE, stderr=sp.PIPE)
        sp.run(["pkill", "jserver-watcher"],
               stdout=sp.PIPE, stderr=sp.PIPE)


def matrix_centroid(matrix):
    """
    Returns the centroid of the coordinates held in `matrix`.

    Parameters
    ----------
    matrix : :class:`np.ndarray`
        A ``n x 3`` matrix. Each row holds the x, y and z
        coordinate of some point, respectively.

    Returns
    -------
    :class:`numpy.ndarray`
        A numpy array which holds the x, y and z
        coordinates of the centroid of the coordinates in `matrix`.

    """

    return np.sum(matrix, axis=0) / len(matrix)


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

    with open(mae_path, 'r') as mae:
        content = re.split(r'[{}]', mae.read())

    prev_block = deque([''], maxlen=1)
    for block in content:
        if 'm_atom[' in prev_block[0]:
            atom_block = block
        if 'm_bond[' in prev_block[0]:
            bond_block = block
        prev_block.append(block)

    labels, data_block, *_ = atom_block.split(':::')
    labels = [label for label in labels.split('\n')
              if not label.isspace() and label != '']

    data_block = [a.split() for a in data_block.split('\n') if
                  not a.isspace() and a != '']

    for line in data_block:
        line = [word for word in line if word != '"']
        if len(labels) != len(line):
            raise RuntimeError(('Number of labels does'
                                ' not match number of columns'
                                ' in .mae file.'))

        for label, data in zip(labels, line):
            if 'x_coord' in label:
                x = float(data)
            if 'y_coord' in label:
                y = float(data)
            if 'z_coord' in label:
                z = float(data)
            if 'atomic_number' in label:
                atom_num = int(data)

        atom_sym = periodic_table[atom_num]
        atom_coord = Point3D(x, y, z)
        atom_id = mol.AddAtom(rdkit.Atom(atom_sym))
        conf.SetAtomPosition(atom_id, atom_coord)

    labels, data_block, *_ = bond_block.split(':::')
    labels = [label for label in labels.split('\n')
              if not label.isspace() and label != '']
    data_block = [a.split() for a in data_block.split('\n')
                  if not a.isspace() and a != '']

    for line in data_block:
        if len(labels) != len(line):
            raise RuntimeError(('Number of labels does'
                                ' not match number of '
                                'columns in .mae file.'))

        for label, data in zip(labels, line):
            if 'from' in label:
                atom1 = int(data) - 1
            if 'to' in label:
                atom2 = int(data) - 1
            if 'order' in label:
                bond_order = str(int(data))
        mol.AddBond(atom1, atom2, bond_dict[bond_order])

    mol = mol.GetMol()
    mol.AddConformer(conf)
    return mol


def mol_from_mol_file(mol_file):
    """
    Creates a rdkit molecule from a ``.mol`` (V3000) file.

    Parameters
    ----------
    mol_file : :class:`str`
        The full of the .mol file from which an rdkit molecule should
        be instantiated.

    Returns
    -------
    :class:`rdkit.Mol`
        An rdkit instance of the molecule held in `mol2_file`.

    Raises
    ------
    :class:`ChargedMolError`
        If an atom row has more than 8 coloumns it is usually because
        there is a 9th coloumn indicating atomic charge. Such molecules
        are not currently supported, so an error is raised.

    :class:`MolFileError`
        If the file is not a V3000 ``.mol`` file.

    """

    e_mol = rdkit.EditableMol(rdkit.Mol())
    conf = rdkit.Conformer()

    with open(mol_file, 'r') as f:
        take_atom = False
        take_bond = False
        v3000 = False

        for line in f:
            if 'V3000' in line:
                v3000 = True

            if 'M  V30 BEGIN ATOM' in line:
                take_atom = True
                continue

            if 'M  V30 END ATOM' in line:
                take_atom = False
                continue

            if 'M  V30 BEGIN BOND' in line:
                take_bond = True
                continue

            if 'M  V30 END BOND' in line:
                take_bond = False
                continue

            if take_atom:
                words = line.split()
                if len(words) > 8:
                    raise ChargedMolError(mol_file,
                                          ('Atom row has more'
                                           ' than 8 coloumns. Likely '
                                           'due to a charged atom.'))
                _, _, _, atom_sym, *coords, _ = words
                coords = [float(x) for x in coords]
                atom_coord = Point3D(*coords)
                atom_id = e_mol.AddAtom(rdkit.Atom(atom_sym))
                conf.SetAtomPosition(atom_id, atom_coord)
                continue

            if take_bond:
                *_, bond_id, bond_order, atom1, atom2 = line.split()
                e_mol.AddBond(int(atom1)-1, int(atom2)-1,
                              bond_dict[bond_order])
                continue
    if not v3000:
        raise MolFileError(mol_file, 'Not a V3000 .mol file.')

    mol = e_mol.GetMol()
    mol.AddConformer(conf)
    return mol


def move_generated_macromodel_files(basename, output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    for filename in iglob(f'{basename}*'):
        # Do not move the output_dir.
        if filename == output_dir:
            continue
        os.rename(filename, f'{output_dir}/{filename}')


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


def remake(mol):
    """
    Remakes a molecule from scratch.

    Parameters
    ----------
    mol : :class:`rdkit.Mol`
        The molecule to be remade.

    Returns
    -------
    :class:`rdkit.Mol`
        The remade molecule.

    """

    emol = rdkit.EditableMol(rdkit.Mol())
    for atom in mol.GetAtoms():
        new = rdkit.Atom(atom.GetAtomicNum())
        new.SetFormalCharge(atom.GetFormalCharge())
        emol.AddAtom(new)

    for bond in mol.GetBonds():
        emol.AddBond(
            beginAtomIdx=bond.GetBeginAtomIdx(),
            endAtomIdx=bond.GetEndAtomIdx(),
            order=bond.GetBondType()
        )

    m = emol.GetMol()
    m.AddConformer(mol.GetConformer())
    return m


def orthogonal_vector(vector):
    ortho = [0, 0, 0]
    for m, val in enumerate(vector):
        if not np.allclose(val, 0, atol=1e-8):
            n = (m+1) % 3
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
            angle=np.pi,
            axis=orthogonal_vector(vector1)
        )

    v = np.cross(vector1, vector2)
    vx = np.array([
        [0, -v[2], v[1]],
        [v[2], 0, -v[0]],
        [-v[1], v[0], 0]
    ])
    s = np.linalg.norm(v)
    c = np.dot(vector1, vector2)
    i = np.identity(3)
    mult_factor = (1-c)/np.square(s)
    return i + vx + np.multiply(np.dot(vx, vx), mult_factor)


def rotation_matrix_arbitrary_axis(angle, axis):
    """
    Returns a rotation matrix of `angle` radians about `axis`.

    Parameters
    ----------
    angle : :class:`float`
        The size of the rotation in radians.

    axis : :class:`numpy.ndarray`
        A 3 element aray which represents a vector. The vector is the
        axis about which the rotation is carried out.

    Returns
    -------
    :class:`numpy.ndarray`
        A ``3x3`` array representing a rotation matrix.

    """

    axis = normalize_vector(axis)

    a = np.cos(angle/2)
    b, c, d = axis * np.sin(angle/2)

    e11 = np.square(a) + np.square(b) - np.square(c) - np.square(d)
    e12 = 2*(b*c - a*d)
    e13 = 2*(b*d + a*c)

    e21 = 2*(b*c + a*d)
    e22 = np.square(a) + np.square(c) - np.square(b) - np.square(d)
    e23 = 2*(c*d - a*b)

    e31 = 2*(b*d - a*c)
    e32 = 2*(c*d + a*b)
    e33 = np.square(a) + np.square(d) - np.square(b) - np.square(c)

    return np.array([[e11, e12, e13],
                     [e21, e22, e23],
                     [e31, e32, e33]])


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


def quaternion(u):
    """
    Returns a translation + rotation quaternion.

    Parameters
    ----------
    u : :class:`list` of :class:`float`
        A :class:`list` of length 3 holding the parameter for the
        quarternion.

    References
    ----------
    K. Shoemake, Uniform random rotations, Graphics Gems III,
    pages 124-132. Academic, New York, 1992.

    """

    a, b, c = u
    q = np.zeros(4, np.float64)
    q[0] = np.sqrt(1. - a) * np.sin(2. * np.pi * b)
    q[1] = np.sqrt(1. - a) * np.cos(2. * np.pi * b)
    q[2] = np.sqrt(a) * np.sin(2. * np.pi * c)
    q[3] = np.sqrt(a) * np.cos(2. * np.pi * c)
    return q


def translation_component(q):
    """
    Extracts translation vector from quaternion.

    Parameters
    ----------
    q : :class:`numpy.ndarray`
        A length 4 quaternion.

    Returns
    -------
    :class:`numpy.ndarray`
        The translation vector encoded within `q`.

    """

    rot_epsilon = 1e-6
    q = np.copy(q)
    if q[0] < 0.:
        q = -q
    if q[0] > 1.0:
        q /= np.sqrt(np.dot(q, q))

    theta = 2. * np.arccos(q[0])
    s = np.sqrt(1. - q[0] * q[0])
    if s < rot_epsilon:
        p = 2. * q[1:4]
    else:
        p = q[1:4] / s * theta
    return p


def tar_output():
    """
    Places all the content in the ``output`` folder into a .tgz file.

    Returns
    -------
    None : :class:`NoneType`

    """

    tname = os.path.join('output', 'output.tgz')
    with tarfile.open(tname, 'w:gz') as tar:
        tar.add('output')


@contextmanager
def time_it():
    """
    Times the code executed within the indent.

    This is a context manager so it should be used as:

    .. code-block:: python

        with time_it():
            something1()
            something2()
            something3()

    After all 3 functions are finished and the nested block is exited
    the time taken to process the entire block is printed.

    """

    start = time.time()
    yield
    time_taken = time.time() - start
    m, s = divmod(time_taken, 60)
    h, m = divmod(m, 60)
    print(f'\nTime taken was {int(h)} : {int(m)} : {int(s)}.\n\n')


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

    numerator = np.dot(vector1, vector2)
    denominator = (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    # This if statement prevents returns of NaN due to floating point
    # incurracy.
    term = numerator/denominator
    if np.isclose(term, 1, atol=1e-12):
        return 0.0
    if np.isclose(term, -1, atol=1e-12):
        return np.pi
    return np.arccos(term)


class XTBInvalidSolventError(Exception):
    ...


def is_valid_xtb_solvent(gfn_version, solvent):
    """
    Check if solvent is valid for the given GFN version.

    Parameters
    ----------
    gfn_version : :class:`int`
        GFN parameterization version. Can be: ``0``, ``1`` or ``2``.

    solvent : :class:`str`
        Solvent being tested [1]_.

    Returns
    -------
    :class:`bool`
        ``True`` if solvent is valid.

    References
    ----------
    .. [1] https://xtb-docs.readthedocs.io/en/latest/gbsa.html

    """
    if gfn_version == 0:
        return False
    elif gfn_version == 1:
        valid_solvents = {
            'acetone', 'acetonitrile', 'benzene',
            'CH2Cl2'.lower(), 'CHCl3'.lower(), 'CS2'.lower(),
            'DMSO'.lower(), 'ether', 'H2O'.lower(),
            'methanol', 'THF'.lower(), 'toluene'
        }
        return solvent in valid_solvents
    elif gfn_version == 2:
        valid_solvents = {
            'acetone', 'acetonitrile', 'CH2Cl2'.lower(),
            'CHCl3'.lower(), 'CS2'.lower(), 'DMF'.lower(),
            'DMSO'.lower(), 'ether', 'H2O'.lower(), 'methanol',
            'n-hexane'.lower(), 'THF'.lower(), 'toluene'
        }
        return solvent in valid_solvents


class XTBExtractor:
    """
    Extracts properties from GFN-xTB output files.

    Attributes
    ----------
    output_file : :class:`str`
        Output file to extract properties from.

    output_lines : :class:`list` : :class:`str`
        :class:`list` of all lines in as :class:`str` in the output
        file.

    total_energy : :class:`float`
        The total energy in the :attr:`output_file` as
        :class:`float`. The energy is in units of a.u..

    homo_lumo_gap : :class:`float`
        The HOMO-LUMO gap in the :attr:`output_file` as
        :class:`float`. The gap is in units of eV.

    fermi_level : :class:`float`
        The Fermi level in the :attr:`output_file` as
        :class:`float` in units of eV.

    qonly_dipole_moment : :class:`list`
        Components of the Q only dipole moment in units
        of Debye in :class:`list` of the form
        ``[x, y, z]``.

    full_dipole_moment : :class:`list`
        Components of the full dipole moment in units
        of Debye in :class:`list` of the form
        ``[x, y, z, total]``.

    qonly_quadrupole_moment : :class:`list`
        Components of the Q only traceless quadrupole moment in units
        of Debye in :class:`list` of the form
        ``[xx, xy, xy, xz, yz, zz]``.

    qdip_quadrupole_moment : :class:`list`
        Components of the Q+Dip traceless quadrupole moment in units of
        Debye in :class:`list` of the form
        ``[xx, xy, xy, xz, yz, zz]``.

    full_quadrupole_moment : :class:`list`
        Components of the full traceless quadrupole moment in units of
        Debye in :class:`list` of the form
        ``[xx, xy, xy, xz, yz, zz]``.

    homo_lumo_occ : :class:`dict`
        :class:`dict` of :class:`list` containing the orbital number,
        energy in eV and occupation of the HOMO and LUMO orbitals in
        the :attr:`output_file`.

    total_free_energy : :class:`float`
        The total free energy in the :attr:`output_file` as
        :class:`float`. The free energy is in units of a.u. and
        calculated at 298.15K.

    frequencies : :class:`list`
        :class:`list` of the vibrational frequencies as :class:`float`
        in the :attr:`output_file`. Vibrational frequencies are in
        units of wavenumber and calculated at 298.15K.

    """
    def __init__(self, output_file):
        """
        Initializes :class:`XTBExtractor`

        Parameters
        ----------
        output_file : :class:`str`
            Output file to extract properties from.

        """
        self.output_file = output_file
        # Explictly set encoding to UTF-8 because default encoding on
        # Windows will fail to read the file otherwise.
        with open(self.output_file, 'r', encoding='UTF-8') as f:
            self.output_lines = f.readlines()

        self._extract_values()

    def _extract_values(self):
        """
        Extract all properties from xTB output file.

        Returns
        -------
        None : :class:`NoneType`

        """

        for i, line in enumerate(self.output_lines):
            if self._check_line(line, 'total_energy'):
                self._extract_total_energy(line)
            elif self._check_line(line, 'homo_lumo_gap'):
                self._extract_homo_lumo_gap(line)
            elif self._check_line(line, 'fermi_level'):
                self._extract_fermi_level(line)
            elif self._check_line(line, 'dipole_moment'):
                self._extract_qonly_dipole_moment(i)
                self._extract_full_dipole_moment(i)
            elif self._check_line(line, 'quadrupole_moment'):
                self._extract_qonly_quadrupole_moment(i)
                self._extract_qdip_quadrupole_moment(i)
                self._extract_full_quadrupole_moment(i)
            elif self._check_line(line, 'homo_lumo_occ_HOMO'):
                self.homo_lumo_occ = {}
                self._extract_homo_lumo_occ(line, 'HOMO')
            elif self._check_line(line, 'homo_lumo_occ_LUMO'):
                self._extract_homo_lumo_occ(line, 'LUMO')
            elif self._check_line(line, 'total_free_energy'):
                self._extract_total_free_energy(line)

        # Frequency formatting requires loop through full file.
        self._extract_frequencies()

    def _check_line(self, line, option):
        """
        Checks a line for a string based on option.

        All formatting based on the 190418 version of xTB.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to check.

        option : :class:`str`
            Define which property and string being checked for.
            Can be one of ``'total_energy'``, ``'homo_lumo_gap'``,
            ``'fermi_level'``, ``'dipole_moment'``,
            ``'quadrupole_moment'``, ``'homo_lumo_occ_HOMO'``,
            ``'homo_lumo_occ_LUMO'``,
            ``'total_free_energy'``.

        Returns
        -------
        :class:`bool`
            Returns ``True`` if the desired string is present.

        """
        options = {
            'total_energy': '          | TOTAL ENERGY  ',
            'homo_lumo_gap': '          | HOMO-LUMO GAP   ',
            'fermi_level': '             Fermi-level        ',
            'dipole_moment': 'molecular dipole:',
            'quadrupole_moment': 'molecular quadrupole (traceless):',
            'homo_lumo_occ_HOMO': '(HOMO)',
            'homo_lumo_occ_LUMO': '(LUMO)',
            'total_free_energy': '          | TOTAL FREE ENERGY  ',
        }

        if options[option] in line:
            return True

    def _extract_total_energy(self, line):
        """
        Updates :attr:`total_energy`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.total_energy = float(string)

    def _extract_homo_lumo_gap(self, line):
        """
        Updates :attr:`homo_lumo_gap`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.homo_lumo_gap = float(string)

    def _extract_fermi_level(self, line):
        """
        Updates :attr:`fermi_level`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        part2 = line.split('Eh')
        string = nums.search(part2[1].rstrip()).group(0)
        self.fermi_level = float(string)

    def _extract_qonly_dipole_moment(self, index):
        """
        Updates :attr:`qonly_dipole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+2].rstrip()

        if 'q only:' in sample_set:
            self.qonly_dipole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_full_dipole_moment(self, index):
        """
        Updates :attr:`full_dipole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+3].rstrip()

        if 'full:' in sample_set:
            self.full_dipole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_qonly_quadrupole_moment(self, index):
        """
        Updates :attr:`qonly_quadrupole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+2].rstrip()

        if 'q only:' in sample_set:
            self.qonly_quadrupole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_qdip_quadrupole_moment(self, index):
        """
        Updates :attr:`qdip_quadrupole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+3].rstrip()

        if 'q+dip:' in sample_set:
            self.qdip_quadrupole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_full_quadrupole_moment(self, index):
        """
        Updates :attr:`full_quadrupole_moment`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        index : :class:`int`
            Index of line in :attr:`output_lines`.

        Returns
        -------
        None : :class:`NoneType`

        """

        sample_set = self.output_lines[index+4].rstrip()

        if 'full:' in sample_set:
            self.full_quadrupole_moment = [
                float(i)
                for i in sample_set.split(':')[1].split(' ') if i
            ]

    def _extract_homo_lumo_occ(self, line, orbital):
        """
        Updates :attr:`homo_lumo_occ`.

        Parameters
        ----------
        line : :class:`str`
            Line of output file to extract property from.

        orbital : :class:`str`
            Can be 'HOMO' or 'LUMO'.

        Returns
        -------
        None : :class:`NoneType`

        """

        if orbital == 'HOMO':
            split_line = [i for i in line.rstrip().split(' ') if i]
            # The line is:
            #   Number, occupation, energy (Ha), energy (ev), label
            # Extract:
            #   Number, occupation, energy (eV)
            orbital_val = [
                int(split_line[0]),
                float(split_line[1]),
                float(split_line[3])
            ]
        elif orbital == 'LUMO':
            split_line = [i for i in line.rstrip().split(' ') if i]
            # The line is:
            #   Number, energy (Ha), energy (ev), label
            # Extract:
            #   Number, occupation (zero), energy (eV)
            orbital_val = [
                int(split_line[0]),
                0,
                float(split_line[2])
            ]

        self.homo_lumo_occ[orbital] = orbital_val

    def _extract_total_free_energy(self, line):
        """
        Updates :attr:`total_free_energy`.

        Returns
        -------
        None : :class:`NoneType`

        """

        nums = re.compile(r"[+-]?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?")
        string = nums.search(line.rstrip()).group(0)
        self.total_free_energy = float(string)

    def _extract_frequencies(self):
        """
        Updates :attr:`frequencies`.

        Returns
        -------
        None : :class:`NoneType`

        """

        test = '|               Frequency Printout                |'

        # Use a switch to make sure we are extracting values after the
        # final property readout.
        switch = False

        frequencies = []
        for i, line in enumerate(self.output_lines):
            if test in line:
                # Turn on reading as final frequency printout has
                # begun.
                switch = True
            if ' reduced masses (amu)' in line:
                # Turn off reading as frequency section is done.
                switch = False
            if 'eigval :' in line and switch is True:
                samp = line.rstrip().split(':')[1].split(' ')
                split_line = [i for i in samp if i]
                for freq in split_line:
                    frequencies.append(freq)

        self.frequencies = [float(i) for i in frequencies]
