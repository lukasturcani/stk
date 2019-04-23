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

# Holds global stk options.
OPTIONS = {}

# Holds the elements Van der Waals radii in Angstroms.
atom_vdw_radii = {
              'Al': 2, 'Sb': 2, 'Ar': 1.88, 'As': 1.85, 'Ba': 2,
              'Be': 2, 'Bi': 2, 'B': 2, 'Br': 1.85, 'Cd': 1.58,
              'Cs': 2, 'Ca': 2, 'C': 1.7, 'Ce': 2, 'Cl': 1.75,
              'Cr': 2, 'Co': 2, 'Cu': 1.4, 'Dy': 2, 'Er': 2,
              'Eu': 2, 'F':  1.47, 'Gd': 2, 'Ga': 1.87, 'Ge': 2,
              'Au': 1.66, 'Hf': 2, 'He': 1.4, 'Ho': 2, 'H': 1.09,
              'In': 1.93, 'I': 1.98, 'Ir': 2, 'Fe': 2, 'Kr': 2.02,
              'La': 2, 'Pb': 2.02, 'Li': 1.82, 'Lu': 2, 'Mg': 1.73,
              'Mn': 2, 'Hg': 1.55, 'Mo': 2, 'Nd': 2, 'Ne': 1.54,
              'Ni': 1.63, 'Nb': 2, 'N':  1.55, 'Os': 2, 'O':  1.52,
              'Pd': 1.63, 'P': 1.8, 'Pt': 1.72, 'K': 2.75, 'Pr': 2,
              'Pa': 2, 'Re': 2, 'Rh': 2, 'Rb': 2, 'Ru': 2, 'Sm': 2,
              'Sc': 2, 'Se': 1.9, 'Si': 2.1, 'Ag': 1.72, 'Na': 2.27,
              'Sr': 2, 'S': 1.8, 'Ta': 2, 'Te': 2.06, 'Tb': 2,
              'Tl': 1.96, 'Th': 2, 'Tm': 2, 'Sn': 2.17, 'Ti': 2,
              'W': 2, 'U':  1.86, 'V':  2, 'Xe': 2.16, 'Yb': 2,
              'Y': 2, 'Zn': 1.29, 'Zr': 2, 'X':  1.0, 'D':  1.0
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
              7: 'N', 8: 'O',  9: 'F', 10: 'Ne', 11: 'Na', 12: 'Mg',
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
              117: 'Uus', 118: 'Uuo'}


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

    Initializing this class with a MacroMolecule finds that
    MacroMolecules -out.maegz file and converts it to a .mae file. It
    then creates and additional .mae file holding only the lowest
    energy conformer found.

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

    def __init__(self, file, n=1):
        name, ext = os.path.splitext(file)
        self.maegz_path = name + '-out.maegz'
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
                                            f'EXTRACTED_{num}.mae')
            else:
                new_name = self.mae_path.replace(
                              '.mae',
                              f'EXTRACTED_{num}_conf_{i}.mae')

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


class AtomicPeriodicBond:
    """
    Represents a periodic bond between two atoms.

    Attributes
    ----------
    bonder1 : :class:`int`
        The id of the first atom involved in the bond.

    bonder2 : :class:`int`
        The id of the second atom involved in the bond.

    bond_type : :class:`rdkit.BondType`
        The bond type.

    direction : :class:`list` of :class:`int`
        A 3 member list describing the axes along which the bond is
        periodic, when going from :attr:`bonder1` to :attr:`bonder2`.
        For example, ``[1, 0, 0]`` means that the bond is periodic
        along the x axis in the positive direction.

    """

    def __init__(self, bonder1, bonder2, bond_type, direction):
        self.bonder1 = bonder1
        self.bonder2 = bonder2
        self.bond_type = bond_type
        self.direction = direction


class PeriodicBond:
    """
    Represents a periodic bond.

    Parameters
    ----------
    fg1 : :class:`.FunctionalGroup`
        The first functional group involved in the bond.

    fg2 : :class:`.FunctionalGroup`
        The second functional group involved in the bond.

    direction : :class:`list` of :class:`int`
        A 3 member list describing the axes along which the bond is
        periodic, when going from :attr:`fg1` to :attr:`fg2`. For
        example, ``[1, 0, 0]`` means that the bond is periodic along
        the x axis in the positive direction.

    Attributes
    ----------
    fg1 : :class:`.FunctionalGroup`
        The first functional group involved in the bond.

    fg2 : :class:`.FunctionalGroup`
        The second functional group involved in the bond.

    direction : :class:`numpy.ndarray` of :class:`int`
        A 3 member list describing the axes along which the bond is
        periodic, when going from `fg1` to `fg2`. For example,
        ``[1, 0, 0]`` means that the bond is periodic along the x axis
        in the positive direction.

    """

    def __init__(self, fg1, fg2, direction):
        self.fg1 = fg1
        self.fg2 = fg2
        self.direction = np.array(direction)

    def __str__(self):
        return (
            f"PeriodicBond({self.fg1}, {self.fg2}, {self.direction})"
        )


class StopLogging:
    ...


def add_fragment_props(mol, bb_index, mol_index):
    """
    Adds properties to `mol` of `bb_index` and `mol_index`.

    Properties called 'bb_index' and 'mol_index' are added to every
    atom in `mol`.

    Parameters
    ----------
    mol : :class:`rdkit.Mol`
        A molecule which needs to have its atoms tagged.

    bb_index : :class:`int`
        The index of `mol` within `building_blocks` of some
        macromolecule.

    mol_index : :class:`int`
        If `mol` is the 5th molecule of building block with `bb_index`
        to be added to a macromolecule, `mol_index` is 4.

    Returns
    -------
    None : NoneType

    """

    for atom in mol.GetAtoms():
        atom.SetIntProp('bb_index', bb_index)
        atom.SetIntProp('mol_index', mol_index)


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


def centroid(*coords):
    """
    Calculates the centroid of a group of coordinates.

    Parameters
    ----------
    *coords : :class:`numpy.ndarray`
        Any number of numpy arrays holding x, y and z positions.

    Returns
    -------
    :class:`numpy.array`
        The centroid of the coordinates `coords`.

    """

    return sum(coords) / len(coords)


def dedupe(iterable, seen=None, key=None):
    """
    Yields items from `iterable` barring duplicates.

    If `seen` is provided it contains elements which are not to be
    yielded at all.

    Parameters
    ----------
    iterable : :class:`iterable`
        An iterable of elements which are to be yielded, only once.

    seen : :class:`set`, optional
        Holds items which are not to be yielded.

    key : :class:`callable`
        A function which gets applied to every member of `iterable`.
        The return of this function is checked for duplication rather
        than the member itself.

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
            yield from flatten(x)
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
                *_, bond_id,  bond_order, atom1, atom2 = line.split()
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

    v = np.divide(vector, np.linalg.norm(vector))
    return np.round(v, decimals=4)


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
    for a in mol.GetAtoms():
        new_atom = rdkit.Atom(a.GetAtomicNum())
        # Set properties.
        for pname, pval in a.GetPropsAsDict(False, False).items():
            if isinstance(pval, int):
                new_atom.SetIntProp(pname, pval)
            elif isinstance(pval, bool):
                new_atom.SetBoolProp(pname, pval)
            else:
                new_atom.SetProp(pname, pval)

        new_atom.SetFormalCharge(a.GetFormalCharge())
        emol.AddAtom(new_atom)

    for bond in mol.GetBonds():
        emol.AddBond(bond.GetBeginAtomIdx(),
                     bond.GetEndAtomIdx(),
                     bond.GetBondType())

    m = emol.GetMol()
    m.AddConformer(rdkit.Conformer(mol.GetConformer()))

    for a in m.GetAtoms():
        a.UpdatePropertyCache()

    rdkit.GetSSSR(m)

    return m


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
    if np.array_equal(vector1, vector2):
        return np.identity(3)

    # Handle the case where the rotation is 180 degrees.
    if np.array_equal(vector1, np.multiply(vector2, -1)):
        # Get a vector orthogonal to `vector1` by finding the smallest
        # component of `vector1` and making that a vector.
        ortho = [0, 0, 0]
        ortho[list(vector1).index(min(abs(vector1)))] = 1
        axis = np.cross(vector1, ortho)
        return rotation_matrix_arbitrary_axis(np.pi, axis)

    v = np.cross(vector1, vector2)

    vx = np.array([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])

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


def vector_theta(vector1, vector2):
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
    denominator = (np.linalg.norm(vector1) *
                   np.linalg.norm(vector2))
    # This if statement prevents returns of NaN due to floating point
    # incurracy.
    if np.isclose(numerator, denominator, atol=1e-8):
        return 0.0
    return np.arccos(numerator/denominator)
