"""
Turbomole Updating Utilities
============================

"""

import numpy as np

from stk.utilities import periodic_table


class _CoordSection:
    """
    Holds the coord section of a turbomole file.

    """

    def __init__(self, lines):
        """
        Initialize a :class:`._CoordSection`.

        This initializer assumes that coordinates used in `lines` are
        defined in angstroms.

        Parameters
        ----------
        lines : :class:`iterable` of :class:`str`
            The lines of the coord section.

        """

        elements = []
        position_matrix = []
        for line in lines:
            x, y, z, element = line.split()

            if element.isnumeric():
                element = periodic_table[int(element)]

            elements.append(element)
            position_matrix.append([float(x), float(y), float(z)])

        self._position_matrix = np.array(position_matrix)
        self._elements = tuple(elements)

    @classmethod
    def init_bohr(cls, lines):
        """
        Initialize a :class:`._CoordSection`.

        This initializer assumes that coordinates used in `lines`
        are defined in bohr units.

        Parameters
        ----------
        lines : :class:`iterable` of :class:`str`
            The lines of the coord section.

        Returns
        -------
        :class:`._CoordSection`
            The coord section.

        """

        bohr_to_ang = 0.5291772105638411

        obj = cls.__new__(cls)
        elements = []
        position_matrix = []
        for line in lines:
            x, y, z, element = line.split()

            if element.isnumeric():
                element = periodic_table[int(element)]

            elements.append(element)
            position_matrix.append([
                float(x)*bohr_to_ang,
                float(y)*bohr_to_ang,
                float(z)*bohr_to_ang,
            ])

        obj._elements = tuple(elements)
        obj._position_matrix = np.array(position_matrix)
        return obj

    def get_position_matrix(self):
        """
        Get the position matrix defined in the coord section.

        Returns
        -------
        :class:`numpy.ndarray`
            The position matrix defined in the coord section.

        """

        return np.array(self._position_matrix)

    def get_num_atoms(self):
        """
        Get the number of atoms defined in the coord section.

        Returns
        -------
        :class:`int`
            The number of atoms defined in the coord section.

        """

        return len(self._elements)

    def get_elements(self):
        """
        Yield the elements of the atoms in the coord section.

        The elements are yielded in order of atom id.

        Yields
        ------
        :class:`str`
            The chemical symbol of an atom.

        """

        yield from self._elements


def _get_coord_section(path, num_atoms):
    """
    Get the coord section defined in `path`.

    Parameters
    ----------
    path : :class:`str`
        The path to the turbomole file.

    num_atoms : :class:`int`
        The number of atoms in the molecule.

    Returns
    -------
    :class:`._CoordSection`
        The coord section defined in `path`.

    Raises
    ------
    :class:`RuntimeError`
        If `path` uses fractional coordinates, as they are not
        currently supported.

    :class:`RuntimeError`
        If no coord section in found in `path`.

    """

    with open(path, 'r') as f:
        content = f.readlines()

    for line_number, line in enumerate(content):
        if '$coord' in line:

            lines = content[line_number+1:line_number+1+num_atoms]

            if 'angs' in line:
                return _CoordSection(lines)

            elif 'frac' in line:
                raise RuntimeError(
                    f'{path} uses fractional coordinates, which are '
                    'not currently supported.'
                )

            return _CoordSection.init_bohr(lines)

    raise RuntimeError(f'No coord section found in {path}.')


def _with_structure_from_turbomole(self, path):
    """
    Update the structure with one taken from a Turbomole file.

    Note that coordinates in ``.coord`` files can be given in Bohr or
    Angstrom, which is handled. Fractional coordinates are not
    currently handled.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.coord`` file from which the
        structure should be updated.

    Returns
    -------
    :class:`.Molecule`
        The molecule is returned.

    Raises
    ------
    :class:`RuntimeError`
        If the number of atoms in the file does not match the
        number of atoms in the molecule.

    :class:`RuntimeError`
        If atom elements in the file do not agree with the atom
        elements in the molecule.

    """

    num_atoms = len(self._atoms)
    section = _get_coord_section(path, num_atoms)

    if section.get_num_atoms() != num_atoms:
        raise RuntimeError(
            f'The number of atoms in {path}, '
            f'{section.get_num_atoms()}, does not match the number '
            f'of atoms in the molecule, {num_atoms}.'
        )

    for atom_id, (element, atom) in enumerate(zip(
        section.get_elements(),
        self._atoms,
    )):
        if element != atom.__class__.__name__:
            raise RuntimeError(
                f'The element of atom {atom_id} in {path}, '
                f'{element}, does not match the element in the '
                f'molecule, {atom.__class__.__name__}.'
            )

    return self._with_position_matrix(section.get_position_matrix())
