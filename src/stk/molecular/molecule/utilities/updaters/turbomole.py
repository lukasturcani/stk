"""
Turbomole Updating Utilities
============================

"""

from __future__ import annotations

import typing
import pathlib
import numpy as np
from collections import abc
from stk.utilities import periodic_table


class _CoordSection:
    """
    Holds the coord section of a turbomole file.

    """

    def __init__(
        self,
        lines: abc.Iterable[str],
    ):
        """
        Initialize a :class:`._CoordSection`.

        This initializer assumes that coordinates used in `lines` are
        defined in angstroms.

        Parameters:
            lines:
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
    def init_bohr(
        cls,
        lines: abc.Iterable[str],
    ) -> _CoordSection:
        """
        Initialize a :class:`._CoordSection`.

        This initializer assumes that coordinates used in `lines`
        are defined in bohr units.

        Parameters:

            lines:
                The lines of the coord section.

        Returns:

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

    def get_position_matrix(self) -> np.ndarray:
        """
        Get the position matrix defined in the coord section.

        Returns:

            The position matrix defined in the coord section.

        """

        return np.array(self._position_matrix)

    def get_num_atoms(self) -> int:
        """
        Get the number of atoms defined in the coord section.

        Returns:

            The number of atoms defined in the coord section.

        """

        return len(self._elements)

    def get_elements(self) -> abc.Iterable[str]:
        """
        Yield the elements of the atoms in the coord section.

        The elements are yielded in order of atom id.

        Yields:

            The chemical symbol of an atom.

        """

        yield from self._elements


def _get_coord_section(
    path: typing.Union[pathlib.Path, str],
    num_atoms: int,
) -> _CoordSection:
    """
    Get the coord section defined in `path`.

    Parameters:

        path:
            The path to the turbomole file.

        num_atoms:
            The number of atoms in the molecule.

    Returns:

        The coord section defined in `path`.

    Raises:

        :class:`RuntimeError`:
            If `path` uses fractional coordinates, as they are not
            currently supported.

        :class:`RuntimeError`:
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


def get_position_matrix_from_turbomole(
    num_atoms: int,
    path: typing.Union[pathlib.Path, str],
) -> np.ndarray:
    """
    Get the position matrix from a Turbomole file.

    Note that coordinates in ``.coord`` files can be given in Bohr or
    Angstrom, which is handled. Fractional coordinates are not
    currently handled.

    Parameters:

        num_atoms:
            The number of atoms in the molecule.

        path:
            The full path to the ``.coord`` file which holds the
            position matrix.

    Returns:

        The position matrix.

    """

    section = _get_coord_section(path, num_atoms)
    return section.get_position_matrix()
