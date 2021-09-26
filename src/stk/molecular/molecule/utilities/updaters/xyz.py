"""
XYZ Updating Utilities
======================

"""

import numpy as np

from stk.utilities import periodic_table


def _with_structure_from_xyz(self, path):
    """
    Return a clone, with its structure taken from an ``.xyz`` file.

    Parameters
    ----------
    path : :class:`str`
        The full path of the ``.mol`` file from which the structure
        should be updated.

    Returns
    -------
    :class:`.Molecule`
        A clone with atomic positions found in `path`.

    Raises
    ------
    :class:`RuntimeError`
        If the number of atoms in the file does not match the
        number of atoms in the molecule or if atom elements in the
        file do not agree with the atom elements in the molecule.

    """

    with open(path, 'r') as f:
        atom_count, _, *content = f.readlines()

    # Check the atom count is correct.
    num_atoms = len(self._atoms)
    if int(atom_count) != num_atoms:
        raise RuntimeError(
            f'The number of atoms in the xyz file, {atom_count}, '
            'does not match the number of atoms in the '
            f'molecule, {num_atoms}.'
        )

    # Save all the coords in the file.
    new_coords = []
    for i, line in enumerate(content):
        element, *coords = line.split()
        # Handle XYZ files with capitilisation of element symbols.
        element = element.title()
        if element.isnumeric():
            element = periodic_table[int(element)]

        if element != self._atoms[i].__class__.__name__:
            raise RuntimeError(
                f'Atom {i} element does not match file.'
            )

        new_coords.append([float(i) for i in coords])

    # Check that the correct number of atom
    # lines was present in the file.
    if i+1 != num_atoms:
        raise RuntimeError(
            f'The number of atom lines in the xyz file, {i+1}, '
            'does not match the number of atoms in the '
            f'molecule, {num_atoms}.'
        )

    # Update the structure.
    new_coords = np.array(new_coords)
    return self._with_position_matrix(new_coords)
