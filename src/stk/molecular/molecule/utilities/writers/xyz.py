"""
XYZ Writing Utilities
=====================

"""

from __future__ import annotations

import typing
import pathlib
import numpy as np

from stk.utilities import typing as _typing
from ....atoms import Atom
from ....bonds import Bond


def write_xyz_file(
    atoms: tuple[Atom, ...],
    bonds: tuple[Bond, ...],
    position_matrix: np.ndarray,
    path: typing.Union[pathlib.Path, str],
    atom_ids: typing.Optional[_typing.OneOrMany[int]],
) -> None:
    """
    Write to a ``.xyz`` file.

    This function should not be used directly, only via
    :meth:`write`.

    Parameters:

        atoms:
            The atoms of the molecule to write.

        bonds:
            The bonds of the molecule to write.

        position_matrix:
            The ``3 x N`` position of the molecule to write.

        path:
            The full path to the file being written.

        atom_ids:
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

    """

    if atom_ids is None:
        atom_ids = range(len(atoms))
    elif isinstance(atom_ids, int):
        atom_ids = (atom_ids, )

    content = ['']
    for i, atom_id in enumerate(atom_ids, 1):
        x, y, z = position_matrix[:, atom_id]
        symbol = atoms[atom_id].__class__.__name__
        content.append(f'{symbol} {x:f} {y:f} {z:f}\n')
    # Set first line to the atom_count.
    content[0] = f'{i}\n\n'

    with open(path, 'w') as xyz:
        xyz.write(''.join(content))
