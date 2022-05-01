"""
XYZ Writer
==========

"""

from __future__ import annotations

import typing

from ...utilities import OneOrMany
from ..molecules import Molecule


class XyzWriter:
    """
    A writer class for ``.xyz`` files.

    Examples:

        *Writing to a File*

        .. testcode:: writing-to-a-file

            import stk

            bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])

            writer = stk.XyzWriter()
            writer.write(molecule=bb1, path='bb1.xyz')

        .. testcode:: writing-to-a-file
            :hide:

            import os

            assert os.path.exists('bb1.xyz')

        .. testcleanup:: writing-to-a-file

            os.remove('bb1.xyz')

    """

    def _write_content(
        self,
        molecule: Molecule,
        atom_ids: typing.Optional[OneOrMany[int]],
    ) -> list[str]:

        if atom_ids is None:
            atom_ids = range(molecule.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        coords = molecule.get_position_matrix()
        content = ['']
        for i, atom_id in enumerate(atom_ids, 1):
            x, y, z = (i for i in coords[atom_id])
            atom, = molecule.get_atoms(atom_ids=atom_id)
            symbol = atom.__class__.__name__
            content.append(f'{symbol} {x:f} {y:f} {z:f}\n')
        # Set first line to the atom_count.
        content[0] = f'{i}\n\n'

        return content

    def to_string(
        self,
        molecule: Molecule,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> str:
        """
        Get the ``.xyz`` string of  `molecule`.

        Parameters:

            molecule:
                Molecule to write to `.xyz` format.

            atom_ids:
                The atom ids of atoms to write. Can be a single
                :class:`int`, if a single atom is to be used, or
                ``None``, if all atoms are to be used.

        Returns:

            String in ``.xyz`` file format.

        """

        content = self._write_content(molecule, atom_ids)

        return ''.join(content)

    def write(
        self,
        molecule: Molecule,
        path: str,
        atom_ids: typing.Optional[OneOrMany[int]] = None,
    ) -> None:
        """
        Write `molecule` to ``.xyz`` file format.

        Parameters:

            molecule:
                Molecule to write to ``.xyz`` format.

            path:
                The full path to the file being written.

            atom_ids:
                The atom ids of atoms to write. Can be a single
                :class:`int`, if a single atom is to be used, or
                ``None``, if all atoms are to be used.

        """

        content = self._write_content(molecule, atom_ids)

        with open(path, 'w') as f:
            f.write(''.join(content))
