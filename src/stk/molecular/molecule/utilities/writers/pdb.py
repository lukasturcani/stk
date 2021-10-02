"""
PDB Writing Utilities
=====================

"""

from __future__ import annotations
import pathlib
import typing
import numpy as np

from stk.utilities import typing as _typing
from .... import atoms as _atoms
from .... import bonds as _bonds


def write_pdb_file(
    atoms: tuple[_atoms.Atom, ...],
    bonds: tuple[_bonds.Bond, ...],
    position_matrix: np.ndarray,
    path: typing.Union[pathlib.Path, str],
    atom_ids: typing.Optional[_typing.OneOrMany[int]],
) -> None:
    """
    Write to a ``.pdb`` file.

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

    lines = []
    atom_counts: dict[str, int] = {}
    hetatm = 'HETATM'
    alt_loc = ''
    res_name = 'UNL'
    chain_id = ''
    res_seq = '1'
    i_code = ''
    occupancy = '1.00'
    temp_factor = '0.00'

    coords = position_matrix
    # This set will be used by bonds.
    seen_atoms: set[int] = set()
    for atom_id in atom_ids:
        seen_atoms.add(atom_id)

        serial = atom_id+1
        element = atoms[atom_id].__class__.__name__
        atom_counts[element] = atom_counts.get(element, 0) + 1
        name = f'{element}{atom_counts[element]}'
        # Make sure the coords are no more than 8 columns wide
        # each.
        x, y, z = (i for i in coords[:, atom_id])
        lines.append(
            f'{hetatm:<6}{serial:>5} {name:<4}'
            f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
            f'{res_seq:>4}{i_code:<1}   '
            f' {x:>7.3f} {y:>7.3f} {z:>7.3f}'
            f'{occupancy:>6}{temp_factor:>6}          '
            f'{element:>2}{atoms[atom_id].get_charge():>2}\n'
        )

    conect = 'CONECT'
    for bond_ in bonds:
        a1 = bond_.get_atom1().get_id()
        a2 = bond_.get_atom2().get_id()
        if a1 in seen_atoms and a2 in seen_atoms:
            lines.append(
                f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
            )

    lines.append('END\n')
    with open(path, 'w') as f:
        f.write(''.join(lines))
