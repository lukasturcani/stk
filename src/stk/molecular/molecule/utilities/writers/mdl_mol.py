"""
MDL Writing Utilities
=====================

"""


def _write_mdl_mol_file(self, path, atom_ids):
    """
    Write to a V3000 ``.mol`` file.

    This function should not be used directly, only via
    :meth:`write`.

    Parameters
    ----------
    path : :class:`str`
        The full path to the file being written.

    atom_ids : :class:`iterable` of :class:`int`
        The atom ids of atoms to write. Can be a single
        :class:`int`, if a single atom is to be used, or ``None``,
        if all atoms are to be used.

    Returns
    -------
    None : :class:`NoneType`

    """

    with open(path, 'w') as f:
        f.write(_to_mdl_mol_block(self, atom_ids))


def _to_mdl_mol_block(self, atom_ids=None):
    """
    Return a V3000 mol block of the molecule.

    Parameters
    ----------
    atom_ids : :class:`iterable` of :class:`int`, optional
        The atom ids of atoms to write. Can be a single
        :class:`int`, if a single atom is to be used, or
        ``None``, if all atoms are to be used.

    Returns
    -------
    :class:`str`
        The V3000 mol block representing the molecule.

    """

    if atom_ids is None:
        atom_ids = range(len(self._atoms))
    elif isinstance(atom_ids, int):
        atom_ids = (atom_ids, )

    atom_lines = []
    # This set gets used by bonds.
    id_map = {}
    for new_atom_id, old_atom_id in enumerate(atom_ids, 1):
        id_map[old_atom_id] = new_atom_id

        x, y, z = self._position_matrix[:, old_atom_id]
        atom = self._atoms[old_atom_id]
        symbol = atom.__class__.__name__
        charge = atom.get_charge()
        charge = f' CHG={charge}' if charge else ''
        atom_lines.append(
            f'M  V30 {new_atom_id} {symbol} {x:.4f} '
            f'{y:.4f} {z:.4f} 0{charge}\n'
        )
    atom_block = ''.join(atom_lines)

    bond_lines = []
    for bond in self._bonds:
        a1 = bond.get_atom1().get_id()
        a2 = bond.get_atom2().get_id()
        if a1 in id_map and a2 in id_map:
            bond_lines.append(
                f'M  V30 {len(bond_lines)+1} '
                f'{int(bond.get_order())} '
                f'{id_map[a1]} {id_map[a2]}\n'
            )

    num_bonds = len(bond_lines)
    bond_block = ''.join(bond_lines)
    return (
        '\n'
        '     RDKit          3D\n'
        '\n'
        '  0  0  0  0  0  0  0  0  0  0999 V3000\n'
        'M  V30 BEGIN CTAB\n'
        f'M  V30 COUNTS {len(id_map)} {num_bonds} 0 0 0\n'
        'M  V30 BEGIN ATOM\n'
        f'{atom_block}'
        'M  V30 END ATOM\n'
        'M  V30 BEGIN BOND\n'
        f'{bond_block}'
        'M  V30 END BOND\n'
        'M  V30 END CTAB\n'
        'M  END\n'
        '\n'
        '$$$$\n'
    )
