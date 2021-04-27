"""
PDB Writing Utilities
=====================

"""


def _write_pdb_file(self, path, atom_ids):
    """
    Write to a ``.pdb`` file.

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

    if atom_ids is None:
        atom_ids = range(len(self._atoms))
    elif isinstance(atom_ids, int):
        atom_ids = (atom_ids, )

    lines = []
    atom_counts = {}
    hetatm = 'HETATM'
    alt_loc = ''
    res_name = 'UNL'
    chain_id = ''
    res_seq = '1'
    i_code = ''
    occupancy = '1.00'
    temp_factor = '0.00'

    coords = self._position_matrix
    # This set will be used by bonds.
    atoms = set()
    for atom in atom_ids:
        atoms.add(atom)

        serial = atom+1
        element = self._atoms[atom].__class__.__name__
        atom_counts[element] = atom_counts.get(element, 0) + 1
        name = f'{element}{atom_counts[element]}'
        # Make sure the coords are no more than 8 columns wide
        # each.
        x, y, z = (i for i in coords[:, atom])
        lines.append(
            f'{hetatm:<6}{serial:>5} {name:<4}'
            f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
            f'{res_seq:>4}{i_code:<1}   '
            f' {x:>7.3f} {y:>7.3f} {z:>7.3f}'
            f'{occupancy:>6}{temp_factor:>6}          '
            f'{element:>2}{self._atoms[atom].get_charge():>2}\n'
        )

    conect = 'CONECT'
    for bond in self._bonds:
        a1 = bond.get_atom1().get_id()
        a2 = bond.get_atom2().get_id()
        if a1 in atoms and a2 in atoms:
            lines.append(
                f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
            )

    lines.append('END\n')
    with open(path, 'w') as f:
        f.write(''.join(lines))
