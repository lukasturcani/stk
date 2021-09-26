"""
XYZ Writing Utilities
=====================

"""


def _write_xyz_file(self, path, atom_ids):
    """
    Write to a ``.xyz`` file.

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

    content = [0]
    for i, atom_id in enumerate(atom_ids, 1):
        x, y, z = self._position_matrix[:, atom_id]
        symbol = self._atoms[atom_id].__class__.__name__
        content.append(f'{symbol} {x:f} {y:f} {z:f}\n')
    # Set first line to the atom_count.
    content[0] = f'{i}\n\n'

    with open(path, 'w') as xyz:
        xyz.write(''.join(content))
