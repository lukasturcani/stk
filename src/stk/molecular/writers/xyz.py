"""
XYZ Writer
==========

"""


class XyzWriter:
    """
    A writer class for ``.xyz`` files.

    Examples
    --------
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

    def _write_content(self, molecule, atom_ids):

        if atom_ids is None:
            atom_ids = range(molecule.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        coords = molecule.get_position_matrix()
        content = [0]
        for i, atom_id in enumerate(atom_ids, 1):
            x, y, z = (i for i in coords[atom_id])
            atom, = molecule.get_atoms(atom_ids=atom_id)
            symbol = atom.__class__.__name__
            content.append(f'{symbol} {x:f} {y:f} {z:f}\n')
        # Set first line to the atom_count.
        content[0] = f'{i}\n\n'

        return content

    def to_string(self, molecule, atom_ids=None):
        """
        Get the ``.xyz`` string of  `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to `.xyz` format.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        Returns
        -------
        :class:`str`
            String in ``.xyz`` file format.

        """

        content = self._write_content(molecule, atom_ids)

        return ''.join(content)

    def write(self, molecule, path, atom_ids=None):
        """
        Write `molecule` to ``.xyz`` file format.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to ``.xyz`` format.

        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        Returns
        -------
        None : :class:`NoneType`
            A file is written.

        """

        content = self._write_content(molecule, atom_ids)

        with open(path, 'w') as f:
            f.write(''.join(content))
