"""
Mol Writer
==========

"""


class MolWriter:
    """
    A writer class for V3000 ``.mol`` files.

    Examples
    --------
    *Writing to a File*

    .. testcode:: writing-to-a-file

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])

        writer = stk.MolWriter()
        writer.write(molecule=bb1, path='bb1.mol')

    .. testcode:: writing-to-a-file
        :hide:

        import os

        assert os.path.exists('bb1.mol')

    .. testcleanup:: writing-to-a-file

        os.remove('bb1.mol')

    """

    def _write_content(self, molecule, atom_ids):

        if atom_ids is None:
            atom_ids = range(molecule.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        atom_lines = []
        coords = molecule.get_position_matrix()
        # This set gets used by bonds.
        id_map = {}
        for new_atom_id, old_atom_id in enumerate(atom_ids, 1):
            id_map[old_atom_id] = new_atom_id

            x, y, z = (i for i in coords[old_atom_id])
            atom, = molecule.get_atoms(atom_ids=old_atom_id)
            symbol = atom.__class__.__name__
            charge = atom.get_charge()
            charge = f' CHG={charge}' if charge else ''
            atom_lines.append(
                f'M  V30 {new_atom_id} {symbol} {x:.4f} '
                f'{y:.4f} {z:.4f} 0{charge}\n'
            )
        atom_block = ''.join(atom_lines)

        bond_lines = []
        for bond in molecule.get_bonds():
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

    def to_string(self, molecule, atom_ids=None):
        """
        Get a V3000 ``.mol`` file format string of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to V3000 ``.mol`` format.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        Returns
        -------
        :class:`str`
            String in V3000 ``.mol`` file format.

        """

        return self._write_content(molecule, atom_ids)

    def write(self, molecule, path, atom_ids=None):
        """
        Write `molecule` to V3000 ``.mol`` file format.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to V3000 ``.mol`` format.

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

        with open(path, 'w') as f:
            f.write(self._write_content(molecule, atom_ids))
