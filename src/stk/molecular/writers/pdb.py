"""
PDB Writer
==========

"""


class PdbWriter:
    """
    A writer class for ``.pdb`` files.

    Examples
    --------
    *Writing to a File with a Unit Cell*

    This writer can write to a file with the unit
    cell included for periodic molecules. Note that this always assumes
    P1 space group.

    .. testcode:: writing-to-a-file-with-a-unit-cell

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
        topology_graph = stk.cof.PeriodicHoneycomb(
            building_blocks=(bb1, bb2),
            lattice_size=(3, 3, 1),
        )
        construction_result = topology_graph.construct()
        cof = stk.ConstructedMolecule.init_from_construction_result(
            construction_result=construction_result,
        )
        writer = stk.PdbWriter()
        writer.write(
            molecule=cof,
            path='cof.pdb',
            periodic_info=construction_result.get_periodic_info(),
        )

    .. testcode:: writing-to-a-file-with-a-unit-cell
        :hide:

        import os

        assert os.path.exists('cof.pdb')

    .. testcleanup:: writing-to-a-file-with-a-unit-cell

        os.remove('cof.pdb')

    """

    def _write_content(self, molecule, atom_ids, periodic_info=None):

        if atom_ids is None:
            atom_ids = range(molecule.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        content = []
        if periodic_info is not None:
            # Input unit cell information.
            a = periodic_info.get_a()
            b = periodic_info.get_b()
            c = periodic_info.get_c()
            alpha = periodic_info.get_alpha()
            beta = periodic_info.get_beta()
            gamma = periodic_info.get_gamma()
            content.append(
                f'CRYST1 {a:>8.3f} {b:>8.3f} {c:>8.3f}'
                f' {alpha:>6.2f} {beta:>6.2f} {gamma:>6.2f} '
                f'P 1\n'
            )

        atom_counts = {}
        hetatm = 'HETATM'
        alt_loc = ''
        res_name = 'UNL'
        chain_id = ''
        res_seq = '1'
        i_code = ''
        occupancy = '1.00'
        temp_factor = '0.00'

        coords = molecule.get_position_matrix()
        # This set will be used by bonds.
        atoms = set()
        for atom_id in atom_ids:
            atoms.add(atom_id)
            atom, = molecule.get_atoms(atom_ids=atom_id)
            serial = atom_id+1
            element = atom.__class__.__name__
            charge = atom.get_charge()
            atom_counts[element] = atom_counts.get(element, 0) + 1
            name = f'{element}{atom_counts[element]}'
            # Make sure the coords are no more than 8 columns wide
            # each.
            x, y, z = (i for i in coords[atom_id])

            content.append(
                f'{hetatm:<6}{serial:>5} {name:<4}'
                f'{alt_loc:<1}{res_name:<3} {chain_id:<1}'
                f'{res_seq:>4}{i_code:<1}   '
                f' {x:>7.3f} {y:>7.3f} {z:>7.3f}'
                f'{occupancy:>6}{temp_factor:>6}          '
                f'{element:>2}{charge:>2}\n'
            )

        conect = 'CONECT'
        for bond in molecule.get_bonds():
            a1 = bond.get_atom1().get_id()
            a2 = bond.get_atom2().get_id()
            if a1 in atoms and a2 in atoms:
                content.append(
                    f'{conect:<6}{a1+1:>5}{a2+1:>5}               \n'
                )

        content.append('END\n')

        return content

    def to_string(
        self,
        molecule,
        atom_ids=None,
        periodic_info=None
    ):
        """
        Get a ``.pdb`` file format string of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to ``.pdb`` format.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        periodic_info : :class:`.PeriodicInfo`
            Information about the periodic cell.

        Returns
        -------
        :class:`string`
            A string holding the content of a ``.pdf`` file.

        """

        content = self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_info=periodic_info,
        )

        return ''.join(content)

    def write(
        self,
        molecule,
        path,
        atom_ids=None,
        periodic_info=None
    ):
        """
        Write `molecule` to ``.pdb`` file format.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to ``.pdb`` format.

        path : :class:`str`
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        periodic_info : :class:`.PeriodicInfo`
            Information about the periodic cell.

        Returns
        -------
        None : :class:`NoneType`
            A file is written.

        """

        content = self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_info=periodic_info,
        )

        with open(path, 'w') as f:
            f.write(''.join(content))
