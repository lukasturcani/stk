"""
Turbomole Writer
================

"""


class TurbomoleWriter:
    """
    A writer class for ``Turbomole`` files.

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
        writer = stk.TurbomoleWriter()
        writer.write(
            molecule=cof,
            path='cof.coord',
            periodic_info=construction_result.get_periodic_info(),
        )

    .. testcode:: writing-to-a-file-with-a-unit-cell
        :hide:

        import os

        assert os.path.exists('cof.coord')

    .. testcleanup:: writing-to-a-file-with-a-unit-cell

        os.remove('cof.coord')

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
                '$periodic 3\n'
                '$cell angs\n'
                f' {a:>8.3f} {b:>8.3f} {c:>8.3f} '
                f'{alpha:>6.2f} {beta:>6.2f} {gamma:>6.2f}\n'
            )

        coords = molecule.get_position_matrix()
        content.append('$coord angs\n')
        for atom_id in atom_ids:
            atom, = molecule.get_atoms(atom_ids=atom_id)
            element = atom.__class__.__name__
            x, y, z = (i for i in coords[atom_id])
            content.append(
                f' {round(x, 4)} {round(y, 4)} {round(z, 4)} '
                f'{element}\n'
            )

        content.append('$end\n')

        return content

    def to_string(
        self,
        molecule,
        atom_ids=None,
        periodic_info=None
    ):
        """
        Get a ``Turbomole`` file format string of `molecule`.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to ``Turbomole`` format.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        periodic_info : :class:`.PeriodicInfo`
            Information about the periodic cell.

        Returns
        -------
        :class:`string`
            The content of a `.coord` file.

        """

        content = self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_info=periodic_info,
        )

        return ''.join(content)

    def write(self, molecule, path, atom_ids=None, periodic_info=None):
        """
        Write `molecule` to ``Turbomole`` file format.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to ``Turbomole`` format.

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
