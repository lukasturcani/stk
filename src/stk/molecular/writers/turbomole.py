"""
TurbomoleWriter
===============

.. toctree::
    :maxdepth: 2

    TurbomoleWriter <stk.molecular.writers.turbomole>

"""

from .utilities import cell_matrix_to_lengths_angles


class TurbomoleWriter:
    """
    A writer class for ``Turbomole`` files.

    Examples
    --------
    *Writing to File with Unit Cell*

    This writer can write to file for visualisation with the unit cell
    included for periodic molecules. Note that this always assumes P1
    space group.

    .. code-block:: python

        import stk

        bb1 = stk.BuildingBlock('BrCCBr', [stk.BromoFactory()])
        bb2 = stk.BuildingBlock('BrCC(CBr)CBr', [stk.BromoFactory()])
        topology_graph = stk.cof.Honeycomb(
            building_blocks=(bb1, bb2),
            lattice_size=(3, 3, 1),
            periodic=True,
        )
        cof = stk.ConstructedMolecule(topology_graph)
        writer = stk.TurbomoleWriter()
        writer.write(
            molecule=cof,
            path='cof.coord',
            periodic_cell=topology_graph.get_periodic_cell()
        )

    """

    def _write_content(self, molecule, atom_ids, periodic_cell=None):

        if atom_ids is None:
            atom_ids = range(molecule.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        content = []
        if periodic_cell is not None:
            # Input unit cell information.
            lengths_and_angles = cell_matrix_to_lengths_angles(
                periodic_cell
            )
            a, b, c, alpha, beta, gamma = lengths_and_angles
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

    def write_string(
        self,
        molecule,
        atom_ids=None,
        periodic_cell=None
    ):
        """
        Write `molecule` to ``Turbomole`` file format as string.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            Molecule to write to ``Turbomole`` format.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        periodic_cell : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        Returns
        -------
        string :class:`string`
            The string containing all information for a `.coord` file.

        """

        content = self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_cell=periodic_cell,
        )

        return ''.join(content)

    def write(self, molecule, path, atom_ids=None, periodic_cell=None):
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

        periodic_cell : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        Returns
        -------
        None : :class:`NoneType`
            A file is written.

        """

        content = self._write_content(
            molecule=molecule,
            atom_ids=atom_ids,
            periodic_cell=periodic_cell,
        )

        with open(path, 'w') as f:
            f.write(''.join(content))
