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
    A writer class for `Turbomole` files.

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
            mol=cof,
            file='cof.coord',
            periodic_cell=topology_graph.get_periodic_cell()
        )

    """

    def _write_string(self, mol, atom_ids, periodic_cell=None):

        if atom_ids is None:
            atom_ids = range(mol.get_num_atoms())
        elif isinstance(atom_ids, int):
            atom_ids = (atom_ids, )

        lines = []
        if periodic_cell is not None:
            # Input unit cell information.
            lengths_and_angles = cell_matrix_to_lengths_angles(
                periodic_cell
            )
            a, b, c, alpha, beta, gamma = lengths_and_angles
            print(lengths_and_angles)
            lines.append(
                '$periodic 3\n'
                '$cell angs\n'
                f' {a:>8.3f} {b:>8.3f} {c:>8.3f} '
                f'{alpha:>6.2f} {beta:>6.2f} {gamma:>6.2f}\n'
            )

        coords = mol.get_position_matrix()
        lines.append('$cood angs\n')
        for atom_id in atom_ids:
            atom, = mol.get_atoms(atom_ids=atom_id)
            element = atom.__class__.__name__
            x, y, z = (i for i in coords[atom_id])
            lines.append(
                f' {round(x, 4)} {round(y, 4)} {round(z, 4)} '
                f'{element}\n'
            )

        lines.append('$end\n')

        return ''.join(lines)

    def write(self, mol, file=None, atom_ids=None, periodic_cell=None):
        """
        Write `mol` to Turbomole file format.

        Parameters
        ----------
        mol : :class:`.Molecule` `Turbomole` format.

        file : :class:`str`, optional
            The full path to the file being written.

        atom_ids : :class:`iterable` of :class:`int`
            The atom ids of atoms to write. Can be a single
            :class:`int`, if a single atom is to be used, or ``None``,
            if all atoms are to be used.

        periodic_cell : :class:`tuple` of :class:`np.array`
            Tuple of cell lattice vectors (shape: (3,)) in Angstrom.

        Returns
        -------
        string :class:`string`
            The string containing all information for a `.pdb` file.

        None : :class:`NoneType`
            A file is written if :attr:`file` is not `None`.

        """

        string = self._write_string(mol, atom_ids, periodic_cell)

        if file is None:
            return string
        else:
            with open(file, 'w') as f:
                f.write(string)
            return None
