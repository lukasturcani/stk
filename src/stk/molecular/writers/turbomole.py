"""
TurbomoleWriter
===============

.. toctree::
    :maxdepth: 2

    TurbomoleWriter <stk.molecular.writers.turbomole>

"""


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

        raise NotImplementedError()

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

        raise NotImplementedError()
