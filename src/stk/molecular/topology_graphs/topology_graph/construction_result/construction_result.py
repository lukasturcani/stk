"""
Construction Result
===================

"""


class ConstructionResult:
    """
    The result of :meth:`.TopologyGraph.construct`.

    """

    __slots__ = [
        '_atoms',
        '_bonds',
        '_atom_infos',
        '_bond_infos',
        '_position_matrix',
        '_num_building_blocks',
    ]

    def __init__(self, construction_state):
        """
        Initialize a :class:`.ConstructionResult`.

        Parameters
        ----------
        construction_state : :class:`.ConstructionState`
            The state from which the result is initialized.

        """

        self._position_matrix = (
            construction_state.get_position_matrix()
        )
        self._position_matrix.setflags(write=False)
        self._atoms = tuple(construction_state.get_atoms())
        self._bonds = tuple(construction_state.get_bonds())
        self._atom_infos = tuple(construction_state.get_atom_infos())
        self._bond_infos = tuple(construction_state.get_bond_infos())
        self._num_building_blocks = {
            building_block: construction_state.get_num_building_block(
                building_block=building_block,
            )
            for building_block
            in construction_state.get_building_blocks()
        }

    def get_position_matrix(self):
        """
        Get the position matrix of the constructed molecule.

        Returns
        -------
        :class:`numpy.ndarray`
            The position matrix of the constructed molecule.

        """

        return self._position_matrix

    def get_atoms(self):
        """
        Get the atoms of the constructed molecule.

        Returns
        -------
        :class:`tuple` of :class:`.Atom`
            The atoms of the constructed molecule.

        """

        return self._atoms

    def get_bonds(self):
        """
        Get the bonds of the constructed molecule.

        Returns
        -------
        :class:`tuple` of :class:`.Bond`
            The bonds of the constructed molecule.

        """

        return self._bonds

    def get_atom_infos(self):
        """
        Get the atom infos of the constructed molecule.

        Returns
        -------
        :class:`tuple` of :class:`.AtomInfo`
            The atom infos of the constructed molecule.

        """

        return self._atom_infos

    def get_bond_infos(self):
        """
        Get the bond infos of the constructed molecule.

        Returns
        -------
        :class:`tuple` of :class:`.BondInfo`
            The bond infos of the constructed molecule.

        """

        return self._bond_infos

    def get_num_building_block(self, building_block):
        """
        Get the number of times `building_block` is present.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block whose frequency in the constructed
            molecule is desired.

        Returns
        -------
        :class:`int`
            The number of times `building_block` was used in the
            construction of the constructed molecule.

        """

        return self._num_building_blocks[building_block]

    def get_building_blocks(self):
        """
        Yield the building blocks.

        Building blocks are yielded in an order based on their
        position in the constructed molecule. For two topologically
        equivalent constructed molecules, but with different building
        blocks, equivalently positioned building blocks will be
        yielded at the same time.

        Yields
        ------
        :class:`.BuildingBlock`
            A building block of the topology graph.

        """

        yield from self._num_building_blocks
