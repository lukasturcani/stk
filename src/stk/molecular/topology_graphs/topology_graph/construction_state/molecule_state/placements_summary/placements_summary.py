import numpy as np
from collections import defaultdict

from .atom_batch import _AtomBatch
from .bond_batch import _BondBatch


class _PlacementsSummary:
    """
    A summary of placement results.

    """

    __slots__ = [
        '_atoms',
        '_atom_infos',
        '_bonds',
        '_bond_infos',
        '_edge_functional_groups',
        '_position_matrices',
        '_num_atoms',
    ]

    def __init__(
        self,
        building_blocks,
        placement_results,
        num_atoms,
        num_previous_placements,
    ):
        """
        Initialize a :class:`._PlacementsSummary` instance.

        Parameters
        ----------
        building_blocks : :class:`iterable` of :class:`.BuildingBlock`
            The building blocks which were placed.

        placement_results : :class:`iterable`
            Holds a :class:`_PlacementResults` instance for each
            building block in `building_blocks`.

        num_atoms : :class:`int`
            The number of atoms in molecule being constructed,
            before this summary is taken into account.

        num_previous_placements : :class:`int`
            The total number of building block placements done
            previously.

        """

        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = defaultdict(list)
        self._position_matrices = []
        # This will get updated as placement results are added to the
        # summary.
        self._num_atoms = num_atoms

        for id_, (building_block, result) in enumerate(
            zip(building_blocks, placement_results),
            num_previous_placements,
        ):
            self._with_placement_result(building_block, id_, result)
            self._num_atoms += building_block.get_num_atoms()

    def _with_placement_result(
        self,
        building_block,
        building_block_id,
        result,
    ):
        """
        Add the placement result to the summary.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block which was placed.

        building_block_id : :class:`int`
            A unique for the building block placement.

        result : :class:`._PlacementResult`
            The result of the placement.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._position_matrices.append(result.position_matrix)

        atom_batch = _AtomBatch(
            atoms=building_block.get_atoms(),
            num_atoms=self._num_atoms,
            building_block=building_block,
            building_block_id=building_block_id,
        )
        self._with_atom_batch(atom_batch)
        atom_map = atom_batch.get_atom_map()

        bond_batch = _BondBatch(
            bonds=building_block.get_bonds(),
            atom_map=atom_map,
            building_block=building_block,
            building_block_id=building_block_id,
        )
        self._with_bond_batch(bond_batch)

        self._with_functional_group_edges(
            building_block=building_block,
            functional_group_edges=result.functional_group_edges,
            atom_map=atom_map,
        )

    def _with_atom_batch(self, batch):
        """
        Add a batch of atoms to the summary.

        Parameters
        ----------
        batch : :class:`._AtomBatch`
            A batch of atoms.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._atoms.extend(batch.get_atoms())
        self._atom_infos.extend(batch.get_atom_infos())

    def _with_bond_batch(self, batch):
        """
        Add a batch of bonds to the summary.

        Parameters
        ----------
        batch : :class:`.BondBatch`
            A batch of bonds.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._bonds.extend(batch.get_bonds())
        self._bond_infos.extend(batch.get_bond_infos())

    def _with_functional_group_edges(
        self,
        building_block,
        functional_group_edges,
        atom_map,
    ):
        """
        Add the mapping from functional groups to edges.

        Parameters
        ----------
        building_block : :class:`.BuildingBlock`
            The building block which owns the functional groups.

        functional_group_edges : :class:`dict`
            Maps the id of each functional group of `building_block`
            the the id of an edge.

        atom_map : :class:`dict`
            Maps the ids of atoms in `building_block` to the new atoms
            of the molecule being constructed.

        Returns
        -------
        None : :class:`NoneType`

        """

        functional_groups = building_block.get_functional_groups(
            fg_ids=functional_group_edges,
        )
        edge_ids = functional_group_edges.values()
        functional_group_edges = zip(functional_groups, edge_ids)
        for functional_group, edge_id in functional_group_edges:
            self._edge_functional_groups[edge_id].append(
                functional_group.with_atoms(atom_map)
            )

    def get_atoms(self):
        """
        Yield the atoms in the summary.

        Yields
        ------
        :class:`.Atom`
            An atom.

        """

        yield from self._atoms

    def get_atom_infos(self):
        """
        Yield infos about atoms in the summary.

        Yields
        ------
        :class:`.AtomInfo`
            Info about an atom.

        """

        yield from self._atom_infos

    def get_bonds(self):
        """
        Yield the bonds in the summary.

        Yields
        ------
        :class:`.Bond`
            A bond.

        """

        yield from self._bonds

    def get_bond_infos(self):
        """
        Yield infos about the bonds in the summary.

        Yields
        ------
        :class:`.BondInfo`
            Info about a bond.

        """

        yield from self._bond_infos

    def get_position_matrix(self):
        """
        Get a position matrix for the atoms in the summary.

        Returns
        -------
        :class:`numpy.ndarray`
            The position matrix.

        """

        return np.vstack(self._position_matrices)

    def get_edge_functional_groups(self):
        """
        Yield the edge ids and functional groups associated with them.

        Yields
        ------
        :class:`tuple`
            Holds the edge id as the first element, and a
            :class:`tuple` of :class:`.FunctionalGroup` instances
            as the second element.

        """

        yield from self._edge_functional_groups.items()
