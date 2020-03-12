import numpy as np
from collections import defaultdict

from .reactions_summary import _ReactionsSummary
from .deletions_summary import _DeletionsSummary
from .placements_summary import _PlacementsSummary


class _MoleculeState:
    """
    Represents the state of a molecule under construction.

    """

    def __init__(self):
        """
        Initialize a :class:`._MoleculeState` instance.

        """

        self._position_matrix = np.empty((0, 3), dtype=np.float64)
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = defaultdict(list)
        self._num_building_blocks = 0

    def clone(self):
        """
        Return a clone.

        Returns
        -------
        :class:`._MoleculeState`
            The clone. Has the same type as the original instance.

        """

        clone = self.__class__.__new__(self.__class__)
        clone._position_matrix = np.array(self._position_matrix)
        clone._atoms = list(self._atoms)
        clone._atom_infos = list(self._atom_infos)
        clone._bonds = list(self._bonds)
        clone._bond_infos = list(self._bond_infos)
        clone._num_building_blocks = self._num_building_blocks
        clone._edge_functional_groups = defaultdict.copy(
            self._edge_functional_groups
        )
        return clone

    def _with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        """
        Modify the instance.

        """

        summary = _PlacementsSummary(
            building_blocks=building_blocks,
            placement_results=results,
            next_id=len(self._num_building_blocks),
        )
        self._num_building_blocks += len(building_blocks)
        self._atoms.extend(summary.get_atoms())
        self._atom_infos.extend(summary.get_atom_infos())
        self._bonds.extend(summary.get_bonds())
        self._bond_infos.extend(summary.get_bond_infos())
        self._edge_functional_groups.update(
            summary.get_edge_functional_groups()
        )
        self._position_matrix = np.concatenate([
            self._position_matrix,
            summary.get_position_matrix(),
        ])

    def with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        """
        Return a clone holding the placement results.

        Parameters
        ----------
        vertices : :class:`tuple` of :class:`.Vertex`
            The vertices used for placement.

        edges : :class:`tuple`
            For each vertex in `vertices`, a :class:`tuple` of
            :class:`.Edge` instances connected to it.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            For each vertex in `vertices`, the building block placed
            on it.

        results : :class:`tuple` of :class:`._PlacementResult`
            For every vertex in `vertices`, the result of the
            placement.

        Returns
        -------
        :class:`._MoleculeState`
            The clone holding the placement results. Has the same
            type as the original instance.

        """

        return self.clone()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

    def get_position_matrix(self):
        """
        Get the position matrix of the molecule.

        Returns
        -------
        :class:`numpy.ndarray`
            The position matrix.

        """

        return np.array(self._position_matrix)

    def get_atoms(self):
        """
        Yield the atoms of the molecule.

        Yields
        ------
        :class:`.Atom`
            An atom of the molecule.

        """

        yield from self._atoms

    def get_bonds(self):
        """
        Yield the bonds of the molecule.

        Yields
        ------
        :class:`.Bond`
            A bond of the molecule.

        """

        yield from self._bonds

    def get_atom_infos(self):
        """
        Yield the atom infos of the molecule.

        Yields
        ------
        :class:`.AtomInfo`
            An atom info of the molecule.

        """

        yield from self._atom_infos

    def get_bond_infos(self):
        """
        Yield the bond infos of the molecule.

        Yields
        ------
        :class:`.BondInfo`
            The bond info of the molecule.

        """

        yield from self._bond_infos

    def get_edge_group_functional_groups(self, edge_group):
        """

        """

        for edge_id in edge_group.get_edge_ids():
            yield from self._edge_functional_groups[edge_id]

    def with_reaction_results(self, reactions, results):
        """

        """
        return self.clone()._with_reaction_results(reactions, results)

    def _with_reaction_results(self, reactions, results):
        """
        Modify the instance.

        """

        reactions_summary = _ReactionsSummary(
            next_id=len(self._atoms),
            reaction_results=results,
        )
        self._with_reactions_summary(reactions_summary)
        self._with_deletions_summary(_DeletionsSummary(
            atoms=self._atoms,
            atom_infos=self._atom_infos,
            bonds=self._bonds,
            bond_infos=self._bond_infos,
            position_matrix=self._position_matrix,
            deleted_ids=reactions_summary.get_deleted_ids(),
        ))

    def _with_reactions_summary(self, summary):
        """

        Parameters
        ----------
        summary : :class:`._ReactionsSummary`
            A summary of the reaction results.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._atoms.extend(summary.get_atoms())
        self._atom_infos.extend(summary.get_atom_infos())
        self._bonds.extend(summary.get_bonds())
        self._bond_infos.extend(summary.get_bond_infos())

        positions = summary.get_positions()
        if positions:
            self._position_matrix = np.vstack([
                self._position_matrix,
                positions,
            ])

    def _with_deletions_summary(self, summary):
        """

        Parameters
        ----------
        summary : :class:`._DeletionsSummary`
            A summary of the atom deletion results.

        Returns
        -------
        None : :class:`NoneType`

        """

        self._atoms = list(summary.get_atoms())
        self._atom_infos = list(summary.get_atom_infos())
        self._bonds = list(summary.get_bonds())
        self._bond_infos = list(summary.get_bond_infos())
        self._position_matrix = np.array(summary.get_positions())
