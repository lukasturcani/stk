"""
Molecule State
==============

"""

import numpy as np

from .deletions_summary import _DeletionsSummary
from .placements_summary import _PlacementsSummary
from .reactions_summary import _ReactionsSummary


class _MoleculeState:
    """
    Represents the state of a molecule under construction.

    """

    __slots__ = [
        '_position_matrix',
        '_atoms',
        '_atom_infos',
        '_bonds',
        '_bond_infos',
        '_edge_functional_groups',
        '_num_placements',
    ]

    def __init__(self):
        """
        Initialize a :class:`._MoleculeState` instance.

        """

        self._position_matrix = np.empty((0, 3), dtype=np.float64)
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = {}
        self._num_placements = 0

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
        clone._num_placements = self._num_placements
        clone._edge_functional_groups = {
            edge_id: list(functional_groups)
            for edge_id, functional_groups
            in self._edge_functional_groups.items()
        }
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
            num_atoms=len(self._atoms),
            num_previous_placements=self._num_placements,
        )
        self._num_placements += len(vertices)
        self._atoms.extend(summary.get_atoms())
        self._atom_infos.extend(summary.get_atom_infos())
        self._bonds.extend(summary.get_bonds())
        self._bond_infos.extend(summary.get_bond_infos())
        self._position_matrix = np.concatenate([
            self._position_matrix,
            summary.get_position_matrix(),
        ])
        for edge_id, functional_groups in (
            summary.get_edge_functional_groups()
        ):
            self._edge_functional_groups[edge_id] = (
                self._edge_functional_groups.get(edge_id, [])
            )
            self._edge_functional_groups[edge_id].extend(
                functional_groups
            )
        return self

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
        Yield the functional groups associated with `edge_group`.

        Parameters
        ----------
        edge_group : :class:`.EdgeGroup`
            The edge group, whose functional groups are desired.

        Yields
        ------
        :class:`.FunctionalGroup`
            A functional group which belongs to `edge_group`.

        """

        for edge_id in edge_group.get_edge_ids():
            yield from self._edge_functional_groups[edge_id]

    def with_position_matrix(self, position_matrix):
        """
        Return a clone holding the `position_matrix`.

        Parameters
        ----------
        position_matrix : :class:`numpy.ndarray`
            The position matrix of the clone. The shape of the matrix
            is ``(n, 3)``.

        Returns
        -------
        :class:`._MoleculeState`
            The clone holding the new position matrix. Has the same
            type as the original instance.

        """

        return self.clone()._with_position_matrix(position_matrix)

    def _with_position_matrix(self, position_matrix):
        """
        Modify the instance.

        """

        self._position_matrix = np.array(position_matrix)
        return self

    def with_reaction_results(self, reactions, results):
        """
        Return a clone holding the reaction results.

        Parameters
        ----------
        reactions : :class:`tuple` of :class:`.Reaction`
            The reactions.

        results : :class:`.ReactionResult`
            For each reaction in `reactions`, its result.

        Returns
        -------
        :class:`._MoleculeState`
            The clone holding the reaction results. Has the same type
            as the original instance.

        """

        return self.clone()._with_reaction_results(reactions, results)

    def _with_reaction_results(self, reactions, results):
        """
        Modify the instance.

        """

        reactions_summary = _ReactionsSummary(
            num_atoms=len(self._atoms),
            reaction_results=results,
        )
        self._with_reactions_summary(reactions_summary)
        self._with_deletions_summary(_DeletionsSummary(
            atoms=self._atoms,
            atom_infos=self._atom_infos,
            bonds=self._bonds,
            bond_infos=self._bond_infos,
            position_matrix=self._position_matrix,
            deleted_atom_ids=reactions_summary.get_deleted_atom_ids(),
            deleted_bond_ids=reactions_summary.get_deleted_bond_ids(),
        ))
        return self

    def _with_reactions_summary(self, summary):
        """
        Add the results held in `summary`.

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

        positions = tuple(summary.get_positions())
        if positions:
            self._position_matrix = np.vstack([
                self._position_matrix,
                positions,
            ])

    def _with_deletions_summary(self, summary):
        """
        Add the results held in `summary`.

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
        self._position_matrix = np.array(list(summary.get_positions()))
        return self
