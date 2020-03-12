import numpy as np
from collections import defaultdict
from functools import partial

from stk.utilities import flatten
from ....atoms import AtomInfo
from ....bonds import BondInfo


class _PlacementsSummary:
    ...


class _ReactionsSummary:
    def __init__(self, num_atoms):
        self._num_atoms = num_atoms
        self._new_atoms = []

    def update(self, result):

        new_atoms = self._new_atoms
        num_new_atoms = len(self._atoms)
        new_atoms.extend(map(self._with_id, result.get_new_atoms()))

        atom_infos = map(self._get_atom_info, new_atoms[num_new_atoms:])
        self._atom_infos.extend(atom_infos)

        old_ids = (atom.get_id() for atom in result.get_new_atoms())
        atom_map = dict(zip(old_ids, new_atoms))

        new_bonds.extend(map(partial(with_atoms, atom_map), new_bonds))

    def _with_id(self, new_atom):
        id_ = -new_atom.get_atom().get_id() - 1
        return new_atom.with_id(self._num_atoms + id_)

    @staticmethod
    def _get_atom_info(new_atom):
        return AtomInfo(new_atom.get_atom(), None, None)


class _MoleculeState:
    def __init__(self):
        self._position_matrix = np.empty((0, 3), dtype=np.float64)
        self._atoms = []
        self._atom_infos = []
        self._bonds = []
        self._bond_infos = []
        self._edge_functional_groups = defaultdict(list)

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        clone._position_matrix = np.array(self._position_matrix)
        clone._atoms = list(self._atoms)
        clone._atom_infos = list(self._atom_infos)
        clone._bonds = list(self._bonds)
        clone._bond_infos = list(self._bond_infos)
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
        Modify the state.

        """

        # Doing a vstack after the loop should be faster than doing one
        # within the loop.
        position_matrices = [self._position_matrix]
        results_ = zip(building_blocks, results)
        for index, (building_block, result) in enumerate(results_):
            position_matrices.append(result.position_matrix)
            atom_map = self._add_atoms(index, building_block)
            self._add_bonds(index, building_block, atom_map)
            self._add_edge_functional_groups(
                building_block=building_block,
                result=result,
                atom_map=atom_map,
            )
        self._position_matrix = np.vstack(position_matrices)
        return self

    def _add_atoms(self, index, building_block):
        atom_map = {}
        for atom in building_block.get_atoms():
            new_atom = atom.with_id(len(self._atoms))
            atom_map[atom.get_id()] = new_atom
            self._atoms.append(new_atom)
            self._atom_infos.append(
                AtomInfo(
                    atom=new_atom,
                    building_block=building_block,
                    building_block_id=index,
                )
            )
        return atom_map

    def _add_bonds(self, index, building_block, atom_map):
        for bond in building_block.get_bonds():
            new_bond = bond.with_atoms(atom_map)
            self._bonds.append(new_bond)
            self._bond_infos.append(
                BondInfo(
                    bond=new_bond,
                    building_block=building_block,
                    building_block_id=index,
                )
            )

    def _add_edge_functional_groups(
        self,
        building_block,
        result,
        atom_map,
    ):
        # Add edge to functional group mappings.
        functional_groups = building_block.get_functional_groups(
            fg_ids=result.functional_group_edges,
        )
        edge_ids = result.functional_group_edges.values()
        functional_group_edges = zip(functional_groups, edge_ids)
        for functional_group, edge_id in functional_group_edges:
            self._edge_functional_groups[edge_id].append(
                functional_group.with_atoms(atom_map)
            )

    def with_placement_results(
        self,
        vertices,
        edges,
        building_blocks,
        results,
    ):
        return self.clone()._with_placement_results(
            vertices=vertices,
            edges=edges,
            building_blocks=building_blocks,
            results=results,
        )

    def get_position_matrix(self):
        """

        """

        return np.array(self._position_matrix)

    def get_atoms(self):
        yield from self._atoms

    def get_bonds(self):
        yield from self._bonds

    def get_atom_infos(self):
        yield from self._atom_infos

    def get_bond_infos(self):
        yield from self._bond_infos

    def get_edge_group_functional_groups(self, edge_group):
        for edge_id in edge_group.get_edge_ids():
            yield from self._edge_functional_groups[edge_id]

    def with_reaction_results(self, reactions, results):
        return self.clone()._with_reaction_results(reactions, results)

    def _with_reaction_results(self, reactions, results):
        summary = _ReactionsSummary(self._atoms, self._bonds)
        for reaction in reactions:
            summary.update(reaction)

        atoms = summary.get_atoms()
        atom_infos = summary.get_atom_infos()
        bonds = summary.get_bonds()
        bond_infos = summary.get_bond_infos()
        deleted_ids = summary.get_deleted_ids()
