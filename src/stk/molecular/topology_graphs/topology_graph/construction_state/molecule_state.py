import numpy as np
from collections import defaultdict

from ....atoms import AtomInfo
from ....bonds import BondInfo


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
        atoms = self._atoms
        atom_infos = self._atom_infos
        bonds = self._bonds
        bond_infos = self._bond_infos
        positions = []
        deleted_ids = set()

        def get_id(item):
            return item.get_id()

        def with_result(result):
            atom_map = {}

            def with_new_atom(atom):
                atoms.append(atom.with_id(len(atoms)))
                atom_infos.append(AtomInfo(atoms[-1], None, None))
                atom_map[atom.get_id()] = atoms[-1]

            def with_new_bond(bond):
                bonds.append(bond.with_atoms(atom_map))
                bond_infos.append(BondInfo(bonds[-1], None, None))

            for atom, position in result.get_new_atoms():
                with_new_atom(atom)
                positions.append(position)

            for bond in result.get_new_bonds():
                with_new_bond(bond)

            deleted_ids.update(map(get_id, result.get_deleted_atoms()))

        for reaction in reactions:
            with_result(reaction)

        positions = (
            np.vstack([self._position_matrix, positions])
            if positions
            else self._position_matrix
        )

        def valid_atom(atom):
            return atom.get_id() not in deleted_ids

        valid_atoms = []
        valid_atom_infos = []
        valid_positions = []
        atom_map = {}

        def with_valid_atom(atom):
            atom_id = atom.get_id()
            valid_atoms.append(atom.with_id(len(valid_atoms)))
            valid_positions.append(positions[atom_id])
            atom_map[atom_id] = valid_atoms[-1]

            info = atom_infos[atom_id]
            valid_atom_infos.append(
                AtomInfo(
                    atom=valid_atoms[-1],
                    building_block=info.get_building_block(),
                    building_block_id=info.get_building_block_id(),
                )
            )

        for atom in filter(valid_atom, atoms):
            with_valid_atom(atom)

        self._atoms = valid_atoms
        self._atom_infos = valid_atom_infos
        self._position_matrix = np.vstack(valid_positions)

        def valid_bond(bond_data):
            index, bond = bond_data
            return (
                bond.get_atom1().get_id() not in deleted_ids
                and bond.get_atom2.get_id() not in deleted_ids
            )

        valid_bonds = []
        valid_bond_infos = []

        def with_valid_bond(index, bond):
            valid_bonds.append(bond.with_atoms(atom_map))
            info = bond_infos[index]
            valid_bond_infos.append(
                BondInfo(
                    bond=valid_bonds[-1],
                    building_block=info.get_building_block(),
                    building_block_id=info.get_building_block_id(),
                )
            )

        for index, bond in filter(valid_bond, enumerate(bonds)):
            with_valid_bond(index, bond)

        self._bonds = valid_bonds
        self._bond_infos = valid_bond_infos
