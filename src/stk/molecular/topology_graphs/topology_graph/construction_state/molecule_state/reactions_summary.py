from .....atoms import AtomInfo
from .....bonds import BondInfo


class _ReactionsSummary:
    def __init__(self, next_id, reaction_results):
        self._next_id = next_id
        self._atoms = []
        self._atom_infos = []
        self._positions = []
        self._bonds = []
        self._bond_infos = []
        self._deleted_ids = set()

        for result in reaction_results:
            self._with_reaction_result(result)

    def _with_reaction_result(self, result):
        atoms = self._atoms
        atom_infos = self._atom_infos
        positions = self._positions
        atom_map = {}

        def with_atom(atom, id):
            atoms.append(atom.with_id(id))
            atom_infos.append(AtomInfo(atoms[-1], None, None))
            atom_map[atom.get_id()] = atoms[-1]
            positions.append(position)

        new_atoms = result.get_new_atoms()
        next_id = self._next_id
        for id_, (atom, position) in enumerate(new_atoms, next_id):
            with_atom(atom, id_)
        self._next_id += len(new_atoms)

        bonds = self._bonds
        bond_infos = self._bond_infos

        def with_bond(bond):
            bonds.append(bond.with_atoms(atom_map))
            bond_infos.append(BondInfo(bonds[-1], None, None))

        for bond in result.get_new_bonds():
            with_bond(bond)

    def get_atoms(self):
        yield from self._atoms

    def get_atom_infos(self):
        yield from self._atom_infos

    def get_bonds(self):
        yield from self._bonds

    def get_bond_infos(self):
        yield from self._bond_infos

    def get_deleted_ids(self):
        yield from self._deleted_ids
