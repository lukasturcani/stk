class ReactionResult:
    def __init__(self, new_atoms, new_bonds, deleted_atoms):
        self._new_atoms = tuple(new_atoms)
        self._new_bonds = tuple(new_bonds)
        self._deleted_atoms = tuple(deleted_atoms)

    def with_new_ids_from(self, id):
        return ReactionResult(
            new_atoms=map(self._update_atom(id), self._new_atoms),
            new_bonds=map(self._update_bond(id), self._new_bonds),
            deleted_atoms=self._deleted_atoms,
        )

    @staticmethod
    def _update_atom(ids_from):

        def update_atom(atom):
            if atom.get_id() < 0:
                return atom.with_id(ids_from - atom.get_id() - 1)
            else:
                return atom

        return update_atom

    def _update_bond(self, ids_from):

        update_atom = self._update_atom(ids_from)

        def update_bond(bond):
            atom1 = update_atom(bond.get_atom1())
            atom2 = update_atom(bond.get_atom2())
            atom_map = {}
            if atom1 is not bond.get_atom1():
                atom_map[atom1.get_id()] = atom1
            if atom2 is not bond.get_atom2():
                atom_map[atom2.get_id()] = atom2

            return bond.with_atoms(atom_map) if atom_map else bond

        return update_bond

    def get_new_atoms(self):
        yield from self._new_atoms

    def get_new_bonds(self):
        yield from self._bonds

    def get_deleted_atoms(self):
        yield from self._deleted_atoms

    def get_num_new_atoms(self):
        return len(self._new_atoms)
