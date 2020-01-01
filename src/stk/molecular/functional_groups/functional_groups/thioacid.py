from .. import FunctionalGroup_


class Thioacid(FunctionalGroup_):
    """
    Represents a thioacid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[sulfur][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        sulfur,
        hydrogen,
        atom,
        bonders,
        deleters
    ):
        self._carbon = carbon
        self._oxygen = oxygen
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

    def get_sulfur(self):
        return self._sulfur.clone()

    def get_hydrogen(self):
        return self._hydrogen.clone()

    def get_atom(self):
        return self._atom.clone()

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}
        else:
            atom_map = dict(atom_map)

        atoms = (
            self._carbon,
            self._oxygen,
            self._sulfur,
            self._hydrogen,
            self._atom,
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.id]
        clone._oxygen = atom_map[self._oxygen.id]
        clone._sulfur = atom_map[self._sulfur.id]
        clone._hydrogen = atom_map[self._hydrogen.id]
        clone._atom = atom_map[self._atom.id]
        return clone
