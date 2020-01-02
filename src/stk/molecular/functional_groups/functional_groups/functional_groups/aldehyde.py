from .. import FunctionalGroup_


class Aldehyde(FunctionalGroup_):
    """
    Represents an aldehyde functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        atom_map = {
            carbon.id: carbon.clone(),
            oxygen.id: oxygen.clone(),
            hydrogen.id: hydrogen.clone(),
            atom.id: atom.clone(),
        }
        self._carbon = atom_map[carbon.id]
        self._oxygen = atom_map[oxygen.id]
        self._hydrogen = atom_map[hydrogen.id]
        self._atom = atom_map[atom.id]
        atoms = (carbon, oxygen, hydrogen, atom)
        bonders = tuple(atom_map[a.id] for a in bonders)
        deleters = tuple(atom_map[a.id] for a in deleters)
        super()._init(atoms, bonders, deleters)

    def get_carbon(self):
        return self._carbon.clone()

    def get_oxygen(self):
        return self._oxygen.clone()

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
            self._carbon, self._oxygen, self._hydrogen, self._atom
        )
        for atom in atoms:
            if atom.id not in atom_map:
                atom_map[atom.id] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.id]
        clone._oxygen = atom_map[self._oxygen.id]
        clone._hydrogen = atom_map[self._hydrogen.id]
        clone._atom = atom_map[self._atom.id]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._hydrogen}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters})'
        )
