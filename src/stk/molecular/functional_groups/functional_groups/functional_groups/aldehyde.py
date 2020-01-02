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
            carbon.get_id(): carbon.clone(),
            oxygen.get_id(): oxygen.clone(),
            hydrogen.get_id(): hydrogen.clone(),
            atom.get_id(): atom.clone(),
        }
        self._carbon = atom_map[carbon.get_id()]
        self._oxygen = atom_map[oxygen.get_id()]
        self._hydrogen = atom_map[hydrogen.get_id()]
        self._atom = atom_map[atom.get_id()]
        super()._init(
            atoms=tuple(atom_map.values()),
            bonders=tuple(atom_map[a.get_id()] for a in bonders),
            deleters=tuple(atom_map[a.get_id()] for a in deleters),
        )

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
            if atom.get_id() not in atom_map:
                atom_map[atom.get_id()] = atom.clone()

        clone = super().clone(atom_map)
        clone._carbon = atom_map[self._carbon.get_id()]
        clone._oxygen = atom_map[self._oxygen.get_id()]
        clone._hydrogen = atom_map[self._hydrogen.get_id()]
        clone._atom = atom_map[self._atom.get_id()]
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._hydrogen}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters})'
        )
