from .. import FunctionalGroup_


class Bromo(FunctionalGroup_):
    """
    Represents a bromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][bromine]``.

    """

    def __init__(self, bromine, atom, bonders, deleters):
        self._bromine = bromine
        self._atom = atom
        super().__init__((bromine, atom), bonders, deleters)

    def get_bromine(self):
        return self._bromine.clone()

    def get_atom(self):
        return self._atom.clone()
