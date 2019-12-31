from .. import FunctionalGroup_


class Fluoro(FunctionalGroup_):
    """
    Represents a fluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][fluorine]``.

    """

    def __init__(self, fluorine, atom, bonders, deleters):
        self._fluorine = fluorine
        self._atom = atom
        super().__init__((fluorine, atom), bonders, deleters)

    def get_fluorine(self):
        return self._fluorine.clone()

    def get_atom(self):
        return self._atom.clone()
