from .. import FunctionalGroup_


class Iodo(FunctionalGroup_):
    """
    Represents an iodo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][iodine]``.

    """

    def __init__(self, iodine, atom, bonders, deleters):
        self._iodine = iodine
        self._atom = atom
        super().__init__((iodine, atom), bonders, deleters)

    def get_iodine(self):
        return self._iodine.clone()

    def get_atom(self):
        return self._atom.clone()
