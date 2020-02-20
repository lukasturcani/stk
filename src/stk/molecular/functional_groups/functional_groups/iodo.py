from .generic_functional_group import GenericFunctionalGroup


class Iodo(GenericFunctionalGroup):
    """
    Represents an iodo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[iodine][atom]``.

    """

    def __init__(self, iodine, atom, bonders, deleters):
        """
        Initialize a :class:`.Iodo` instance.

        Parameters
        ----------
        iodine : :class:`.I`
            The ``[iodine]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._iodine = iodine
        self._atom = atom
        super().__init__((iodine, atom), bonders, deleters)

    def get_iodine(self):
        """
        Get the ``[iodine]`` atom.

        Returns
        -------
        :class:`.I`
            The ``[iodine]`` atom.

        """

        return self._iodine

    def get_atom(self):
        """
        Get the ``[atom]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self):
        clone = super().clone()
        clone._iodine = self._iodine
        clone._atom = self._atom
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {
            atom.get_id(): index
            for index, atom in enumerate(self._atoms)
        }
        d.update({
            'iodine': indices[self._iodine.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._iodine = obj._atoms[functional_group['iodine']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._iodine = atom_map.get(
            self._iodine.get_id(),
            self._iodine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._iodine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
