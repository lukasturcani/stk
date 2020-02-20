from .generic_functional_group import GenericFunctionalGroup


class Fluoro(GenericFunctionalGroup):
    """
    Represents a fluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine][atom]``.

    """

    def __init__(self, fluorine, atom, bonders, deleters):
        """
        Initialize a :class:`.Fluoro` instance.

        Parameters
        ----------
        fluorine : :class:`.F`
            The ``[fluorine]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._fluorine = fluorine
        self._atom = atom
        super().__init__((fluorine, atom), bonders, deleters)

    def get_fluorine(self):
        """
        Get the ``[fluorine]`` atom.

        Returns
        -------
        :class:`.F`
            The ``[fluorine]`` atom.

        """

        return self._fluorine

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
        clone._fluorine = self._fluorine
        clone._atom = self._atom
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {
            atom.get_id(): index
            for index, atom in enumerate(self._atoms)
        }
        d.update({
            'fluorine': indices[self._fluorine.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._fluorine = obj._atoms[functional_group['fluorine']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._fluorine = atom_map.get(
            self._fluorine.get_id(),
            self._fluorine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
