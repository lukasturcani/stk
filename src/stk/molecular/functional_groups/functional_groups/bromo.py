from .generic_functional_group import GenericFunctionalGroup


class Bromo(GenericFunctionalGroup):
    """
    Represents a bromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine][atom]``.

    """

    def __init__(self, bromine, atom, bonders, deleters):
        """
        Initialize a :class:`.Bromo` instance.

        Parameters
        ----------
        bromine : :class:`.Br`
            The ``[bromine]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._bromine = bromine
        self._atom = atom
        super().__init__((bromine, atom), bonders, deleters)

    def get_bromine(self):
        """
        Get the ``[bromine]`` atom.

        Returns
        -------
        :class:`.Br`
            The ``[bromine]`` atom.

        """

        return self._bromine

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
        clone._bromine = self._bromine
        clone._atom = self._atom
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {
            atom.get_id(): index
            for index, atom in enumerate(self._atoms)
        }
        d.update({
            'bromine': indices[self._bromine.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._bromine = obj._atoms[functional_group['bromine']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._bromine = atom_map.get(
            self._bromine.get_id(),
            self._bromine,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
