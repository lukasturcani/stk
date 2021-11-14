"""
Fluoro
======

"""

from .generic_functional_group import GenericFunctionalGroup


class Fluoro(GenericFunctionalGroup):
    """
    Represents a fluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine][atom]``.

    """

    def __init__(
        self,
        fluorine,
        atom,
        bonders,
        deleters,
        placers=None,
    ):
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

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._fluorine = fluorine
        self._atom = atom
        super().__init__(
            atoms=(fluorine, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

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
