"""
Difluoro
========

"""

from .generic_functional_group import GenericFunctionalGroup


class Difluoro(GenericFunctionalGroup):
    """
    Represents a difluoro functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[fluorine1][atom1][atom2][fluorine2]``.

    """

    def __init__(
        self,
        fluorine1,
        atom1,
        fluorine2,
        atom2,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Difluoro` instance.

        Parameters
        ----------
        fluorine1 : :class:`.F`
            The ``[fluorine1]`` atom.

        atom1 : :class:`.Atom`
            The ``[atom1]`` atom.

        fluorine2 : :class:`.F`
            The ``[fluorine2]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom2]`` atom.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._fluorine1 = fluorine1
        self._atom1 = atom1
        self._fluorine2 = fluorine2
        self._atom2 = atom2
        super().__init__(
            atoms=(fluorine1, atom1, fluorine2, atom2),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_atom1(self):
        """
        Get the ``[atom1]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_fluorine1(self):
        """
        Get the ``[fluorine1]`` atom.

        Returns
        -------
        :class:`.F`
            The ``[fluorine1]`` atom.

        """

        return self._fluorine1

    def get_atom2(self):
        """
        Get the ``[atom2]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_fluorine2(self):
        """
        Get the ``[fluorine2]`` atom.

        Returns
        -------
        :class:`.F`
            The ``[fluorine2]`` atom.

        """

        return self._fluorine2

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._fluorine1 = atom_map.get(
            self._fluorine1.get_id(),
            self._fluorine1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._fluorine2 = atom_map.get(
            self._fluorine2.get_id(),
            self._fluorine2,
        )
        return clone

    def clone(self):
        clone = super().clone()
        clone._atom1 = self._atom1
        clone._fluorine1 = self._fluorine1
        clone._atom2 = self._atom2
        clone._fluorine2 = self._fluorine2
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._fluorine1}, {self._atom1}, {self._fluorine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
