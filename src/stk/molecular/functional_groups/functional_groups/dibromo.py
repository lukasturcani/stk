"""
Dibormo
=======

"""

from .generic_functional_group import GenericFunctionalGroup


class Dibromo(GenericFunctionalGroup):
    """
    Represents a dibromo functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[bromine1][atom1][atom2][bromine2]``.

    """

    def __init__(
        self,
        bromine1,
        atom1,
        bromine2,
        atom2,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Dibromo` instance.

        Parameters
        ----------
        bromine1 : :class:`.Br`
            The ``[bromine1]`` atom.

        atom1 : :class:`.Atom`
            The ``[atom1]`` atom.

        bromine2 : :class:`.Br`
            The ``[bromine2]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._bromine1 = bromine1
        self._atom1 = atom1
        self._bromine2 = bromine2
        self._atom2 = atom2
        super().__init__(
            atoms=(bromine1, atom1, bromine2, atom2),
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

    def get_bromine1(self):
        """
        Get the ``[bromine1]`` atom.

        Returns
        -------
        :class:`.Br`
            The ``[bromine1]`` atom.

        """

        return self._bromine1

    def get_atom2(self):
        """
        Get the ``[atom2]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_bromine2(self):
        """
        Get the ``[bromine2]`` atom.

        Returns
        -------
        :class:`.Br`
            The ``[bromine2]`` atom.

        """

        return self._bromine2

    def clone(self):
        clone = super().clone()
        clone._atom1 = self._atom1
        clone._bromine1 = self._bromine1
        clone._atom2 = self._atom2
        clone._bromine2 = self._bromine2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._bromine1 = atom_map.get(
            self._bromine1.get_id(),
            self._bromine1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._bromine2 = atom_map.get(
            self._bromine2.get_id(),
            self._bromine2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._bromine1}, {self._atom1}, {self._bromine2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
