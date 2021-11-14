"""
Alkyne
======

"""

from .generic_functional_group import GenericFunctionalGroup


class Alkyne(GenericFunctionalGroup):
    """
    Represents an alkyne functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon1]([atom1])#[carbon2][atom2]``.

    """

    def __init__(
        self,
        carbon1,
        atom1,
        carbon2,
        atom2,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Alkyne` instance.

        Parameters
        ----------
        carbon1 : :class:`.C`
            The ``[carbon1]`` atom.

        atom1 : :class:`.Atom`
            The ``[atom1]`` atom.

        carbon2 : :class:`.C`
            The ``[carbon2]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom2]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
                The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._carbon1 = carbon1
        self._atom1 = atom1
        self._carbon2 = carbon2
        self._atom2 = atom2
        atoms = (carbon1, atom1, carbon2, atom2)
        super().__init__(
            atoms=atoms,
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

    def get_carbon1(self):
        """
        Get the ``[carbon1]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon1]`` atom.

        """

        return self._carbon1

    def get_carbon2(self):
        """
        Get the ``[carbon2]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon2]`` atom.

        """

        return self._carbon2

    def get_atom2(self):
        """
        Get the ``[atom2]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom2]`` atom.

        """

        return self._atom2

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._atom1 = self._atom1
        clone._carbon2 = self._carbon2
        clone._atom2 = self._atom2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._atom1}, {self._carbon2}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
