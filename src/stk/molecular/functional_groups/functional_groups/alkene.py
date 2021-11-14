"""
Alkene
======

"""

from .generic_functional_group import GenericFunctionalGroup


class Alkene(GenericFunctionalGroup):
    """
    Represents an alkene functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[carbon1]([atom1])([atom2])=[carbon2]([atom3])[atom4]``.

    """

    def __init__(
        self,
        carbon1,
        atom1,
        atom2,
        carbon2,
        atom3,
        atom4,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Alkene` instance.

        Parameters
        ----------
        carbon1 : :class:`.C`
            The ``[carbon1]`` atom.

        atom1 : :class:`.Atom`
            The ``[atom1]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom2]`` atom.

        carbon2 : :class:`.C`
            The ``[carbon2]`` atom.

        atom3 : :class:`.Atom`
            The ``[atom3]`` atom.

        atom4 : :class:`.Atom`
            The ``[atom4]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._carbon1 = carbon1
        self._atom1 = atom1
        self._atom2 = atom2
        self._carbon2 = carbon2
        self._atom3 = atom3
        self._atom4 = atom4
        atoms = (carbon1, atom1, atom2, carbon2, atom3, atom4)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_carbon1(self):
        """
        Get the ``[carbon1]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon1]`` atom.

        """

        return self._carbon1

    def get_atom1(self):
        """
        Get the ``[atom1]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom1]`` atom.

        """

        return self._atom1

    def get_atom2(self):
        """
        Get the ``[atom2]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_carbon2(self):
        """
        Get the ``[carbon2]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon2]`` atom.

        """

        return self._carbon2

    def get_atom3(self):
        """
        Get the ``[atom3]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom3]`` atom.

        """

        return self._atom3

    def get_atom4(self):
        """
        Get the ``[atom4]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom4]`` atom.

        """

        return self._atom4

    def clone(self):
        clone = super().clone()
        clone._carbon1 = self._carbon1
        clone._atom1 = self._atom1
        clone._atom2 = self._atom2
        clone._carbon2 = self._carbon2
        clone._atom3 = self._atom3
        clone._atom4 = self._atom4
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
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        clone._atom3 = atom_map.get(
            self._atom3.get_id(),
            self._atom3,
        )
        clone._atom4 = atom_map.get(
            self._atom4.get_id(),
            self._atom4,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon1}, {self._atom1}, {self._atom2}, '
            f'{self._carbon2}, {self._atom3}, {self._atom4}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
