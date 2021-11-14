"""
Alcohol
=======

"""

from .generic_functional_group import GenericFunctionalGroup


class Alcohol(GenericFunctionalGroup):
    """
    Represents an alcohol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][oxygen][hydrogen]``.

    """

    def __init__(
        self,
        oxygen,
        hydrogen,
        atom,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Alcohol` instance.

        Parameters
        ----------
        oxygen : :class:`.O`
            The oxygen atom.

        hydrogen : :class:`.H`
            The hydrogen atom.

        atom : :class:`.Atom`
            The atom to which the alcohol is attached.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (oxygen, hydrogen, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_oxygen(self):
        """
        Get the oxygen atom.

        Returns
        -------
        :class:`.O`
            The oxygen atom.

        """

        return self._oxygen

    def get_hydrogen(self):
        """
        Get the hydrogen atom.

        Returns
        -------
        :class:`.H`
            The hydrogen atom.

        """

        return self._hydrogen

    def get_atom(self):
        """
        Get the atom to which the functional group is attached.

        Returns
        -------
        :class:`.Atom`
            The atom to which the functional group is attached.

        """

        return self._atom

    def clone(self):
        clone = super().clone()
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._oxygen = atom_map.get(
            self._oxygen.get_id(),
            self._oxygen,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._oxygen}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
