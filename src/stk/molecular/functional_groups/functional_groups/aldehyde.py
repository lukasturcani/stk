"""
Aldehyde
========

"""

from .generic_functional_group import GenericFunctionalGroup


class Aldehyde(GenericFunctionalGroup):
    """
    Represents an aldehyde functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        hydrogen,
        atom,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Aldehyde` instance.

        Parameters
        ----------
        carbon : :class:`.C`
            The carbon atom.

        oxygen : :class:`.O`
            The oxygen atom.

        hydrogen : :class:`.H`
            The hydrogen atom.

        atom : :class:`.Atom`
            The atom to which the functional group is attached.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._carbon = carbon
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, hydrogen, atom)
        super().__init__(
            atoms=atoms,
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_carbon(self):
        """
        Get the carbon atom.

        Returns
        -------
        :class:`.C`
            The carbon atom.

        """

        return self._carbon

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

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
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

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._hydrogen}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters})'
        )
