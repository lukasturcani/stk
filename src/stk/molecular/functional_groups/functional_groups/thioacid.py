"""
Thioacid
========

"""

from .generic_functional_group import GenericFunctionalGroup


class Thioacid(GenericFunctionalGroup):
    """
    Represents a thioacid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[sulfur][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        sulfur,
        hydrogen,
        atom,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Thioacid` functional group.

        Parameters
        ----------
        carbon : :class:`.C`
            The ``[carbon]`` atom.

        oxygen : :class:`.O`
            The ``[oxygen]`` atom.

        sulfur : :class:`.S`
            The ``[sulfur]`` atom.

        hydrogen : :class:`.H`
            The ``[hydrogen]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._carbon = carbon
        self._oxygen = oxygen
        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        super().__init__(
            atoms=(carbon, oxygen, sulfur, hydrogen, atom),
            bonders=bonders,
            deleters=deleters,
            placers=bonders if placers is None else placers,
        )

    def get_carbon(self):
        """
        Get the ``[carbon]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon]`` atom.

        """

        return self._carbon

    def get_oxygen(self):
        """
        Get the ``[oxygen]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen]`` atom.

        """

        return self._oxygen

    def get_sulfur(self):
        """
        Get the ``[sulfur]`` atom.

        Returns
        -------
        :class:`.S`
            The ``[sulfur]`` atom.

        """

        return self._sulfur

    def get_hydrogen(self):
        """
        Get the ``[hydrogen]`` atom.

        Returns
        -------
        :class:`.H`
            Get the ``[hydrogen]`` atom.

        """

        return self._hydrogen

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
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._sulfur = self._sulfur
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

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
        clone._sulfur = atom_map.get(
            self._sulfur.get_id(),
            self._sulfur,
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
            f'{self._carbon}, {self._oxygen}, {self._sulfur}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
