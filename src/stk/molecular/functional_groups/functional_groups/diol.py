"""
Diol
====

"""

from .generic_functional_group import GenericFunctionalGroup


class Diol(GenericFunctionalGroup):
    """
    Represents a diol functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[hydrogen1][oxygen1][atom1][atom2][oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        atom1,
        oxygen1,
        hydrogen1,
        atom2,
        oxygen2,
        hydrogen2,
        bonders,
        deleters,
        placers=None,
    ):
        """
        Initialize a :class:`.Diol` instance.

        Parameters
        ----------
        atom1 : :class:`.Atom`
            The ``[atom1]`` atom.

        oxygen1 : :class:`.O`
            The ``[oxygen1]`` atom.

        hydrogen1 : :class:`.H`
            The ``[hydrogen1]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom2]`` atom.

        oxygen2 : :class:`.O`
            The ``[oxygen2]`` atom.

        hydrogen2 : :class:`.H`
            The ``[hydrogen2]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        placers : :class:`tuple` of :class:`.Atom`, optional
            The placer atoms. If ``None`` the `bonders` will be used.

        """

        self._atom1 = atom1
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._atom2 = atom2
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        atoms = (atom1, oxygen1, hydrogen1, atom2, oxygen2, hydrogen2)
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

    def get_oxygen1(self):
        """
        Get the ``[oxygen1]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen1]`` atom.

        """

        return self._oxygen1

    def get_hydrogen1(self):
        """
        Get the ``[hydrogen1]`` atom.

        Returns
        -------
        :class:`.H`
            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

    def get_atom2(self):
        """
        Get the ``[atom2]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom2]`` atom.

        """

        return self._atom2

    def get_oxygen2(self):
        """
        Get the ``[oxygen2]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen2]`` atom.

        """

        return self._oxygen2

    def get_hydrogen2(self):
        """
        Get the ``[hydrogen2]`` atom.

        Returns
        -------
        :class:`.H`
            The ``[hydrogen2]`` atom.

        """

        return self._hydrogen2

    def clone(self):
        clone = super().clone()
        clone._atom1 = self._atom1
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._atom2 = self._atom2
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._atom1}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._atom2}, {self._oxygen2}, {self._hydrogen2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
