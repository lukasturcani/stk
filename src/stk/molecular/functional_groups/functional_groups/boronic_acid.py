from .generic_functional_group import GenericFunctionalGroup


class BoronicAcid(GenericFunctionalGroup):
    """
    Represents a boronic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][boron]([oxygen1][hydrogen1])[oxygen2][hydrogen2]``.

    """

    def __init__(
        self,
        boron,
        oxygen1,
        hydrogen1,
        oxygen2,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        """
        Initialize a :class:`.BoronicAcid` instance.

        Parameters
        ----------
        boron : :class:`.B`
            The ``[boron]`` atom.

        oxygen1 : :class:`.O`
            The ``[oxygen1]`` atom.

        hydrogen1 : :class:`.H`
            The ``[hydrogen]`` atom.

        oxygen2 : :class:`.O`
            The ``[oyxgen2]`` atom.

        hydrogen2 : :class:`.H`
            The ``[hydrogen2]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            the bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._boron = boron
        self._oxygen1 = oxygen1
        self._hydrogen1 = hydrogen1
        self._oxygen2 = oxygen2
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (boron, oxygen1, hydrogen1, oxygen2, hydrogen2, atom)
        super().__init__(atoms, bonders, deleters)

    def get_boron(self):
        """
        Get the ``[boron]`` atom.

        Returns
        -------
        :class:`.B`
            The ``[boron]`` atom.

        """

        return self._boron

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
        clone._boron = self._boron
        clone._oxygen1 = self._oxygen1
        clone._hydrogen1 = self._hydrogen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen2 = self._hydrogen2
        clone._atom = self._atom
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {
            atom.get_id(): index
            for index, atom in enumerate(self._atoms)
        }
        d.update({
            'boron': indices[self._boron.get_id()],
            'oxygen1': indices[self._oxygen1.get_id()],
            'hydrogen1': indices[self._hydrogen1.get_id()],
            'oxygen2': indices[self._oxygen2.get_id()],
            'hydrogen2': indices[self._hydrogen2.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._boron = obj._atoms[functional_group['boron']]
        obj._oxygen1 = obj._atoms[functional_group['oxygen1']]
        obj._hydrogen1 = obj._atoms[functional_group['hydrogen1']]
        obj._oxygen2 = obj._atoms[functional_group['oxygen2']]
        obj._hydrogen2 = obj._atoms[functional_group['hydrogen2']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._boron = atom_map.get(
            self._boron.get_id(),
            self._boron,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
        )
        clone._hydrogen2 = atom_map.get(
            self._hydrogen2.get_id(),
            self._hydrogen2,
        )
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._boron}, {self._oxygen1}, {self._hydrogen1}, '
            f'{self._oxygen2}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
