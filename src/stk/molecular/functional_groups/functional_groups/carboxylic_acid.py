from .generic_functional_group import GenericFunctionalGroup


class CarboxylicAcid(GenericFunctionalGroup):
    """
    Represents a carboxylic acid functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen1])[oxygen2][hydrogen]``.

    """

    def __init__(
        self,
        carbon,
        oxygen1,
        oxygen2,
        hydrogen,
        atom,
        bonders,
        deleters,
    ):
        """
        Initialize a :class:`.CarboxylicAcid` instance.

        Parameters
        ----------
        carbon : :class:`.C`
            The ``[carbon]`` atom.

        oxygen1 : :class:`.O`
            The ``[oxygen1]`` atom.

        oxygen2 : :class:`.O`
            The ``[oxygen2]`` atom.

        hydrogen : :class:`.H`
            The ``[hydrogen]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon = carbon
        self._oxygen1 = oxygen1
        self._oxygen2 = oxygen2
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen1, oxygen2, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

    def get_carbon(self):
        """
        Get the ``[carbon]`` atom.

        Returns
        -------
        :class:`.C`
            The ``[carbon]`` atom.

        """

        return self._carbon

    def get_oxygen1(self):
        """
        Get the ``[oxygen1]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen]`` atom.

        """

        return self._oxygen1

    def get_oxygen2(self):
        """
        Get the ``[oxygen2]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen2]`` atom.

        """

        return self._oxygen2

    def get_hydrogen(self):
        """
        Get the ``[hydrogen]`` atom.

        Returns
        -------
        :class:``
            The ``[hydrogen]`` atom.

        """

        return self._hydrogen

    def get_atom(self):
        """
        Get the ``[atom]`` atom.

        Returns
        -------
        :class:``
            The ``[atom]`` atom.

        """

        return self._atom

    def clone(self):
        clone = super().clone()
        clone._carbon = self._carbon
        clone._oxygen1 = self._oxygen1
        clone._oxygen2 = self._oxygen2
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {
            atom.get_id(): index
            for index, atom in enumerate(self._atoms)
        }
        d.update({
            'carbon': indices[self._carbon.get_id()],
            'oxygen1': indices[self._oxygen1.get_id()],
            'oxygen2': indices[self._oxygen2.get_id()],
            'hydrogen': indices[self._hydrogen.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._carbon = obj._atoms[functional_group['carbon']]
        obj._oxygen1 = obj._atoms[functional_group['oxygen1']]
        obj._oxygen2 = obj._atoms[functional_group['oxygen2']]
        obj._hydrogen = obj._atoms[functional_group['hydrogen']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._carbon = atom_map.get(
            self._carbon.get_id(),
            self._carbon,
        )
        clone._oxygen1 = atom_map.get(
            self._oxygen1.get_id(),
            self._oxygen1,
        )
        clone._oxygen2 = atom_map.get(
            self._oxygen2.get_id(),
            self._oxygen2,
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
            f'{self._carbon}, {self._oxygen1}, {self._oxygen2}, '
            f'{self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
