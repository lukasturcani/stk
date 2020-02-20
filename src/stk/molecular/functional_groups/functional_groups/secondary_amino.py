from .generic_functional_group import GenericFunctionalGroup


class SecondaryAmino(GenericFunctionalGroup):
    """
    Represents a secondary amino functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom1][nitrogen]([hydrogen])[atom2]``.

    """

    def __init__(
        self,
        nitrogen,
        hydrogen,
        atom1,
        atom2,
        bonders,
        deleters,
    ):
        """
        Initialize a :class:`.SecondaryAmine` instance.

        Parameters
        ----------
        nitrogen : :class:`.N`
            The ``[nitrogen]`` atom

        hydrogen : :class:`.H`
            The ``[hydrogen]`` atom.

        atom1 : :class:`.Atom`
            The ``[atom]`` atom.

        atom2 : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._nitrogen = nitrogen
        self._hydrogen = hydrogen
        self._atom1 = atom1
        self._atom2 = atom2
        atoms = (nitrogen, hydrogen, atom1, atom2)
        super().__init__(atoms, bonders, deleters)

    def get_nitrogen(self):
        """
        Get the ``[nitrogen]`` atom.

        Returns
        -------
        :class:`.N`
            The ``[nitrogen]`` atom.

        """

        return self._nitrogen

    def get_hydrogen(self):
        """
        Get the ``[hydrogen]`` atom.

        Returns
        -------
        :class:`.H`
            The ``[hydrogen]`` atom.

        """

        return self._hydrogen

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

    def clone(self):
        clone = super().clone()
        clone._nitrogen = self._nitrogen
        clone._hydrogen = self._hydrogen
        clone._atom1 = self._atom1
        clone._atom2 = self._atom2
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {atom.get_id(): index for index, atom in self._atoms}
        d.update({
            'nitrogen': indices[self._nitrogen.get_id()],
            'hydrogen': indices[self._hydrogen.get_id()],
            'atom1': indices[self._atom1.get_id()],
            'atom2': indices[self._atom2.get_id()],
        })
        return d

    def _init_from_dict(self, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._nitrogen = self._atoms[functional_group['nitrogen']]
        obj._hydrogen = self._atoms[functional_group['hydrogen']]
        obj._atom1 = self._atoms[functional_group['atom1']]
        obj._atom2 = self._atoms[functional_group['atom2']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen = atom_map.get(
            self._hydrogen.get_id(),
            self._hydrogen,
        )
        clone._atom1 = atom_map.get(
            self._atom1.get_id(),
            self._atom1,
        )
        clone._atom2 = atom_map.get(
            self._atom2.get_id(),
            self._atom2,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen}, {self._atom1}, '
            f'{self._atom2}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
