from .generic_functional_group import GenericFunctionalGroup


class Amide(GenericFunctionalGroup):
    """
    Represents an amide functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][carbon](=[oxygen])[nitrogen]([hydrogen1])[hydrogen2]``.

    """

    def __init__(
        self,
        carbon,
        oxygen,
        nitrogen,
        hydrogen1,
        hydrogen2,
        atom,
        bonders,
        deleters,
    ):
        """
        Initialize a :class:`.Amide` instance.

        Parameters
        ----------
        carbon : :class:`.C`
            The ``[carbon]`` atom.

        oxygen : :class:`.O`
            The ``[oxygen]` atom.

        nitrogen : :class:`.N`
            The ``[nitrogen]`` atom.

        hydrogen1 : :class:`.H`
            The ``[hydrogen1]`` atom.

        hydrogen2 : :class:`.H`
            The ``[hydrogen2]`` atom.

        atom : :class:`.Atom`
            The ``[atom]`` atom.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms.

        """

        self._carbon = carbon
        self._oxygen = oxygen
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._atom = atom
        atoms = (carbon, oxygen, nitrogen, hydrogen1, hydrogen2, atom)
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

    def get_oxygen(self):
        """
        Get the ``[oxygen]`` atom.

        Returns
        -------
        :class:`.O`
            The ``[oxygen]`` atom.

        """

        return self._oxygen

    def get_nitrogen(self):
        """
        Get the ``[nitrogen]`` atom.

        Returns
        -------
        :class:`.N`
            The ``[nitrogen]`` atom.

        """

        return self._nitrogen

    def get_hydrogen1(self):
        """
        Get the ``[hydrogen1]`` atom.

        Returns
        -------
        :class:`.H`
            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

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
        clone._carbon = self._carbon
        clone._oxygen = self._oxygen
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
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
            'carbon': indices[self._carbon.get_id()],
            'oxygen': indices[self._oxygen.get_id()],
            'nitrogen': indices[self._nitrogen.get_id()],
            'hydrogen1': indices[self._hydrogen1.get_id()],
            'hydrogen2': indices[self._hydrogen2.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._carbon = obj._atoms[functional_group['carbon']]
        obj._oxygen = obj._atoms[functional_group['oxygen']]
        obj._nitrogen = obj._atoms[functional_group['nitrogen']]
        obj._hydrogen1 = obj._atoms[functional_group['hydrogen1']]
        obj._hydrogen2 = obj._atoms[functional_group['hydrogen2']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

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
        clone._nitrogen = atom_map.get(
            self._nitrogen.get_id(),
            self._nitrogen,
        )
        clone._hydrogen1 = atom_map.get(
            self._hydrogen1.get_id(),
            self._hydrogen1,
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
            f'{self._carbon}, {self._oxygen}, {self._nitrogen}, '
            f'{self._hydrogen1}, {self._hydrogen2}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
