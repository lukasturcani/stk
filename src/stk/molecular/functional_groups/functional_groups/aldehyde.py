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

        """

        self._carbon = carbon
        self._oxygen = oxygen
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (carbon, oxygen, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

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

    def to_dict(self):
        d = super().to_dict()
        indices = {
            atom.get_id(): index
            for index, atom in enumerate(self._atoms)
        }
        d.update({
            'carbon': indices[self._carbon.get_id()],
            'oxygen': indices[self._oxgyen.get_id()],
            'hydrogen': indices[self._hydrogen.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    @classmethod
    def _init_from_dict(cls, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._carbon = obj._atoms[functional_group['carbon']]
        obj._oxygen = obj._atoms[functional_group['oxygen']]
        obj._hydrogen = obj._atoms[functional_group['hydrogen']]
        obj._atom = obj._atoms[functional_group['atom']]
        return obj

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'{self._carbon}, {self._oxygen}, {self._hydrogen}, '
            f'{self._atom}, bonders={self._bonders}, '
            f'deleters={self._deleters})'
        )
