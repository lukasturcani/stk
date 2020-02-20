from .generic_functional_group import GenericFunctionalGroup


class Thiol(GenericFunctionalGroup):
    """
    Represents a thiol functional group.


    The structure of the functional group is given by the pseudo-SMILES
    ``[atom][sulfur][hydrogen]``.

    """

    def __init__(self, sulfur, hydrogen, atom, bonders, deleters):
        """
        Initialize a :class:`.Thiol` instance.

        Parameters
        ----------
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

        """

        self._sulfur = sulfur
        self._hydrogen = hydrogen
        self._atom = atom
        atoms = (sulfur, hydrogen, atom)
        super().__init__(atoms, bonders, deleters)

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
            The ``[hydrogen]`` atom.

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
        clone._sulfur = self._sulfur
        clone._hydrogen = self._hydrogen
        clone._atom = self._atom
        return clone

    def to_dict(self):
        d = super().to_dict()
        indices = {atom.get_id(): index for index, atom in self._atoms}
        d.update({
            'sulfur': indices[self._sulfur.get_id()],
            'hydrogen': indices[self._hydrogen.get_id()],
            'atom': indices[self._atom.get_id()],
        })
        return d

    def _init_from_dict(self, functional_group):
        obj = super()._init_from_dict(functional_group)
        obj._sulfur = self._atoms[functional_group['sulfur']]
        obj._hydrogen = self._atoms[functional_group['hydrogen']]
        obj._atom = self._atoms[functional_group['atom']]
        return obj

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
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
            f'{self._sulfur}, {self._hydrogen}, {self._atom}, '
            f'bonders={self._bonders}, deleters={self._deleters})'
        )
