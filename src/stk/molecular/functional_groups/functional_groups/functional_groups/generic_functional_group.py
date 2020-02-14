from ..functional_group import FunctionalGroup_


class GenericFunctionalGroup(FunctionalGroup_):
    """
    A functional group which defines general atom classes.

    Bonders are atoms which should have bonds added by a
    :class:`.Reaction`. Deleters are atoms which should be removed
    by a :class:`.Reaction`.

    This interface allows the same reactions to be carried out across
    different functional groups.

    """

    def __init__(self, atoms, bonders, deleters):
        """
        Initialize a :class:`.GenericFunctionalGroup`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms in the functional group.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms in the functional group.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms in the functional group.

        """

        super().__init__(atoms)
        self._bonders = bonders
        self._deleters = deleters

    def clone(self):
        clone = super().clone()
        clone._bonders = self._bonders
        clone._deleters = self._deleters
        return clone

    def _with_atoms(self, atom_map):
        super()._with_atoms(atom_map)
        self._bonders = tuple(
            atom_map.get(a.get_id(), a) for a in self._bonders
        )
        self._deleters = tuple(
            atom_map.get(a.get_id(), a) for a in self._deleters
        )
        return self


    def get_bonders(self):
        """
        Yield bonder atoms in the functional group.

        These are atoms which have bonds added during
        :class:`.ConstructedMolecule` construction.

        Yields
        ------
        :class:`.Atom`
            A bonder atom.

        """

        yield from self._bonders

    def get_num_bonders(self):
        """
        Get the number of bonder atoms.

        Returns
        -------
        :class:`int`
            The number of bonder atoms.

        """

        return len(self._bonders)

    def get_bonder_ids(self):
        """
        Yield the ids of bonder atoms.

        Yields
        ------
        :class:`int`
            The id of a bonder :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._bonders)

    def get_deleters(self):
        """
        Yield the deleter atoms in the functional group.

        These are atoms which are removed during
        :class:`.ConstructedMolecule` construction.

        Yields
        ------
        :class:`.Atom`
            A deleter atom.

        """

        yield from self._deleters

    def get_deleter_ids(self):
        """
        Yield the ids of deleter atoms.

        Yields
        -------
        :class:`int`
            The id of a deleter :class:`.Atom`.

        """

        yield from (a.get_id() for a in self._deleters)

    def get_positional_atom_ids(self):
        yield from (a.get_id() for a in self._bonders)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'atoms={self._atoms}, '
            f'bonders={self._bonders}, '
            f'deleters={self._deleters}'
            ')'
        )
