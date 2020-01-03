from .functional_group import FunctionalGroup


class FunctionalGroup_(FunctionalGroup):
    """
    An implementation of the :class:`.FunctionalGroup` interface.

    """

    def __init__(self, atoms, bonders, deleters):
        """
        Initialize a :class:`.FunctionalGroup`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms in the functional group.

        bonders : :class:`tuple` of :class:`.Atom`
            The bonder atoms in the functional group.

        deleters : :class:`tuple` of :class:`.Atom`
            The deleter atoms in the functional group.

        """

        self._atoms = atoms
        self._bonders = bonders
        self._deleters = deleters

    def clone(self):

        clone = self.__class__.__new__(self.__class__)
        for attr, value in vars(self).items():
            if not attr.startswith('_'):
                setattr(clone, attr, value)

        FunctionalGroup_.__init__(
            self=clone,
            atoms=self._atoms,
            bonders=self._bonders,
            deleters=self._deleters,
        )
        return clone

    def _with_atoms(self, atom_map):
        """
        Modify the functional group.

        """

        self._atoms = tuple(
            atom_map.get(a.get_id(), a) for a in self._atoms
        )
        self._bonders = tuple(
            atom_map.get(a.get_id(), a) for a in self._bonders
        )
        self._deleters = tuple(
            atom_map.get(a.get_id(), a) for a in self._deleters
        )
        return self

    def with_atoms(self, atom_map):
        """
        Return a clone holding different atoms.

        Parameters
        ----------
        atom_map : :class:`dict`
            Maps the id of an atom in the functional group to the new
            atom the clone should hold. If the id of an atom in the
            functional group is not found in `atom_map`, the atom will
            not be replaced in the clone.

        Returns
        -------
        :class:`.FunctionalGroup`
            The clone.

        """

        return self.clone()._with_atoms(atom_map)

    def get_atoms(self):
        yield from self._atoms

    def get_atom_ids(self):
        yield from (a.get_id() for a in self._atoms)

    def get_bonders(self):
        yield from self._bonders

    def get_num_bonders(self):
        return len(self._bonders)

    def get_bonder_ids(self):
        yield from (a.get_id() for a in self._bonders)

    def get_deleters(self):
        yield from self._deleters

    def get_deleter_ids(self):
        yield from (a.get_id() for a in self._deleters)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'atoms={self._atoms}, '
            f'bonders={self._bonders}, '
            f'deleters={self._deleters}'
            ')'
        )

    def __str__(self):
        return repr(self)
