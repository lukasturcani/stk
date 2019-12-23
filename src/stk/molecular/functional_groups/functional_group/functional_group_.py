class FunctionalGroup_:
    """
    An implementation of :class:`.FunctionalGroup` interface.

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

    def clone(self, atom_map=None):
        if atom_map is None:
            atom_map = {}

        atom_map.update(
            (a.id, a.clone()) for a in self._atoms
            if a.id not in atom_map
        )

        clone = self.__class__.__new__(self.__class__)
        for attr, value in vars(self).items():
            if not attr.startswith('_'):
                setattr(clone, attr, value)

        _FunctionalGroup.__init__(
            self=clone,
            atoms=tuple(atom_map[a.id] for a in self._atoms),
            bonders=tuple(atom_map[a.id] for a in self._bonders),
            deleters=tuple(atom_map[a.id] for a in self._deleters),
        )
        return clone

    def get_atoms(self):
        yield from (a.clone() for a in self._atoms)

    def get_atom_ids(self):
        yield from (a.id for a in self._atoms)

    def get_bonders(self):
        yield from (a.clone() for a in self._bonders)

    def get_bonder_ids(self):
        yield from (a.id for a in self._bonders)

    def get_deleters(self):
        yield from (a.clone() for a in self._deleters)

    def get_deleter_ids(self):
        yield from (a.id for a in self._deleters)

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
