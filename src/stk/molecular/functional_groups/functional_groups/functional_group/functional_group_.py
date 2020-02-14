from .functional_group import FunctionalGroup


class FunctionalGroup_(FunctionalGroup):
    """
    A partial implementation of :class:`.FunctionalGroup`.

    """

    def __init__(self, atoms):
        """
        Initialize a :class:`.FunctionalGroup`.

        Parameters
        ----------
        atoms : :class:`tuple` of :class:`.Atom`
            The atoms in the functional group.

        """

        self._atoms = atoms

    def clone(self):
        clone = self.__class__.__new__(self.__class__)
        FunctionalGroup_.__init__(clone, self._atoms)
        return clone

    def _with_atoms(self, atom_map):
        """
        Modify the functional group.

        """

        self._atoms = tuple(
            atom_map.get(a.get_id(), a) for a in self._atoms
        )
        return self

    def with_atoms(self, atom_map):
        return self.clone()._with_atoms(atom_map)

    def get_atoms(self):
        yield from self._atoms

    def get_atom_ids(self):
        yield from (a.get_id() for a in self._atoms)

    def __str__(self):
        return repr(self)
