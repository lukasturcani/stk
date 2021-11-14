"""
Single Atom
===========

"""

from .generic_functional_group import GenericFunctionalGroup


class SingleAtom(GenericFunctionalGroup):
    """
    Represents an abstract single atom functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[atom]``.

    """

    def __init__(self, atom):
        """
        Initialize a :class:`.SingleAtom` instance.

        Parameters
        ----------
        atom : :class:`.Atom`
            Any :class:`.Atom` will work.

        """

        self._atom = atom
        atoms = (atom, )
        super().__init__(atoms, atoms, ())

    def get_atom(self):
        """
        Get the atom that defines the functional group.

        Returns
        -------
        :class:`.Atom`
            The atom.

        """

        return self._atom

    def with_atoms(self, atom_map):

        clone = super().with_atoms(atom_map)
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        return clone

    def clone(self):
        clone = super().clone()
        clone._atom = self._atom
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}({self._atom})'
        )
