"""
Metal Bound Atom
================

"""

from .generic_functional_group import GenericFunctionalGroup


class MetalBoundAtom(GenericFunctionalGroup):
    """
    Represents an atom that is bound to a metal atom functional group.

    The structure of the functional group is given by the pseudo-SMILES
    ``[metal][atom]``.

    """

    def __init__(self, atom, metal):
        """
        Initialize a :class:`.MetalBoundAtom` instance.

        Parameters
        ----------
        atom : :class:`.Atom`
            The ``[atom]`` atom.

        metal : :class:`.Atom`
            The ``[metal]`` atom.

        """

        self._atom = atom
        self._metal = metal
        bonders = (atom, )
        super().__init__((atom, metal), bonders, ())

    def get_atom(self):
        """
        Get the ``[atom]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[atom]`` atom.

        """

        return self._atom

    def get_metal(self):
        """
        Get the ``[metal]`` atom.

        Returns
        -------
        :class:`.Atom`
            The ``[metal]`` atom.

        """

        return self._metal

    def clone(self):
        clone = super().clone()
        clone._atom = self._atom
        clone._metal = self._metal
        return clone

    def with_atoms(self, atom_map):
        clone = super().with_atoms(atom_map)
        clone._atom = atom_map.get(
            self._atom.get_id(),
            self._atom,
        )
        clone._metal = atom_map.get(
            self._metal.get_id(),
            self._metal,
        )
        return clone

    def __repr__(self):
        return (
            f'{self.__class__.__name__}({self._atom}, {self._metal})'
        )
