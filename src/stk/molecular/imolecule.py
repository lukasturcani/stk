"""
IMolecule
=========

"""

from __future__ import annotations

from collections import abc
import typing

__all__ = (
    'IAtom',
    'IBond',
    'IMolecule',
)


class IMolecule(typing.Protocol):
    """
    A generic interface for molecules.

    Any object which obeys this interface can be used for initializing
    a :class:`.BuildingBlock`.

    """

    def get_atoms(self) -> abc.Iterable[IAtom]:
        """
        Get the atoms of the molecule.

        Yields:

            An atom of the molecule.

        """

        pass

    def get_bonds(self) -> abc.Iterable[IBond]:
        """
        Get the bonds of the molecule.

        Yields:

            A bond of the molecule.

        """

        pass


class IAtom(typing.Protocol):
    """
    A generic interface for atoms.

    Objects which obey this interface can be used for initializing
    a :class:`.BuildingBlock`.

    """

    def get_atomic_number(self) -> int:
        """
        Get the atomic number.

        Returns:

            The atomic number of the atom.

        """

        pass

    def get_charge(self) -> int:
        """
        Get the charge.

        Returns:

            The charge of the atom.

        """

        pass


class IBond(typing.Protocol):
    """
    A generic interface for bonds.

    Objects which obey this interface can be used for initializing
    a :class:`.BuildingBlock`.

    """

    def get_atom1_id(self) -> int:
        """
        Get the id of the first atom in the bond.

        """

        pass

    def get_atom2_id(self) -> int:
        """
        Get the id of the second atom in the bond.

        """

        pass

    def get_order(self) -> int:
        """
        Get the bond order.

        Returns:

            The bond order of the bond.

        """

        pass
