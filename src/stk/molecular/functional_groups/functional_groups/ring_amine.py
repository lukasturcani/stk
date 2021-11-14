"""
Ring Amine
==========

"""

from __future__ import annotations

from .utilities import get_atom_map
from .functional_group import FunctionalGroup
from ...atoms import N, H, C


__all__ = (
    'RingAmine',
)


class RingAmine(FunctionalGroup):
    """
    Represents an amine bonded to a ring.

    The structure of the functional group is given by the pseudo-SMILES
    ``[hydrogen1][nitrogen]([hydrogen2])[carbon1][carbon2]
    ([hydrogen3])[carbon3]``.

    """

    def __init__(
        self,
        nitrogen: N,
        hydrogen1: H,
        hydrogen2: H,
        carbon1: C,
        carbon2: C,
        hydrogen3: H,
        carbon3: C,
    ) -> None:
        """
        Initializes a :class:`.RingAmine` instance.

        Parameters:

            nitrogen:
                The ``[nitrogen]`` atom.

            hydrogen1:
                The ``[hydrogen1]`` atom.

            hydrogen2:
                The ``[hydrogen2]`` atom.

            carbon1:
                The ``[carbon1]`` atom.

            carbon2:
                The ``[carbon2]`` atom.

            hydrogen3:
                The ``[hydrogen3]`` atom.

            carbon3:
                The ``[carbon3]`` atom.

        """

        FunctionalGroup.__init__(
            self=self,
            atoms=(
                nitrogen,
                hydrogen1,
                hydrogen2,
                carbon1,
                carbon2,
                hydrogen3,
                carbon3,
            ),
            placers=(nitrogen, ),
            core_atoms=(nitrogen, ),
        )
        self._nitrogen = nitrogen
        self._hydrogen1 = hydrogen1
        self._hydrogen2 = hydrogen2
        self._hydrogen3 = hydrogen3
        self._carbon1 = carbon1
        self._carbon2 = carbon2
        self._carbon3 = carbon3

    def get_nitrogen(self) -> N:
        """
        Get the ``[nitrogen]`` atom.

        Returns:

            The ``[nitrogen]`` atom.

        """

        return self._nitrogen

    def get_hydrogen1(self) -> H:
        """
        Get the ``[hydrogen1]`` atom.

        Returns:

            The ``[hydrogen1]`` atom.

        """

        return self._hydrogen1

    def get_hydrogen2(self) -> H:
        """
        Get the ``[hydrogen2]`` atom.

        Returns:

            The ``[hydrogen2]`` atom.

        """

        return self._hydrogen2

    def get_carbon1(self) -> C:
        """
        Get the ``[carbon1]`` atom.

        Returns:

            The ``[carbon1]`` atom.

        """

        return self._carbon1

    def get_carbon2(self) -> C:
        """
        Get the ``[carbon2]`` atom.

        Returns:

            The ``[carbon2]`` atom.

        """

        return self._carbon2

    def get_hydrogen3(self) -> H:
        """
        Get the ``[hydrogen3]`` atom.

        Returns:

            The ``[hydrogen3]`` atom.

        """

        return self._hydrogen3

    def get_carbon3(self) -> C:
        """
        Get the ``[carbon3]`` atom.

        Returns:

            The ``[carbon3]`` atom.

        """

        return self._carbon3

    def clone(self) -> RingAmine:
        clone = self._clone()
        clone._nitrogen = self._nitrogen
        clone._hydrogen1 = self._hydrogen1
        clone._hydrogen2 = self._hydrogen2
        clone._hydrogen3 = self._hydrogen3
        clone._carbon1 = self._carbon1
        clone._carbon2 = self._carbon2
        clone._carbon3 = self._carbon3
        return clone

    def with_ids(
        self,
        id_map: dict[int, int],
    ):
        atom_map = get_atom_map(
            id_map=id_map,
            atoms=(
                *self._atoms,
                *self._placers,
                *self._core_atoms,
                self._nitrogen,
                self._hydrogen1,
                self._hydrogen2,
                self._hydrogen3,
                self._carbon1,
                self._carbon2,
                self._carbon3,
            ),
        )
        clone = self.__class__.__new__(self.__class__)
        clone._atoms = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._atoms
        )
        clone._placers = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._placers
        )
        clone._core_atoms = tuple(
            atom_map.get(atom.get_id(), atom)
            for atom in self._core_atoms
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
        clone._hydrogen3 = atom_map.get(
            self._hydrogen3.get_id(),
            self._hydrogen3,
        )
        clone._carbon1 = atom_map.get(
            self._carbon1.get_id(),
            self._carbon1,
        )
        clone._carbon2 = atom_map.get(
            self._carbon2.get_id(),
            self._carbon2,
        )
        clone._carbon3 = atom_map.get(
            self._carbon3.get_id(),
            self._carbon3,
        )
        return clone

    def __repr__(self) -> str:
        return (
            f'{self.__class__.__name__}('
            f'{self._nitrogen}, {self._hydrogen1}, {self._hydrogen2}, '
            f'{self._carbon1}, {self._carbon2}, {self._hydrogen3}, '
            f'{self._carbon3})'
        )
