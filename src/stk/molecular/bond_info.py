"""
Bond Info
=========

"""

import typing

from .bond import Bond
from .molecule import Molecule

__all__ = (
    'BondInfo',
)


class BondInfo:
    """
    Holds additional info about :class:`.ConstructedMolecule` bonds.

    """

    def __init__(
        self,
        bond: Bond,
        building_block: typing.Optional[Molecule],
        building_block_id: typing.Optional[int],
    ) -> None:
        """
        Initialize an :class:`.BondInfo` instance.

        Parameters:

            bond:
                The bond about which information is held.

            building_block:
                The building block from which the bond originates.
                Can be ``None``, if the bond was not part of a building
                block, but was added by the construction process
                instead.

            building_block_id:
                A unique id for each :class:`.Molecule` placed during
                the construction of the :class:`.ConstructedMolecule`.
                As a single :class:`.Molecule` can be placed multiple
                times during construction, the `building_block_id`
                allows the user to distinguish between each placement.
                Can be ``None``, if the bond was not part of a building
                block, but was added by the construction process
                instead.

        """
        self._bond = bond
        self._building_block = building_block
        self._building_block_id = building_block_id

    def get_bond(self) -> Bond:
        """
        Get the bond about which information is held.

        Returns:

            The bond.

        """

        return self._bond

    def get_building_block(self) -> typing.Optional[Molecule]:
        """
        Get the building block from which the bond originates.

        Returns:

            The building block or ``None`` if the bond was not
            originally found in a building block, but was added by the
            construction process instead.

        """

        return self._building_block

    def get_building_block_id(self) -> typing.Optional[int]:
        """
        Get the id of the bond's building block.

        A unique id for each :class:`.Molecule` placed during
        the construction of the :class:`.ConstructedMolecule`. As a
        single :class:`.Molecule` can be placed multiple times
        during construction, the building block id  allows
        the user to distinguish between each placement.

        Returns:
            The id or ``None`` if the bond was not originally found in
            a building block, but was added by the construction process
            instead.

        """

        return self._building_block_id
