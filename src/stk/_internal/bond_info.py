from stk._internal.bond import Bond
from stk._internal.molecule import Molecule


class BondInfo:
    """
    Holds additional info about :class:`.ConstructedMolecule` bonds.
    """

    def __init__(
        self,
        bond: Bond,
        building_block: Molecule | None,
        building_block_id: int | None,
    ) -> None:
        """
        Parameters:
            bond:
                The bond about which information is held.

            building_block:
                The building block from which the bond originates.
                Can be ``None``, if the bond was not part of a building
                block, but was added by the construction process instead.

            building_block_id:
                A unique id for each :class:`.Molecule` placed during
                the construction of the :class:`.ConstructedMolecule`. As a
                single :class:`.Molecule` can be placed multiple times
                during construction, the `building_block_id` allows
                the user to distinguish between each placement. Can be
                ``None``, if the bond was not part of a building block, but
                was added by the construction process instead.
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

    def get_building_block(self) -> Molecule | None:
        """
        Get the building block from which the bond originates.

        Returns:
            The building block or ``None``, if the bond was not
            originally found in a building block, but was added by
            the construction process instead.
        """
        return self._building_block

    def get_building_block_id(self) -> int | None:
        """
        Get the id of the bond's building block.

        A unique id for each :class:`.Molecule` placed during
        the construction of the :class:`.ConstructedMolecule`. As a
        single :class:`.Molecule` can be placed multiple times
        during construction, the building block id  allows
        the user to distinguish between each placement.

        Returns:
            The id, if the bond was not originally found in a
            building block, but was added by the construction
            process instead.
        """
        return self._building_block_id
