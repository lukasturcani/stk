"""
BondInfo
========

"""


class BondInfo:
    """
    Holds additional info about :class:`.ConstructedMolecule` bonds.

    """

    def __init__(self, bond, building_block, building_block_id):
        """
        Initialize an :class:`.BondInfo` instance.

        Parameters
        ----------
        bond : :class:`.Bond`
            The bond about which information is held.

        building_block : :class:`.BuildingBlock` or :class:`NoneType`
            The building block from which the bond originates.
            Can be ``None``, if the bond was not part of a building
            block, but was added by the construction process instead.

        building_block_id : :class:`int` or :class:`NoneType`
            A unique id for each :class:`.BuildingBlock` placed during
            the construction of the :class:`.ConstructedMolecule`. As a
            single :class:`.BuildingBlock` can be placed multiple times
            during construction, the `building_block_id` allows
            the user to distinguish between each placement. Can be
            ``None``, if the bond was not part of a building block, but
            was added by the construction process instead.

        """
        self._bond = bond
        self._building_block = building_block
        self._building_block_id = building_block_id

    def get_bond(self):
        return self._bond

    def get_building_block(self):
        return self._building_block

    def get_building_block_id(self):
        return self._building_block_id
