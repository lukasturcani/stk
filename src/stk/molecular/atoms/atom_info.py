"""
AtomInfo
========

"""

from dataclasses import dataclass


@dataclass(init=False, eq=False, frozen=True)
class AtomInfo:
    """
    Holds additional info about :class:`.ConstructedMolecule` atoms.

    Attributes
    ----------
    atom : :class:`.Atom`
        The atom about which information is held.

    building_block : :class:`.BuildingBlock`
        The building block from which the atom originates.

    building_block_index : :class:`int`
        A unique id for each :class:`.BuildingBlock` placed during the
        construction of the :class:`.ConstructedMolecule`. As a single
        :class:`.BuildingBlock` can be placed multiple times during
        construction, the :attr:`building_block_index` allows the user
        to distinguish between each placement.

    """

    atom: object
    building_block: object
    building_block_index: object
