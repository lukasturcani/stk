"""
Bond Batch
==========

"""

from __future__ import annotations

from typing import Iterable

from ......bonds import Bond, BondInfo
from ......molecules import BuildingBlock


class _BondBatch:
    """
    A batch of bonds.

    """

    __slots__ = ['_bonds', '_bond_infos']

    _bonds: list[Bond]
    _bond_infos: list[BondInfo]

    def __init__(
        self,
        bonds: Iterable[Bond],
        id_map: dict[int, int],
        building_block: BuildingBlock,
        building_block_id: int,
    ) -> None:
        """
        Initialize a :class:`.BondBatch` instance.

        Parameters:

            bonds:
                The bonds, which should be added to the batch.

            id_map:
                Maps the ids of atoms held by `bonds`, to the new
                atoms, which the bonds in the :class:`.BondBatch`
                instance should hold.

            building_block:
                The building block from which the bonds originate.

            building_block_id:
                An id, unique to that building block and placement.

        """

        self._bonds = []
        self._bond_infos = []

        for bond in bonds:
            self._bonds.append(bond.with_ids(id_map))
            self._bond_infos.append(
                BondInfo(
                    bond=self._bonds[-1],
                    building_block=building_block,
                    building_block_id=building_block_id,
                )
            )

    def get_bonds(self) -> Iterable[Bond]:
        """
        Yield the bonds in the batch.

        Yields:

            A bond in the batch.

        """

        yield from self._bonds

    def get_bond_infos(self) -> Iterable[BondInfo]:
        """
        Yield info about the bonds in the batch.

        Yields:

            Info about a bond in the batch.

        """

        yield from self._bond_infos
