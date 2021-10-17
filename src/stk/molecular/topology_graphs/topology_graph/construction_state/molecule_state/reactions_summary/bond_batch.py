"""
Bond Batch
==========


"""

from collections import abc

from ......atom import Atom
from ......bond import Bond
from ......bond_info import BondInfo


__all__ = (
    'BondBatch',
)


class BondBatch:
    """
    A batch of bonds.

    """

    _bonds: list[Bond]
    _bond_infos: list[BondInfo]

    __slots__ = ['_bonds', '_bond_infos']

    def __init__(
        self,
        bonds: abc.Iterable[Bond],
        atom_map: dict[int, Atom],
    ) -> None:
        """
        Initialize a :class:`.BondBatch` instance.

        Parameters:

            bonds:
                The bonds, which should be added to the batch.

            atom_map:
                Maps the ids of atoms held by `bonds`, to the new
                atoms, which the bonds in the :class:`.BondBatch`
                instance should hold.

        """

        _bonds = self._bonds = []
        bond_infos = self._bond_infos = []

        for bond in bonds:
            _bonds.append(bond.with_atoms(atom_map))
            bond_infos.append(BondInfo(_bonds[-1], None, None))

    def get_bonds(self) -> abc.Iterator[Bond]:
        """
        Yield the bonds in the batch.

        Yields:

            A bond in the batch.

        """

        yield from self._bonds

    def get_bond_infos(self) -> abc.Iterator[BondInfo]:
        """
        Yield info about the bonds in the batch.

        Yields:

            Info about a bond in the batch.

        """

        yield from self._bond_infos
