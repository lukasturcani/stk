"""
Construction Result
===================

"""

from collections import abc
import numpy as np

from ....atom import Atom
from ....atom_info import AtomInfo
from ....bond import Bond
from ....bond_info import BondInfo
from ....building_block import BuildingBlock
from ..construction_state import ConstructionState


__all__ = (
    'ConstructionResult',
)


class ConstructionResult:
    """
    The result of :meth:`.TopologyGraph.construct`.

    """

    __slots__ = [
        '_atoms',
        '_bonds',
        '_atom_infos',
        '_bond_infos',
        '_position_matrix',
        '_num_building_blocks',
    ]

    def __init__(
        self,
        construction_state: ConstructionState,
    ) -> None:
        """
        Initialize a :class:`.ConstructionResult`.

        Parameters:

            construction_state:
                The state from which the result is initialized.

        """

        self._position_matrix = (
            construction_state.get_position_matrix()
        )
        self._position_matrix.setflags(write=False)
        self._atoms = tuple(construction_state.get_atoms())
        self._bonds = tuple(construction_state.get_bonds())
        self._atom_infos = tuple(construction_state.get_atom_infos())
        self._bond_infos = tuple(construction_state.get_bond_infos())
        self._num_building_blocks = {
            building_block: construction_state.get_num_building_block(
                building_block=building_block,
            )
            for building_block
            in construction_state.get_building_blocks()
        }

    def get_position_matrix(self) -> np.ndarray:
        """
        Get the position matrix of the constructed molecule.

        Returns:

            The position matrix of the constructed molecule.

        """

        return self._position_matrix

    def get_atoms(self) -> tuple[Atom, ...]:
        """
        Get the atoms of the constructed molecule.

        Returns:

            The atoms of the constructed molecule.

        """

        return self._atoms

    def get_bonds(self) -> tuple[Bond, ...]:
        """
        Get the bonds of the constructed molecule.

        Returns:

            The bonds of the constructed molecule.

        """

        return self._bonds

    def get_atom_infos(self) -> tuple[AtomInfo, ...]:
        """
        Get the atom infos of the constructed molecule.

        Returns:

            The atom infos of the constructed molecule.

        """

        return self._atom_infos

    def get_bond_infos(self) -> tuple[BondInfo, ...]:
        """
        Get the bond infos of the constructed molecule.

        Returns:

            The bond infos of the constructed molecule.

        """

        return self._bond_infos

    def get_num_building_block(
        self,
        building_block: BuildingBlock,
    ) -> int:
        """
        Get the number of times `building_block` is present.

        Parameters:

            building_block:
                The building block whose frequency in the constructed
                molecule is desired.

        Returns:

            The number of times `building_block` was used in the
            construction of the constructed molecule.

        """

        return self._num_building_blocks[building_block]

    def get_building_blocks(self) -> abc.Iterator[BuildingBlock]:
        """
        Yield the building blocks.

        Building blocks are yielded in an order based on their
        position in the constructed molecule. For two topologically
        equivalent constructed molecules, but with different building
        blocks, equivalently positioned building blocks will be
        yielded at the same time.

        Yields:

            A building block of the topology graph.

        """

        yield from self._num_building_blocks
