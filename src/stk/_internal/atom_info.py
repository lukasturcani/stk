"""
Atom Info
=========

"""


import typing

from stk._internal.atom import Atom
from stk._internal.molecule import Molecule


class AtomInfo:
    """
    Holds additional info about :class:`.ConstructedMolecule` atoms.

    """

    def __init__(
        self,
        atom: Atom,
        building_block_atom: typing.Optional[Atom],
        building_block: typing.Optional[Molecule],
        building_block_id: typing.Optional[int],
    ) -> None:
        """
        Initialize an :class:`.AtomInfo` instance.

        Parameters:

            atom:
                The atom about which information is held.

            building_block_atom:
                The building block atom from which this atom
                originates. Can be ``None``, if the atom was not part
                of the building block, but was added by the
                construction process instead.

            building_block:
                The building block from which the atom originates.
                Can be ``None``, if the atom was not part of a building
                block, but was added by the construction process
                instead.

            building_block_id:
                A unique id for each :class:`.Molecule` placed during
                the construction of the :class:`.ConstructedMolecule`.
                As a single :class:`.Molecule` can be placed multiple
                times during construction, the `building_block_id`
                allows the user to distinguish between each placement.
                Can be ``None``, if the atom was not part of a
                building block, but was added by the construction
                process instead.

        """

        self._atom = atom
        self._building_block_atom = building_block_atom
        self._building_block = building_block
        self._building_block_id = building_block_id

    def get_atom(self) -> Atom:
        """
        Get the atom about which information is held.

        Returns:

            The atom.

        """

        return self._atom

    def get_building_block_atom(self) -> typing.Optional[Atom]:
        """
        Get the original atom held by the building block.

        Returns:

            The building block atom or ``None`` if the atom was not
            originally found in a building block, but was added by the
            construction process instead.

        """

        return self._building_block_atom

    def get_building_block(self) -> typing.Optional[Molecule]:
        """
        Get the building block from which the atom originates.

        Returns:

            The building block or ``None`` if the atom was not
            originally found in a building block, but was added by the
            construction process instead.

        """

        return self._building_block

    def get_building_block_id(self) -> typing.Optional[int]:
        """
        Get the id of the atom's building block.

        A unique id for each :class:`.Molecule` placed during
        the construction of the :class:`.ConstructedMolecule`. As a
        single :class:`.Molecule` can be placed multiple times
        during construction, the building block id allows
        the user to distinguish between each placement.

        Returns:

            The unique building block id or ``None`` if the atom was
            not originally found in a building block, but was added by
            the construction process instead.

        """

        return self._building_block_id
