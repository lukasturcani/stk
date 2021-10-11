"""
Atom Info
=========

"""


class AtomInfo:
    """
    Holds additional info about :class:`.ConstructedMolecule` atoms.

    """

    def __init__(
        self,
        atom,
        building_block_atom,
        building_block,
        building_block_id,
    ):
        """
        Initialize an :class:`.AtomInfo` instance.

        Parameters
        ----------
        atom : :class:`.Atom`
            The atom about which information is held.

        building_block_atom : :class:`.Atom`
            The building block atom from which this atom originates.
            Can be ``None``, if the atom was not part of the building
            block, but was added by the construction process instead.

        building_block : :class:`.Molecule` or :class:`NoneType`
            The building block from which the atom originates.
            Can be ``None``, if the atom was not part of a building
            block, but was added by the construction process instead.

        building_block_id : :class:`int` or :class:`NoneType`
            A unique id for each :class:`.Molecule` placed during
            the construction of the :class:`.ConstructedMolecule`. As a
            single :class:`.Molecule` can be placed multiple times
            during construction, the `building_block_id` allows
            the user to distinguish between each placement. Can be
            ``None``, if the atom was not part of a building block, but
            was added by the construction process instead.

        """

        self._atom = atom
        self._building_block_atom = building_block_atom
        self._building_block = building_block
        self._building_block_id = building_block_id

    def get_atom(self):
        """
        Get the atom about which information is held.

        Returns
        -------
        :class:`.Atom`
            The atom.

        """

        return self._atom

    def get_building_block_atom(self):
        """
        Get the original atom held by the building block.

        Returns
        -------
        :class:`.Atom`
            The building block atom.

        None : :class:`NoneType`
            If the atom was not originally found in a building block,
            but was added by the construction process instead.

        """

        return self._building_block_atom

    def get_building_block(self):
        """
        Get the building block from which the atom originates.

        Returns
        -------
        :class:`.Molecule`
            The building block.

        None : :class:`NoneType`
            If the atom was not originally found in a building block,
            but was added by the construction process instead.

        """

        return self._building_block

    def get_building_block_id(self):
        """
        Get the id of the atom's building block.

        A unique id for each :class:`.Molecule` placed during
        the construction of the :class:`.ConstructedMolecule`. As a
        single :class:`.Molecule` can be placed multiple times
        during construction, the building block id allows
        the user to distinguish between each placement.

        Returns
        -------
        :class:`int`
            The id.

        None : :class:`NoneType`
            If the atom was not originally found in a building block,
            but was added by the construction process instead.

        """

        return self._building_block_id
