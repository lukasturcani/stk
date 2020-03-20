class CaseData:
    """
    A test case.

    Attributes
    ----------
    constructed_molecule : :class:`.ConstructedMolecule`
        The molecule to test.

    num_new_atoms : :class:`int`
        The number of atoms added by the construction process, which
        do not belong to any building block.

    num_new_bonds : :class:`int`
        The number of bonds added by the construction process, which
        do not belong to any building block.

    num_building_blocks : :class:`dict`
        Maps each building block to the number of times it should
        be present in the constructed molecule.

    building_blocks : :class:`tuple` of :class:`.BuildingBlock`
        The building blocks of the constructed molecule.

    Notes
    -----
    You may wonder why only the number of new atoms and bonds is tested
    and not which specific atoms and bonds are new. This is because the
    atomic and bond ordering of the final constructed molecule is an
    implementation detail, and therefore it is extremely hard specify
    such a test.

    """

    def __init__(
        self,
        constructed_molecule,
        num_new_atoms,
        num_new_bonds,
        num_building_blocks,
        building_blocks,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        constructed_molecule : :class:`.ConstructedMolecule`
            The molecule to test.

        num_new_atoms : :class:`int`
            The number of atoms added by the construction process,
            which do not belong to any building block.

        num_new_bonds : :class:`int`
            The number of bonds added by the construction process,
            which do not belong to any building block.

        num_building_blocks : :class:`dict`
            Maps each building block to the number of times it should
            be present in the constructed molecule.

        building_blocks : :class:`tuple` of :class:`.BuildingBlock`
            The building blocks of the constructed molecule.

        """

        self.constructed_molecule = constructed_molecule
        self.num_new_atoms = num_new_atoms
        self.num_new_bonds = num_new_bonds
        self.num_building_blocks = num_building_blocks
        self.building_blocks = building_blocks
