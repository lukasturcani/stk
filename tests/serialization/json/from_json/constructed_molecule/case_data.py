class CaseData:
    """
    A test case.

    Attributes
    ----------
    dejsonizer : :class:`.ConstructedMoleculeDejsonizer`
        The dejsonizer to test.

    molecule_json : :class:`dict`
        A JSON of the molecular information of the constructed
        molecule.

    constructed_molecule_json : :class:`dict`
        A JSON of the constructed molecule information of the
        constructed molecule.

    position_matrix : :class:`numpy.ndarray`
        The position matrix of the constructed molecule.

    building_blocks : :class:`tuple` of :class:`.Molecule`
        The building blocks of the constructed molecule.

    molecule : :class:`.ConstructedMolecule`
        The correct dejsonized molecule.

    """

    def __init__(
        self,
        dejsonizer,
        molecule_json,
        constructed_molecule_json,
        position_matrix,
        building_blocks,
        molecule,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        dejsonizer : :class:`.ConstructedMoleculeDejsonizer`
            The dejsonizer to test.

        molecule_json : :class:`dict`
            A JSON of the molecular information of the constructed
            molecule.

        constructed_molecule_json : :class:`dict`
            A JSON of the constructed molecule information of the
            constructed molecule.

        position_matrix : :class:`numpy.ndarray`
            The position matrix of the constructed molecule.

        building_blocks : :class:`tuple` of :class:`.Molecule`
            The building blocks of the constructed molecule.

        molecule : :class:`.ConstructedMolecule`
            The correct dejsonized molecule.

        """

        self.dejsonizer = dejsonizer
        self.molecule_json = molecule_json
        self.constructed_molecule_json = constructed_molecule_json
        self.position_matrix = position_matrix
        self.building_blocks = building_blocks
        self.molecule = molecule
