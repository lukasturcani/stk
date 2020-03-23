class CaseData:
    """
    A test case.

    Attributes
    ----------
    jsonizer : :class:`.MoleculeJsonizer` or \
            :class:`.ConstructedMoleculeJsonizer`
        The jsonizer to test.

    molecule : :class:`.Molecule`
        The molecule to JSONize.

    json : :class:`dict`
        The JSON of :attr:`.molecule`.

    """

    def __init__(self, jsonizer, molecule, json):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        jsonizer : :class:`.MoleculeJsonizer` or \
                :class:`.ConstructedMoleculeJsonizer`
            The jsonizer to test.

        molecule : :class:`.Molecule`
            The molecule to JSONize.

        json : :class:`dict`
            The correct JSON of `molecule`.

        """

        self.jsonizer = jsonizer
        self.molecule = molecule
        self.json = json
