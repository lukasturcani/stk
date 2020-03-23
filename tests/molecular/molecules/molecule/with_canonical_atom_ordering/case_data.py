class CaseData:
    """
    A test case.

    Attributes
    ----------
    molecule : :class:`.Molecule`
        The molecule which should be canonically ordered.

    result : :class:`.Molecule`
        What the molecule should look like after canonical ordering.

    """

    def __init__(self, molecule, result):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        molecule : :class:`.Molecule`
            The molecule which should be canonically ordered.

        result : :class:`.Molecule`
            What the molecule should look like after canonical
            ordering.

        """

        self.molecule = molecule
        self.result = result
