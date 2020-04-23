class CaseData:
    """
    A test case.

    Attributes
    ----------
    key_maker : :class:`.MoleculeKeyMaker`
        The key maker to test.

    molecule : :class:`.Molecule`
        The molecule to test.

    key_name : :class:`str`
        The name of the key made by :attr:`.key_maker`.

    key : :class:`object`
        The correct key.

    """

    def __init__(
        self,
        key_maker,
        molecule,
        key_name,
        key,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        key_maker : :class:`.MoleculeKeyMaker`
            The key maker to test.

        molecule : :class:`.Molecule`
            The molecule to test.

        key_name : :class:`str`
            The name of the key made by `key_maker`.

        key : :class:`object`
            The correct key.

        """

        self.key_maker = key_maker
        self.molecule = molecule
        self.key_name = key_name
        self.key = key
