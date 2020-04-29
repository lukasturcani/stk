class CaseData:
    """
    A :class:`.Bond` test case.

    Attributes
    ----------
    bond : :class:`.Bond`
        The bond to test.

    atom1 : :class:`.Atom`
        The correct first atom of the bond.

    atom2 : :class:`.Atom`
        The correct second atom of the bond.

    order : :class:`int`
        The correct bond order of the bond.

    periodicity : :class:`tuple` of :class:`int`
        The correct periodicity of the bond.

    is_dative : :class:`bool`
        Correct assignment of whether the bond is dative.

    """

    def __init__(
        self,
        bond,
        atom1,
        atom2,
        order,
        periodicity,
        is_dative
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        bond : :class:`.Bond`
            The bond to test.

        atom1 : :class:`.Atom`
            The correct first atom of the bond.

        atom2 : :class:`.Atom`
            The correct second atom of the bond.

        order : :class:`int`
            The correct bond order of the bond.

        periodicity : :class:`tuple` of :class:`int`
            The correct periodicity of the bond.

        is_dative : :class:`bool`
            Correct assignment of whether the bond is dative.

        """

        self.bond = bond
        self.atom1 = atom1
        self.atom2 = atom2
        self.order = order
        self.periodicity = periodicity
        self.is_dative = is_dative
