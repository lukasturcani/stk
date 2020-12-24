class CaseData:
    """
    A :class:`.Reaction` test case.

    Attributes
    ----------
    reaction : :class:`.Reaction`
        The reaction to test.

    new_atoms : :class:`tuple` of :class:`.NewAtom`
        The atoms, which should be added by :attr:`reaction`.

    new_bonds : :class:`tuple` of :class:`.Bond`
        The bonds, which should be added by :attr:`reaction`.

    deleted_atoms : :class:`tuple` of :class:`.Atom`
        The atoms, which should be deleted by :attr:`reaction`.

    deleted_bonds : :class:`tuple` of :class:`.Bond`
        The bonds, which should be deleted by :attr:`reaction`.

    """

    def __init__(
        self,
        reaction,
        new_atoms,
        new_bonds,
        deleted_atoms,
        deleted_bonds,
    ):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        reaction : :class:`.Reaction`
            The reaction to test.

        new_atoms : :class:`tuple` of :class:`.NewAtom`
            The atoms, which should be added by `reaction`.

        new_bonds : :class:`tuple` of :class:`.Bond`
            The bonds, which should be added by `reaction`.

        deleted_atoms : :class:`tuple` of :class:`.Atom`
            The atoms, which should be deleted by `reaction`.

        deleted_bonds : :class:`tuple` of :class:`.Bond`
            The bonds, which should be deleted by :attr:`reaction`.

        """

        self.reaction = reaction
        self.new_atoms = new_atoms
        self.new_bonds = new_bonds
        self.deleted_atoms = deleted_atoms
        self.deleted_bonds = deleted_bonds
