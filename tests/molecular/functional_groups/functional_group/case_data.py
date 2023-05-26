class CaseData:
    """
    A test case for :class:`.FunctionalGroup` instances.

    Attributes
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to test.

    atoms : :class:`tuple` of :class:`.Atom`
        The correct atoms held by the functional group.

    placers : :class:`tuple` of :class:`.Atom`
        The correct placer atoms held by the functional group.

    core_atoms : :class:`tuple` of :class:`.Atom`
        The correct core atoms held by the functional group.

    """

    def __init__(self, functional_group, atoms, placers, core_atoms):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        functional_group : :class:`.FunctionalGroup`
            The functional group to test.

        atoms : :class:`tuple` of :class:`.Atom`
            The correct atoms held by the functional group.

        placers : :class:`tuple` of :class:`.Atom`
            The correct placer atoms held by the functional group.

        core_atoms : :class:`tuple` of :class:`.Atom`
            The correct core atoms held by the functional group.

        """

        self.functional_group = functional_group
        self.atoms = atoms
        self.placers = placers
        self.core_atoms = core_atoms


class GenericCaseData(CaseData):
    """
    A test case for :class:`.GenericFunctionalGroup` instances.

    Attributes
    ----------
    functional_group : :class:`.FunctionalGroup`
        The functional group to test.

    atoms : :class:`tuple` of :class:`.Atom`
        The correct atoms held by the functional group.

    placers : :class:`tuple` of :class:`.Atom`
        The correct placer atoms held by the functional group.

    core_atoms : :class:`tuple` of :class:`.Atom`
        The correct core atoms held by the functional group.

    bonders : :class:`tuple` of :class:`.Atom`
        The correct bonder atoms held by the functional group.

    deleters : :class:`tuple` of :class:`.Atom`
        The correct deleter atoms held by the functional group.

    """

    def __init__(self, functional_group, atoms, bonders, deleters):
        """
        Initialize a :class:`.GenericCaseData` instance.

        Parameters
        ----------
        functional_group : :class:`.FunctionalGroup`
            The functional group to test.

        atoms : :class:`tuple` of :class:`.Atom`
            The correct atoms held by `functional_group`.

        bonders : :class:`tuple` of :class:`.Atom`
            The correct bonder atoms of `functional_group`.

        deleters : :class:`tuple` of :class:`.Atom`
            The correct deleter atoms of `functional_group`.

        """

        self.bonders = bonders
        self.deleters = deleters
        deleter_ids = set(atom.get_id() for atom in deleters)
        super().__init__(
            functional_group=functional_group,
            atoms=atoms,
            placers=bonders,
            core_atoms=tuple(
                atom for atom in atoms if atom.get_id() not in deleter_ids
            ),
        )
