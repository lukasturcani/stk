class CaseData:
    """
    A test case for :class:`.FunctionalGroup` instances.

    """

    def __init__(self, functional_group, atoms, placers, core_atoms):
        self.functional_group = functional_group
        self.atoms = atoms
        self.placers = placers
        self.core_atoms = core_atoms


class GenericCaseData(CaseData):
    """
    A test case for :class:`.GenericFunctionalGroup` instances.

    """

    def __init__(self, functional_group, atoms, bonders, deleters):
        self.bonders = bonders
        self.deleters = deleters
        deleter_set = set(deleters)
        super().__init__(
            functional_group=functional_group,
            atoms=atoms,
            placers=bonders,
            core_atoms=tuple(
                atom for atom in atoms if atom not in deleter_set
            ),
        )
