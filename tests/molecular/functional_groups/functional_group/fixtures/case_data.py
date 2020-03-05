class CaseData:
    """
    A test case for :class:`.FunctionalGroup` instances.

    """

    def __init__(self, functional_group, atoms, placers):
        self.functional_group = functional_group
        self.atoms = atoms
        self.placers = placers


class GenericCaseData(CaseData):
    """
    A test case for :class:`.GenericFunctionalGroup` instances.

    """

    def __init__(self, functional_group, atoms, bonders, deleters):
        self.bonders = bonders
        self.deleters = deleters
        super().__init__(functional_group, atoms, bonders)
