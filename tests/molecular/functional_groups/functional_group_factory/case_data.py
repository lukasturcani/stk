class CaseData:
    """
    A test case.

    Attributes
    ----------
    factory : :class:`.FunctionalGroupFactory`
        The factory being tested.

    molecule : :class:`.Molecule`
        The molecule passed to
        :meth:`.FunctionalGroupFactory.get_functional_groups`.

    functional_groups : :class:`tuple` of :class:`.FunctionalGroup`
        The expected functional groups.

    """

    def __init__(self, factory, molecule, functional_groups):
        self.factory = factory
        self.molecule = molecule
        self.functional_groups = functional_groups
