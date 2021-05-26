class CaseData:
    """
    A test case.

    Attributes
    ----------
    fitness_calculator : :class:`.FitnessCalculator`
        The fitness calculator to test.

    molecule : :class:`.Molecule`
        The molecule whose fitness is calculated.

    fitness_value : :class:`object`
        The correct fitness value.

    """

    def __init__(self, fitness_calculator, molecule, fitness_value):
        """
        Initialize a :class:`.CaseData` instance.

        Parameters
        ----------
        fitness_calculator : :class:`.FitnessCalculator`
            The fitness calculator to test.

        molecule : :class:`.Molecule`
            The molecule whose fitness is calculated.

        fitness_value : :class:`object`
            The correct fitness value.

        """

        self.fitness_calculator = fitness_calculator
        self.molecule = molecule
        self.fitness_value = fitness_value
