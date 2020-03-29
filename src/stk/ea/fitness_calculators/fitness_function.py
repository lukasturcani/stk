from .fitness_calculator import FitnessCalculator


class FitnessFunction(FitnessCalculator):
    """
    Takes a function and uses it as a calculator.

    """

    def __init__(self, fitness_fn, use_cache=False):
        """
        Initialize a :class:`.FitnessFunction` instance.

        Parameters
        ----------
        fitness_fn : :class:`callable`
            Take a single parameter, the :class:`.Molecule` whose
            fitness needs to be calculated, and returns its
            fitness value.

        use_cache : :class:`bool`, optional
            If ``True`` a fitness calculation will not be performed on
            the same molecule twice, instead the previously returned
            value will be returned.

        """

        self._fitness_fn = fitness_fn
        super().__init__(use_cache=use_cache)

    def _get_fitness(self, mol):
        return self._fitness_fn(mol)
