"""
Evolutionary Algorithm
======================

"""


from ..fitness_normalizers import NullFitnessNormalizer
from .utilities import get_logger




class EvolutionaryAlgorithm:
    """
    An abstract base class for evolutionary algorithms.

    Notes
    -----


    """

    def __init__(
        self,
        initial_population,
        fitness_calculator,
        mutator,
        crosser,
        generation_selector,
        mutation_selector,
        crossover_selector,
        terminator,
        fitness_normalizer=NullFitnessNormalizer(),
        logger=get_logger(),
    ):
        """
        Initialize a :class:`EvolutionaryAlgorithm` instance.

        Parameters
        ----------

        """

        self._population = initial_population
        self._fitness_calculator = fitness_calculator
        self._fitness_normalizer = fitness_normalizer
        self._mutator = mutator
        self._crosser = crosser
        self._generation_selector = generation_selector
        self._mutation_selector = mutation_selector
        self._crossover_selector = crossover_selector
        self._terminator = terminator
        self._logger = logger

    def get_generations(self):
        """
        Yield the generations of the evolutionary algorithm.

        Yields
        ------
        :class:`.Generation`
            A generation.

        """

        self._logger.info(
            'Calculating the fitness of population members.'
        )

