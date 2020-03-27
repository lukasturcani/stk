from .utilities import get_generations


class Serial:
    """
    A serial implementation of the default evolutionary algorithm.

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
        fitness_normalizer,
        duplicate_key,
        logger,
    ):
        """
        Initialize a :class:`.Serial` instance.

        """

        self._initial_population = initial_population
        self._fitness_calculator = fitness_calculator
        self._fitness_normalizer = fitness_normalizer
        self._mutator = mutator
        self._crosser = crosser
        self._generation_selector = generation_selector
        self._mutation_selector = mutation_selector
        self._crossover_selector = crossover_selector
        self._terminator = terminator
        self._duplicate_key = duplicate_key
        self._logger = logger

    def get_generations(self):
        yield from get_generations(
            initial_population=self._initial_population,
            fitness_calculator=self._fitness_calculator,
            mutator=self._mutator,
            crosser=self._crosser,
            generation_selector=self._generation_selector,
            mutation_selector=self._mutation_selector,
            crossover_selector=self._crossover_selector,
            terminator=self._terminator,
            fitness_normalizer=self._fitness_normalizer,
            duplicate_key=self._duplicate_key,
            logger=self._logger,
            map=map,
        )
