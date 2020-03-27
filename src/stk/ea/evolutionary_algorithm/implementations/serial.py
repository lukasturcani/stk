from .implementation import Implementation


class Serial(Implementation):
    """
    A serial implementation of the default evolutionary algorithm.

    """

    def __init__(
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
        super().__init__(
            initial_population=initial_population,
            fitness_calculator=fitness_calculator,
            mutator=mutator,
            crosser=crosser,
            generation_selector=generation_selector,
            mutation_selector=mutation_selector,
            crossover_selector=crossover_selector,
            terminator=terminator,
            fitness_normalizer=fitness_normalizer,
            duplicate_key=duplicate_key,
            logger=logger,
            map=map,
        )

    def get_generations(self):
        yield from self._get_generations(
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
        )
