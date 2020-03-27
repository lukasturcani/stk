import pathos

from .implementation import Implementation


class Parallel(Implementation):
    """
    A parallel implementation of the default evolutionary algorithm.

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
        num_processes,
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
        )
        self._num_processes = num_processes

    def get_generations(self):
        with pathos.pools.ProcessPool(self._num_processes) as pool:
            yield from self._get_generations(pool.map)
