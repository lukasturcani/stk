"""
Parallel Evolutionary Algorithm
===============================

"""

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
        fitness_normalizer,
        key_maker,
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
            fitness_normalizer=fitness_normalizer,
            key_maker=key_maker,
            logger=logger,
        )
        self._num_processes = num_processes

    def get_generations(self, num_generations):
        with pathos.pools.ProcessPool(self._num_processes) as pool:
            yield from self._get_generations(num_generations, pool.map)
