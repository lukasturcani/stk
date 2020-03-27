"""
Evolutionary Algorithm
======================

"""

import itertools as it

from ..fitness_normalizers import NullFitnessNormalizer
from .utilities import get_logger, get_inchi
from .implementations import Serial, Parallel


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
        duplicate_key=get_inchi,
        logger=get_logger(),
        num_processes=1,
    ):
        """
        Initialize a :class:`EvolutionaryAlgorithm` instance.

        Parameters
        ----------

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
        self._logger = logger
        self._implementation = (
            Serial()
            if num_processes == 1
            else Parallel(num_processes)
        )

    def get_generations(self):
        """
        Yield the generations of the evolutionary algorithm.

        Yields
        ------
        :class:`.Generation`
            A generation.

        """

        yield from self._implementation.get_generations(self)

        logger = self._logger
        population = self._initial_population
        generation = 0
        while not self._terminator.terminate(population):
            generation += 1

            logger.info(f'Starting generation {generation}.')
            logger.debug(f'Population size is {len(population)}.')

            self._logger.info('Calculating fitness values.')
            population = tuple(self._with_fitness_values(population))

            logger.info('Doing crossovers.')
            crossover_parents = self._crossover_selector.select(
                population=population,
            )
            crossover_records = it.starmap(
                self._crosser.cross,
                crossover_parents,
            )
            offspring = (
                offspring
                for offspring_batch in crossover_records
                for offspring in offspring_batch
            )

    def _with_fitness_values(self, population):
        fitness_values = map(
            self._fitness_calculator,
            (record.get_molecule() for record in population),
        )
        for record, fitness_value in zip(population, fitness_values):
            yield record.with_fitness_value(fitness_value)
