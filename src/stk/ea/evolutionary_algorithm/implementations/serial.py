from stk.utilities import dedupe
from ...generation import Generation

import itertools as it


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
        self._logger = logger

    def get_generations(self):

        def get_mutation_record(batch):
            return self._mutator.mutate(batch[0])

        logger = self._logger
        population = self._initial_population
        generation = 0

        self._logger.info(
            'Calculating fitness values of initial population.'
        )
        population = tuple(self._with_fitness_values(population))
        yield Generation(
            id=generation,
            molecule_records=population,
            mutation_records=(),
            crossover_records=(),
        )

        while not self._terminator.terminate(population):
            generation += 1

            logger.info(f'Starting generation {generation}.')
            logger.debug(f'Population size is {len(population)}.')

            logger.info('Doing crossovers.')
            crossover_records = tuple(map(
                self._crosser.cross,
                self._crossover_selector.select(population),
            ))

            logger.info('Doing mutations.')
            mutation_records = tuple(map(
                get_mutation_record,
                self._mutation_selector.select(population),
            ))

            logger.info('Calculating fitness values.')

            offspring = (
                record
                for crossover_record in crossover_records
                for record in crossover_record.get_molecule_records()
            )
            mutants = (
                record.get_molecule_record()
                for record in mutation_records
            )

            population = self._with_fitness_values(dedupe(
                iterable=it.chain(population, offspring, mutants),
                key=self._duplicate_key,
            ))

            population = tuple(
                molecule_record
                for molecule_record,
                in self._generation_selector.select(population)
            )

            yield Generation(
                id=generation,
                molecule_records=population,
                mutation_records=mutation_records,
                crossover_records=crossover_records,
            )

    def _with_fitness_values(self, population):
        fitness_calculator = self._fitness_calculator

        for record in population:
            if record.get_fitness_value() is None:
                fitness_value = fitness_calculator.get_fitness_value(
                    molecule=record.get_molecule(),
                )
                yield record.with_fitness_value(fitness_value)
            else:
                yield record
