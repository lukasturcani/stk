import itertools as it

from stk.utilities import dedupe
from ...generation import Generation


class Implementation:
    """
    An implementation of the default evolutionary algorithm.

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
        Initialize an :class:`.Implementation` instance.

        """

        self._initial_population = initial_population
        self._fitness_calculator = fitness_calculator
        self._mutator = mutator
        self._crosser = crosser
        self._generation_selector = generation_selector
        self._mutation_selector = mutation_selector
        self._crossover_selector = crossover_selector
        self._terminator = terminator
        self._fitness_normalizer = fitness_normalizer
        self._duplicate_key = duplicate_key
        self._logger = logger

    def _get_generations(self, map_):
        def get_mutation_record(batch):
            return self._mutator.mutate(batch[0])

        population = self._initial_population
        generation = 0

        self._logger.info(
            'Calculating fitness values of initial population.'
        )
        population = tuple(self._with_fitness_values(map_, population))
        yield Generation(
            id=generation,
            molecule_records=population,
            mutation_records=(),
            crossover_records=(),
        )

        while not self._terminator.terminate(population):
            generation += 1

            self._logger.info(f'Starting generation {generation}.')
            self._logger.debug(
                f'Population size is {len(population)}.'
            )

            self._logger.info('Doing crossovers.')
            crossover_records = tuple(
                self._get_crossover_records(population)
            )

            self._logger.info('Doing mutations.')
            mutation_records = tuple(map(
                get_mutation_record,
                self._mutation_selector.select(population),
            ))

            self._logger.info('Calculating fitness values.')

            offspring = (
                record.get_molecule_record()
                for record in crossover_records
            )
            mutants = (
                record.get_molecule_record()
                for record in mutation_records
            )

            population = self._with_fitness_values(
                map_=map_,
                population=dedupe(
                    iterable=it.chain(population, offspring, mutants),
                    key=self._duplicate_key,
                ),
            )
            population = self._fitness_normalizer.normalize(population)

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

    def _get_crossover_records(self, population):
        for batch in self._crossover_selector.select(population):
            yield from self._crosser.cross(batch)

    def _with_fitness_values(self, map_, population):
        no_fitness = (
            record for record in population
            if record.get_fitness_value(normalized=False) is None
        )
        molecules = (record.get_molecule() for record in no_fitness)
        fitness_values = map_(
            self._fitness_calculator.get_fitness_value,
            molecules,
        )
        yield from (
            record for record in population
            if record.get_fitness_value(normalized=False) is not None
        )
        for record, fitness_value in zip(no_fitness, fitness_values):
            yield record.with_fitness_value(
                fitness_value=fitness_value,
                normalized=False,
            )