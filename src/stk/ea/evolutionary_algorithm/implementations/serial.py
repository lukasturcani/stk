import itertools as it


class Serial:
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
