import itertools as it

from stk.utilities import dedupe
from ...generation import Generation


class Implementation:
    def get_generations(self):
        def get_mutation_record(batch):
            return mutator.mutate(batch[0])

        population = initial_population
        generation = 0

        logger.info(
            'Calculating fitness values of initial population.'
        )
        population = tuple(
            with_fitness_values(fitness_calculator, population)
        )
        yield Generation(
            id=generation,
            molecule_records=population,
            mutation_records=(),
            crossover_records=(),
        )

        while not terminator.terminate(population):
            generation += 1

            logger.info(f'Starting generation {generation}.')
            logger.debug(f'Population size is {len(population)}.')

            logger.info('Doing crossovers.')
            crossover_records = tuple(map(
                crosser.cross,
                crossover_selector.select(population),
            ))

            logger.info('Doing mutations.')
            mutation_records = tuple(map(
                get_mutation_record,
                mutation_selector.select(population),
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

            population = with_fitness_values(
                fitness_calculator=fitness_calculator,
                population=dedupe(
                    iterable=it.chain(population, offspring, mutants),
                    key=duplicate_key,
                ),
            )
            population = with_normalized_fitness_values(
            )

            population = tuple(
                molecule_record
                for molecule_record,
                in generation_selector.select(population)
            )

            yield Generation(
                id=generation,
                molecule_records=population,
                mutation_records=mutation_records,
                crossover_records=crossover_records,
            )


    def with_fitness_values(fitness_calculator, population):
        for record in population:
            if record.get_fitness_value() is None:
                fitness_value = fitness_calculator.get_fitness_value(
                    molecule=record.get_molecule(),
                )
                yield record.with_fitness_value(fitness_value)
            else:
                yield record
