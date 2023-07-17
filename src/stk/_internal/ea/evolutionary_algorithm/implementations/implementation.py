import itertools
import logging
import typing
from collections.abc import Callable, Iterable, Iterator

from stk._internal.ea.crossover.molecule_crosser import MoleculeCrosser
from stk._internal.ea.crossover.record import CrossoverRecord
from stk._internal.ea.fitness_calculators.fitness_calculator import (
    FitnessCalculator,
)
from stk._internal.ea.fitness_normalizers.fitness_normalizer import (
    FitnessNormalizer,
)
from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.mutation.mutator import MoleculeMutator
from stk._internal.ea.mutation.record import MutationRecord
from stk._internal.ea.selection.batch import Batch
from stk._internal.ea.selection.selectors.selector import Selector
from stk._internal.key_makers.molecule import MoleculeKeyMaker

from ...generation import FitnessValues, Generation

T = typing.TypeVar("T", bound=MoleculeRecord)
A = typing.TypeVar("A")
B = typing.TypeVar("B")
Map: typing.TypeAlias = Callable[[Callable[[A], B], Iterable[A]], B]


class Implementation(typing.Generic[T]):
    """
    An implementation of the default evolutionary algorithm.
    """

    def __init__(
        self,
        initial_population: Iterable[T],
        fitness_calculator: FitnessCalculator[T],
        mutator: MoleculeMutator[T],
        crosser: MoleculeCrosser[T],
        generation_selector: Selector[T],
        mutation_selector: Selector[T],
        crossover_selector: Selector[T],
        fitness_normalizer: FitnessNormalizer[T],
        key_maker: MoleculeKeyMaker,
        logger: logging.Logger,
    ) -> None:
        self._initial_population = tuple(initial_population)
        self._fitness_calculator = fitness_calculator
        self._mutator = mutator
        self._crosser = crosser
        self._generation_selector = generation_selector
        self._mutation_selector = mutation_selector
        self._crossover_selector = crossover_selector
        self._fitness_normalizer = fitness_normalizer
        self._key_maker = key_maker
        self._logger = logger

    def _get_generations(
        self,
        num_generations: int,
        map_: Map,
    ) -> Iterator[Generation[T]]:
        def get_mutation_record(batch: Batch[T]) -> MutationRecord[T] | None:
            return self._mutator.mutate(next(iter(batch)))

        def get_key(record: T) -> str:
            return self._key_maker.get_key(record.get_molecule())

        population, keys = dedupe(self._initial_population, get_key)

        self._logger.info("Calculating fitness values of initial population.")
        fitness_values = dict(
            zip(
                population,
                map_(self._fitness_calculator.get_fitness_value, population),
            )
        )
        normalized_fitness_values = self._fitness_normalizer.normalize(
            fitness_values=fitness_values,
        )
        yield Generation(
            fitness_values={
                record: FitnessValues(raw, normalized_fitness_values[record])
                for record, raw in fitness_values.items()
            },
            mutation_records=(),
            crossover_records=(),
        )

        for generation in range(1, num_generations):
            self._logger.info(f"Starting generation {generation}.")
            self._logger.info(f"Population size is {len(population)}.")

            self._logger.info("Doing crossovers.")
            crossover_records = tuple(
                self._get_crossover_records(normalized_fitness_values)
            )

            self._logger.info("Doing mutations.")
            mutation_records = tuple(
                record
                for record in map(
                    get_mutation_record,
                    self._mutation_selector.select(normalized_fitness_values),
                )
                if record is not None
            )

            self._logger.info("Calculating fitness values.")

            offspring, keys = dedupe(
                items=(
                    record.get_molecule_record()
                    for record in crossover_records
                ),
                get_key=get_key,
                seen=keys,
            )
            mutants, keys = dedupe(
                items=(
                    record.get_molecule_record() for record in mutation_records
                ),
                get_key=get_key,
                seen=keys,
            )
            fitness_values.update(
                zip(
                    itertools.chain(offspring, mutants),
                    map_(
                        self._fitness_calculator.get_fitness_value,
                        itertools.chain(offspring, mutants),
                    ),
                )
            )
            normalized_fitness_values = self._fitness_normalizer.normalize(
                fitness_values=fitness_values,
            )
            population, keys = dedupe(
                items=(
                    molecule_record
                    for molecule_record, in self._generation_selector.select(
                        population=normalized_fitness_values
                    )
                ),
                get_key=get_key,
            )
            fitness_values = {
                record: fitness_values[record] for record in population
            }
            normalized_fitness_values = {
                record: normalized_fitness_values[record]
                for record in population
            }
            yield Generation(
                fitness_values={
                    record: FitnessValues(
                        raw, normalized_fitness_values[record]
                    )
                    for record, raw in fitness_values.items()
                },
                mutation_records=mutation_records,
                crossover_records=crossover_records,
            )

    def _get_crossover_records(
        self,
        population: dict[T, float],
    ) -> Iterator[CrossoverRecord[T]]:
        for batch in self._crossover_selector.select(population):
            yield from self._crosser.cross(tuple(batch))


def dedupe(
    items: Iterable[A],
    get_key: Callable[[A], str],
    seen: set[str] | None = None,
) -> tuple[list[A], set[str]]:
    if seen is None:
        seen = set()
    unique = []
    for item in items:
        if (key := get_key(item)) not in seen:
            unique.append(item)
            seen.add(key)
    return unique, seen
