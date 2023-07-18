import logging
import typing
from collections.abc import Iterable, Iterator

from stk._internal.ea.crossover.molecule_crosser import MoleculeCrosser
from stk._internal.ea.fitness_calculators.fitness_calculator import (
    FitnessCalculator,
)
from stk._internal.ea.fitness_normalizers.fitness_normalizer import (
    FitnessNormalizer,
)
from stk._internal.ea.generation import Generation
from stk._internal.ea.molecule_record import MoleculeRecord
from stk._internal.ea.mutation.mutator import MoleculeMutator
from stk._internal.ea.selection.selectors.selector import Selector
from stk._internal.key_makers.inchi import Inchi
from stk._internal.key_makers.molecule import MoleculeKeyMaker

from ..fitness_normalizers.null import NullFitnessNormalizer
from .implementations.parallel import Parallel
from .implementations.serial import Serial

logger = logging.getLogger(__name__)

T = typing.TypeVar("T", bound=MoleculeRecord)


class EvolutionaryAlgorithm(typing.Generic[T]):
    """
    An abstract base class for evolutionary algorithms.

    Notes:

        You might notice that the public methods of this abstract base
        class are implemented. This is purely for convenience, so that
        there is a default evolutionary algorithm implementation that
        users can use. However, feel free to override the default
        implementation when implementing subclasses.

        If you do want to use the default implementation, here is a
        summary of the roles of the different components:

        .. image:: https://i.imgur.com/hGXboaU.png

    Examples:

        *Subclass Implementation*

        The source code of this class can work as a good example. There
        is only one method that a subclass of
        :class:`.EvolutionaryAlgorithm` needs to implement,
        :meth:`.get_generations`, which yields :class:`.Generation`
        instances. These correspond to the generations of your
        evolutionary algorithm implementation.

        *Usage*

        There are a couple of tutorials on how to use the
        :class:`.EvolutionaryAlgorithm`, which can be found in the sidebar.

    """

    _implementation: Serial | Parallel

    def __init__(
        self,
        initial_population: Iterable[T],
        fitness_calculator: FitnessCalculator[T],
        mutator: MoleculeMutator[T],
        crosser: MoleculeCrosser[T],
        generation_selector: Selector[T],
        mutation_selector: Selector[T],
        crossover_selector: Selector[T],
        fitness_normalizer: FitnessNormalizer[T] = NullFitnessNormalizer(),
        key_maker: MoleculeKeyMaker = Inchi(),
        num_processes: int | None = None,
    ) -> None:
        """
        Parameters:

            initial_population (list[T]):
                The initial population the EA should use.

            fitness_calculator:
                Calculates fitness values.

            mutator:
                Carries out mutation operations.

            crosser:
                Carries out crossover operations.

            generation_selector:
                Selects the next generation.

            mutation_selector:
                Selects molecules for mutation.

            crossover_selector:
                Selects molecules for crossover.

            fitness_normalizer:
                Normalizes fitness values.

            key_maker:
                Used to detect duplicate molecules in the EA. If two
                molecules in a generation return the same key, one of them
                is removed.

            num_processes:
                The number of parallel processes the EA should create.
                If ``None``, all available cores will be used.
        """
        if num_processes == 1:
            self._implementation = Serial(
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

        else:
            self._implementation = Parallel(
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
                num_processes=num_processes,
            )

    def get_generations(self, num_generations: int) -> Iterator[Generation[T]]:
        """
        Yield the generations of the evolutionary algorithm.

        Parameters:
            num_generations:
                The number of generations which should be yielded.
                Note that the initial population counts as a generation.
        Yields:
            A generation.
        """
        yield from self._implementation.get_generations(
            num_generations=num_generations,
        )
