"""
Evolutionary Algorithm
======================

"""

import logging

from stk.molecular import Inchi

from ..fitness_normalizers import NullFitnessNormalizer
from .implementations import Parallel, Serial

logger = logging.getLogger(__name__)


class EvolutionaryAlgorithm:
    """
    An abstract base class for evolutionary algorithms.

    Notes
    -----
    You might notice that the public methods of this abstract base
    class are implemented. This is purely for convenience, so that
    there is a default evolutionary algorithm implementation that
    users can use. However, feel free to override the default
    implementation when implementing subclasses.

    If you do want to use the default implementation, here is a
    summary of the roles of the different components:

    .. image:: https://i.imgur.com/hGXboaU.png

    Examples
    --------
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

    def __init__(
        self,
        initial_population,
        fitness_calculator,
        mutator,
        crosser,
        generation_selector,
        mutation_selector,
        crossover_selector,
        fitness_normalizer=NullFitnessNormalizer(),
        key_maker=Inchi(),
        num_processes=None,
    ):
        """
        Initialize a :class:`EvolutionaryAlgorithm` instance.

        Parameters
        ----------
        initial_population : :class:`tuple` of :class:`.MoleculeRecord`
            The initial population the EA should use.

        fitness_calculator : :class:`.FitnessCalculator`
            Calculates fitness values.

        mutator : :class:`.MoleculeMutator`
            Carries out mutation operations.

        crosser : :class:`.MoleculeCrosser`
            Carries out crossover operations.

        generation_selector : :class:`.Selector`
            Selects the next generation.

        mutation_selector : :class:`.Selector`
            Selects molecules for mutation.

        crossover_selector : :class:`.Selector`
            Selects molecules for crossover.

        fitness_normalizer : :class:`.FitnessNormalizer`
            Normalizes fitness values.

        key_maker : :class:`.MoleculeKeyMaker`, optional
            Used to detect duplicate molecules in the EA. If two
            molecules in a generation return the same key, one of them
            is removed.

        num_processes : :class:`int`, optional
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

    def get_generations(self, num_generations):
        """
        Yield the generations of the evolutionary algorithm.

        Parameters
        ----------
        num_generations : :class:`int`
            The number of generations which should be yielded.
            Note that the initial population counts as a generation.

        Yields
        ------
        :class:`.Generation`
            A generation.

        """

        yield from self._implementation.get_generations(
            num_generations=num_generations,
        )
