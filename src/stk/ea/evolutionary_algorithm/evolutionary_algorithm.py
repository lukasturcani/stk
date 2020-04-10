"""
Evolutionary Algorithm
======================

"""


from ..fitness_normalizers import NullFitnessNormalizer
from .utilities import get_logger, get_inchi
from .implementations import Serial, Parallel


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
        duplicate_key=get_inchi,
        logger=get_logger(),
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

        duplicate_key : :class:`callable`, optional
            Takes a :class:`.MoleculeRecord` and returns some value.
            If two molecules return the same value, they are considered
            duplicates, and one of them will be removed from the
            population. By default
            :func:`~.evolutionary_algorithm.utilities.get_inchi` will
            be used.

        logger : :class:`logging.Logger`, optional
            The logger the EA should use.

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
                duplicate_key=duplicate_key,
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
                duplicate_key=duplicate_key,
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
